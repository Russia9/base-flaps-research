#!/usr/bin/env python3
"""Break the wall forces of a solved case down by body region and fin chord.

Usage:
    python3 scripts/wall_breakdown.py <case-dir> [--time T]

Samples p and wallShearStress on the wall patches with a one-off
`postProcess` (OpenFOAM env must be sourced), then integrates per face:
- axial force split into nose / cylinder / fins / base, so the report's
  C_Ab (base disk) and C_Af = C_A - C_Ab can be compared directly;
- fin rolling moment Mx by 10 mm chord bin and by surface (concave or
  convex side, with the LE/TE bevels split out as *_LE / *_TE);
- base-disk Cp by radius.

Writes results/<case-name>/wall_breakdown.csv. Same normalisation as
post_process.py: q_inf*S and q_inf*S*D, moments about the nose tip.
"""

from __future__ import annotations

import argparse
import math
import re
import subprocess
import tempfile
from pathlib import Path

from post_process import D_REF, S_REF, parse_freestream

SAMPLE_DICT = """FoamFile { version 2.0; format ascii; class dictionary; object wallSample; }
functions
{
    wallSample
    {
        type            surfaces;
        libs            (sampling);
        surfaceFormat   foam;
        fields          (p wallShearStress);
        interpolationScheme cell;
        surfaces
        {
            fuselage    { type patch; patches (fuselage);    interpolate false; }
            stabilizers { type patch; patches (stabilizers); interpolate false; }
        }
    }
}
"""

CHORD_BIN = 0.010  # m
BASE_BIN = 0.005  # m


def list_body(path: Path) -> str:
    text = path.read_text()
    start = re.search(r"^\d+\s*$", text, re.MULTILINE).end()
    return text[text.index("(", start) + 1:text.rindex(")")]


def read_scalars(path: Path) -> list[float]:
    return [float(v) for v in list_body(path).split()]


def read_vectors(path: Path) -> list[tuple[float, ...]]:
    return [tuple(map(float, v.split())) for v in re.findall(r"\(([^()]+)\)", list_body(path))]


def read_faces(path: Path) -> list[list[int]]:
    return [list(map(int, v.split())) for v in re.findall(r"\d+\(([^()]+)\)", list_body(path))]


def patch_faces(sample: Path, patch: str):
    """Yield (centre, area vector, p, wall shear) per face; area vector points into the wall."""
    d = sample / patch
    pts = read_vectors(d / "points")
    faces = read_faces(d / "faces")
    p = read_scalars(d / "scalarField" / "p")
    tau = read_vectors(d / "vectorField" / "wallShearStress")
    if not len(faces) == len(p) == len(tau):
        raise ValueError(f"{patch}: face/field size mismatch")
    for f, pv, tv in zip(faces, p, tau):
        sx = sy = sz = cx = cy = cz = 0.0
        for a, b in zip(f, f[1:] + f[:1]):
            (x1, y1, z1), (x2, y2, z2) = pts[a], pts[b]
            sx += y1 * z2 - z1 * y2
            sy += z1 * x2 - x1 * z2
            sz += x1 * y2 - y1 * x2
            cx += x1
            cy += y1
            cz += z1
        n = len(f)
        yield (cx / n, cy / n, cz / n), (0.5 * sx, 0.5 * sy, 0.5 * sz), pv, tv


def read_case_scalar(case: Path, name: str) -> float:
    text = (case / "constant" / "caseProperties").read_text()
    return float(re.search(rf"^\s*{name}\s+([-+0-9.eE]+)\s*;", text, re.MULTILINE).group(1))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", type=Path)
    ap.add_argument("--time", default="latestTime", help="time directory to sample (default: latest)")
    args = ap.parse_args()

    case = args.case.resolve()
    fs = parse_freestream(case)
    p_inf, q = fs["pInf"], fs["qInf"]
    diameter = read_case_scalar(case, "D") / 1000.0
    nose_len, body_len = 2.0 * diameter, 10.0 * diameter

    with tempfile.TemporaryDirectory() as tmp:
        dict_path = Path(tmp) / "wallSample"
        dict_path.write_text(SAMPLE_DICT)
        time_opt = ["-latestTime"] if args.time == "latestTime" else ["-time", args.time]
        subprocess.run(["postProcess", "-case", str(case), "-dict", str(dict_path),
                        "-fields", "(p wallShearStress)", *time_opt],
                       check=True, stdout=subprocess.DEVNULL)
    sample_root = case / "postProcessing" / "wallSample"
    sample = max(sample_root.iterdir(), key=lambda t: float(t.name))

    # key -> [Cx_p, Cx_v, Mx_p, Mx_v, area]; accumulated dimensional, normalised at the end
    rows: dict[tuple[str, str, str], list[float]] = {}

    def add(key, force_p, force_v, centre, area):
        r = rows.setdefault(key, [0.0] * 5)
        _, y, z = centre
        r[0] += force_p[0]
        r[1] += force_v[0]
        r[2] += y * force_p[2] - z * force_p[1]
        r[3] += y * force_v[2] - z * force_v[1]
        r[4] += area

    patches = ["fuselage"] + (["stabilizers"] if (sample / "stabilizers").is_dir() else [])
    fin_lateral_sum: dict[int, float] = {}
    fin_faces = []
    for patch in patches:
        for c, s, pv, tv in patch_faces(sample, patch):
            area = math.sqrt(s[0] ** 2 + s[1] ** 2 + s[2] ** 2)
            f_p = tuple((pv - p_inf) * si for si in s)
            # wallShearStress is the stress the wall exerts on the fluid; flip it for force on the body
            f_v = tuple(-ti * area for ti in tv)
            if patch == "stabilizers":
                # Fins root on the +Y, +Z, -Y, -Z axes (index 0..3) and bow sideways from them.
                phi = math.atan2(c[2], c[1])
                k = round(phi / (math.pi / 2)) % 4
                ay, az = math.cos(k * math.pi / 2), math.sin(k * math.pi / 2)
                fin_lateral_sum[k] = fin_lateral_sum.get(k, 0.0) + (-az * c[1] + ay * c[2]) * area
                fin_faces.append((k, c, s, area, f_p, f_v))
                add(("region", "fins", ""), f_p, f_v, c, area)
            elif c[0] > body_len - 5e-4 and abs(s[0]) > 0.99 * area:
                add(("region", "base", ""), f_p, f_v, c, area)
                r = math.hypot(c[1], c[2])
                add(("base_r", f"{int(r / BASE_BIN) * BASE_BIN * 1000:.0f}", ""), f_p, f_v, c, area)
            elif c[0] < nose_len:
                add(("region", "nose", ""), f_p, f_v, c, area)
            else:
                add(("region", "cylinder", ""), f_p, f_v, c, area)

    fin_le = min(c[0] for _, c, _, _, _, _ in fin_faces) if fin_faces else 0.0
    for k, c, s, area, f_p, f_v in fin_faces:
        ay, az = math.cos(k * math.pi / 2), math.sin(k * math.pi / 2)
        bow = math.copysign(1.0, fin_lateral_sum[k])
        # Fluid-side normal is -s. The concave side faces away from the bow.
        n_lat = -(-az * s[1] + ay * s[2]) / area * bow
        side = "convex" if n_lat > 0 else "concave"
        # The 45 deg edge wedges tilt the normal sin(22.5 deg) = 0.38 off the flat surface.
        if abs(s[0]) > 0.2 * area:
            side += "_LE" if -s[0] < 0 else "_TE"
        chord = f"{int((c[0] - fin_le) / CHORD_BIN) * CHORD_BIN * 1000:.0f}"
        add(("fin_chord", chord, side), f_p, f_v, c, area)

    q_s, q_s_d = q * S_REF, q * S_REF * D_REF
    out_dir = Path("results") / case.name
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "wall_breakdown.csv"
    with out.open("w") as fh:
        fh.write(f"# sampled at time {sample.name}\n")
        fh.write("group,key,side,Cx_p,Cx_v,Mx_p,Mx_v,area_over_S,Cp_mean\n")
        for (group, key, side), (fp, fv, mp, mv, a) in sorted(
                rows.items(), key=lambda kv: (kv[0][0], float(kv[0][1]) if kv[0][1][0].isdigit() else 0, kv[0][1], kv[0][2])):
            # base-disk normal points upstream into the body, so Cp = -F_x / (q A)
            cp = -fp / a / q if group == "base_r" else float("nan")
            fh.write(f"{group},{key},{side},{fp / q_s:.6f},{fv / q_s:.6f},"
                     f"{mp / q_s_d:.6f},{mv / q_s_d:.6f},{a / S_REF:.6f},{cp:.4f}\n")

    reg = {k[1]: v for k, v in rows.items() if k[0] == "region"}
    base = reg.get("base", [0.0] * 5)
    total_x = sum(v[0] + v[1] for v in reg.values())
    print(f"case  : {case.name} (time {sample.name})")
    for name in ("nose", "cylinder", "fins", "base"):
        if name in reg:
            v = reg[name]
            print(f"{name:9s} Cx_p {v[0] / q_s:+.4f}  Cx_v {v[1] / q_s:+.4f}  "
                  f"Mx {(v[2] + v[3]) / q_s_d:+.5f}")
    print(f"C_A {total_x / q_s:.4f}   C_Ab {(base[0] + base[1]) / q_s:.4f}   "
          f"C_Af {(total_x - base[0] - base[1]) / q_s:.4f}")
    print(f"output: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
