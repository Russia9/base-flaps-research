/*
 * Parametric fuselage with aft-base arc fins
 *
 * Coordinate system (after assembly):
 *   X — axial, nose at origin, positive toward base and beyond
 *   Y — lateral, first fin centered on +Y axis
 *   Z — lateral
 *
 * All dimensions are in units of D (fuselage diameter).
 * Set D = 1 for normalized geometry; override for dimensional output.
 *
 * Command-line override example:
 *   openscad -o out.stl -D "N=2; xi=45; LD=1.0;" model.scad
 */

// ── Parameters (overridable from CLI) ────────────────────────────────────────

D  = 80.0;       // fuselage diameter
N  = 2;         // number of fins [1..4]
xi = 45;        // fin arc angle [degrees]
LD = 1.0;       // fin length / D
TD = 0.02;      // fin thickness / D (default 0.02)

// ── Derived constants ─────────────────────────────────────────────────────────

R          = D / 2;
t          = TD * D;
fin_len    = LD * D;
ogive_rho  = 8.5 * D;
total_len  = 10.0 * D;
ogive_len  = sqrt(ogive_rho * ogive_rho - (ogive_rho - R) * (ogive_rho - R));
cyl_len    = total_len - ogive_len;

echo(str("Ogive length     = ", ogive_len, " D"));
echo(str("Cylinder length  = ", cyl_len,   " D"));

// ── Resolution ────────────────────────────────────────────────────────────────

FN_BODY = 360;  // circumferential facets on fuselage
FN_NOSE = 360;  // profile points along ogive
FN_FIN  = 360;  // circumferential facets on fin arcs

// ── Ogive profile ─────────────────────────────────────────────────────────────

// Radius at axial distance x from nose tip along ogive (0 ≤ x ≤ ogive_len)
function ogive_r(x) =
    sqrt(ogive_rho * ogive_rho - (ogive_len - x) * (ogive_len - x))
    - (ogive_rho - R);

// 2D profile for rotate_extrude: [r, z] pairs
// Traces:  nose tip → ogive curve → cylinder wall → base edge → axis → tip
function body_profile() = concat(
    [[0, 0]],
    [for (i = [1 : FN_NOSE])
        let(x = i * ogive_len / FN_NOSE)
        [ogive_r(x), x]
    ],
    [[R, total_len],
     [0, total_len]]
);

// ── Fuselage ──────────────────────────────────────────────────────────────────

// Z-axis-aligned body; tip at origin, base at Z = total_len
module fuselage() {
    rotate_extrude($fn = FN_BODY)
    polygon(body_profile());
}

// ── Fins ──────────────────────────────────────────────────────────────────────

// Arc cross-section (XY plane), centred on +Y axis.
// Outer arc at radius R (flush with fuselage), inner arc at R - t.
module fin_xsec() {
    a_start = 90 - xi / 2;
    a_end   = 90 + xi / 2;
    n_arc   = max(4, round(xi));   // ~1 point per degree
    r_cut   = R * 1.05;            // sector radius, must exceed R

    intersection() {
        difference() {
            circle(r = R,     $fn = FN_FIN);
            circle(r = R - t, $fn = FN_FIN);
        }
        // Pie-slice sector selects only the desired arc
        polygon(concat(
            [[0, 0]],
            [for (i = [0 : n_arc])
                let(a = a_start + i * xi / n_arc)
                [r_cut * cos(a), r_cut * sin(a)]
            ]
        ));
    }
}

// Single fin: cross-section extruded axially from base
module fin_one() {
    translate([0, 0, total_len])
    linear_extrude(height = fin_len)
    fin_xsec();
}

// All N fins equally spaced; first fin centred on +Y, rest at +k*360/N
module fins() {
    for (k = [0 : N - 1])
        rotate([0, 0, k * 360 / N])
        fin_one();
}

// ── Assembly ──────────────────────────────────────────────────────────────────

// rotate([0, 90, 0]) maps Z → X, leaving Y unchanged,
// so the Z-aligned body becomes X-aligned with the first fin at +Y.
module assembly() {
    rotate([0, 90, 0]) union() {
        fuselage();
        fins();
    }
}

// ── Named boundary regions for snappyHexMesh / OpenFOAM ──────────────────────
//
// Each module below corresponds to one named solid in the exported STL.
// The Makefile renders each module separately and fixes the solid name via sed:
//   openscad -D "EXPORT=\"body\"" -o body.stl model.scad
//   sed -i '' 's/OpenSCAD_Model/body/g' body.stl
//
// snappyHexMeshDict geometry entry:
//   body.stl { type triSurfaceMesh; name body; }

EXPORT = "";   // set from CLI: -D "EXPORT=\"body\""

module body() { assembly(); }

if      (EXPORT == "body") body();
else if (EXPORT == "")     body();   // default interactive render
