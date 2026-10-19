# CFD server setup

How to rebuild the CFD server from a fresh OS, as recorded from the machine on
2026-10-18 (Debian 13 "trixie", kernel 6.12, 24 cores, 62 GB RAM, 467 GB disk).
The versions below are the ones every result in `results/` was produced with.

| Component | Version | Source |
|---|---|---|
| OpenFOAM | v2512 (`openfoam2512-default` 2512.0-2) | dl.openfoam.com apt repo |
| HiSA | git `271b449` ("Fix redistributeHelper for agard445") | https://gitlab.com/hisa/hisa.git, no local changes |
| Open MPI | 5.0.7 | Debian (pulled in by OpenFOAM) |
| Gmsh | 4.15.2 (pip), Python 3.13.5 | venv `~/venvs/gmsh` |
| OpenSCAD | `openscad-nightly` 2026.09.09 (git 33e8fdd8) | OBS home:t-paul repo |
| uv | 0.12.16 | astral.sh installer, `~/.local/bin` |

## 1. Packages

```bash
sudo apt update
sudo apt install -y build-essential cmake gfortran git rsync screen curl python3-venv
sudo apt install -y libxft2   # the pip gmsh wheel links libXft even when headless
sudo apt install -y libjemalloc2   # run-simulation.sh preloads it for HiSA

# OpenFOAM v2512 (ESI): adds /etc/apt/sources.list.d/openfoam.list
curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt install -y openfoam2512-default openfoam2512-dev
```

`openfoam2512-dev` provides the headers and wmake that HiSA needs. It pulls in
`libopenmpi-dev`, `libscotch-dev` and `libptscotch-dev`.

OpenSCAD is only needed for the snappy cases (`openfoam/arc*`). The Gmsh path
doesn't use it.

```bash
# OBS repository home:t-paul, Debian_13
curl -fsSL https://download.opensuse.org/repositories/home:/t-paul/Debian_13/Release.key \
  | sudo tee /etc/apt/trusted.gpg.d/obs-openscad-nightly.asc >/dev/null
echo "deb https://download.opensuse.org/repositories/home:/t-paul/Debian_13/ ./" \
  | sudo tee /etc/apt/sources.list.d/openscad-nightly.list
sudo apt update && sudo apt install -y openscad-nightly
```

`rebuild-mesh.sh` prefers `openscad-nightly` and falls back to `openscad`.
The Manifold backend and the cylinder-rung behaviour noted in `CLAUDE.md` are
the 2026.09 nightly's.

## 2. HiSA

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
git clone https://gitlab.com/hisa/hisa.git ~/hisa
cd ~/hisa && git checkout 271b449
./Allwmake -j 24
```

This installs into `$FOAM_USER_APPBIN` and `$FOAM_USER_LIBBIN`
(`~/OpenFOAM/debian-v2512/platforms/linux64GccDPInt32Opt/`). Check with
`which hisa`.

## 3. Python tools

```bash
python3 -m venv ~/venvs/gmsh
~/venvs/gmsh/bin/pip install gmsh==4.15.2
curl -LsSf https://astral.sh/uv/install.sh | sh     # optional; the repo scripts use GMSH_PYTHON on the server
```

All other scripts are stdlib-only and use the system `python3`. The server has
no numpy, and none is needed.

## 4. Repository

```bash
git clone https://github.com/Russia9/base-flaps-research.git ~/base-flaps-research
```

The server only reads from git (HTTPS, no key). It syncs with
`git fetch && git reset --hard origin/master` and never commits; see
`EXPERIMENTS.md`.

## 5. Running a case

See README §5. Non-interactive shells have no OpenFOAM environment, so a run
is a small script started with `nohup`:

```bash
cat > ~/run_case.sh <<'SCRIPT'
#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc   # before any set -u
cd ~/base-flaps-research
C=openfoam/<case>
GMSH_PYTHON=~/venvs/gmsh/bin/python MAX_CELLS=13000000 ./rebuild-mesh.sh $C || exit 1
./run-simulation.sh $C
python3 scripts/wall_breakdown.py $C
python3 scripts/te_map.py $C openfoam/gm_arc_M1.6 > results/$(basename $C)/te_map.txt
SCRIPT
chmod +x ~/run_case.sh && nohup ~/run_case.sh > ~/run_case.log 2>&1 < /dev/null &
```

`run-simulation.sh` pins one HiSA rank per physical core (`--map-by core
--bind-to core`; the 9900X has 12 cores, SMT siblings stay idle) and preloads
jemalloc with retained, never-decayed memory. Set `BIND=0` or `JEMALLOC=` to
turn either off. Run one case at a time: two pinned runs would share cores 0–11.

Limits measured on this machine:
- **Memory:** HiSA on 12 ranks uses about 4.3–4.5 GB per million cells at
  peak. 12.0 M cells peaked at 50 GB; 13.07 M was OOM-killed.
- **Disk:** a solved 12 M case keeps a ~1.9 GB snapshot per write interval in
  `processor*/`. After the results are committed, delete every
  `processor*/<time>` except `0` and the latest.

## Not needed

- `microfluidica/openfoam:13` (Docker) was used once, to test OpenFOAM
  Foundation v13's stitchMesh. The ESI stitchMesh is the one in use; see the
  stitch notes in `CLAUDE.md`.
- `~/.ssh/id_ed25519` is not used by this project. Git works over HTTPS.
