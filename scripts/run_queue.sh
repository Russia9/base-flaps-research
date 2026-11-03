#!/usr/bin/env bash
# Mesh, solve and post-process several cases one after another, unattended.
#
# Usage (from the repo root, typically under nohup or screen on the CFD server):
#   nohup scripts/run_queue.sh openfoam/<caseA> openfoam/<caseB> ... > ~/queue.log 2>&1 < /dev/null &
#
# Per case: rebuild-mesh.sh -> run-simulation.sh -> wall_breakdown.py -> te_map.py.
# A failing step is recorded and the queue moves on to the next case. The
# outcome of every step goes to queue.log in the repo root (a case's own
# log.* files are wiped by rebuild-mesh.sh). Environment: NP, GMSH_PYTHON,
# MAX_CELLS (default 13300000, the HiSA memory ceiling on 12 ranks),
# OPENFOAM_BASHRC (sourced when checkMesh is not on PATH).

# The OpenFOAM bashrc is not clean under set -u: source it first.
if ! command -v checkMesh >/dev/null 2>&1; then
    # shellcheck disable=SC1090
    source "${OPENFOAM_BASHRC:-/usr/lib/openfoam/openfoam2512/etc/bashrc}"
fi
set -uo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT" || exit 2
export MAX_CELLS=${MAX_CELLS:-13300000}
LOG="$ROOT/queue.log"
REF=${TE_MAP_REF:-openfoam/gm_arc_M1.6}     # second case for te_map.py's side-by-side map

[ $# -ge 1 ] || { sed -n '2,13p' "$0" >&2; exit 2; }

stamp() { date -u +%Y-%m-%dT%H:%M:%SZ; }
note() { echo "$(stamp) $*" | tee -a "$LOG"; }

declare -A STATUS
for case in "$@"; do
    name=$(basename "$case")
    if [ ! -d "$case" ]; then
        note "$name: missing case directory"
        STATUS[$name]="missing"
        continue
    fi
    note "$name: mesh"
    if ! ./rebuild-mesh.sh "$case" > "$ROOT/queue.$name.mesh.log" 2>&1; then
        note "$name: MESH FAILED (queue.$name.mesh.log, $case/log.checkMesh)"
        STATUS[$name]="mesh failed"
        continue
    fi
    cells=$(awk '$1 == "cells:" { print $2; exit }' "$case/log.checkMesh")
    note "$name: mesh done, $cells cells"
    note "$name: solve"
    if ! ./run-simulation.sh "$case" > "$ROOT/queue.$name.solve.log" 2>&1; then
        note "$name: SOLVE FAILED (queue.$name.solve.log, $case/log.hisa)"
        STATUS[$name]="solve failed at it. $(grep -oE '^Time = [0-9]+' "$case/log.hisa" 2>/dev/null | tail -n 1 | tr -dc 0-9)"
        continue
    fi
    note "$name: solve done; $(grep -E '^(Cx|Mx)' "$ROOT/queue.$name.solve.log" | tr '\n' ' ')"
    if python3 scripts/wall_breakdown.py "$case" >> "$ROOT/queue.$name.solve.log" 2>&1; then
        # te_map.py maps the case and, when its wall sample exists, the reference case beside it
        ref=()
        [ -d "$REF/postProcessing/wallSample" ] && ref=("$REF")
        python3 scripts/te_map.py "$case" "${ref[@]}" > "results/$name/te_map.txt" 2>> "$ROOT/queue.$name.solve.log" \
            || note "$name: te_map.py failed (non-fatal)"
    else
        note "$name: wall_breakdown.py failed (non-fatal)"
    fi
    STATUS[$name]="done, $cells cells"
done

note "queue finished"
for case in "$@"; do
    name=$(basename "$case")
    note "  $name: ${STATUS[$name]:-?}"
done
