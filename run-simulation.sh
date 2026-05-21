#!/usr/bin/env bash
set -euo pipefail

# Run rhoCentralFoam on a decomposed case, reconstruct the final time step, and
# emit coefficient CSV output.
#
# Usage:
#   ./run-simulation.sh [--dry-run] [case-dir]

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CASE_ARG=openfoam/test
DRY_RUN=0
NP=${NP:-6}

usage() {
    sed -n '1,10p' "$0" >&2
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --dry-run)
            DRY_RUN=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            echo "error: unknown option: $1" >&2
            usage
            exit 2
            ;;
        *)
            CASE_ARG=$1
            ;;
    esac
    shift
done

case "$CASE_ARG" in
    /*) CASE=$CASE_ARG ;;
    *)  CASE="$ROOT/$CASE_ARG" ;;
esac

need_command() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "error: '$1' is not on PATH; source OpenFOAM v2512 first" >&2
        exit 127
    }
}

require_file() {
    [ -f "$1" ] || {
        echo "error: missing $1" >&2
        exit 1
    }
}

processor_count() {
    local count=0
    local dir
    for dir in "$CASE"/processor*; do
        [ -d "$dir" ] || continue
        count=$((count + 1))
    done
    printf '%s\n' "$count"
}

case "$NP" in
    ''|*[!0-9]*)
        echo "error: NP must be a positive integer" >&2
        exit 2
        ;;
esac
[ "$NP" -gt 0 ] || {
    echo "error: NP must be a positive integer" >&2
    exit 2
}

for exe in mpirun rhoCentralFoam foamListTimes reconstructPar; do
    need_command "$exe"
done

[ -d "$CASE/processor0" ] || {
    echo "error: $CASE is not decomposed; run ./rebuild-mesh.sh $CASE_ARG first" >&2
    exit 1
}

PROC_COUNT=$(processor_count)
[ "$PROC_COUNT" -eq "$NP" ] || {
    echo "error: $CASE has $PROC_COUNT processor directories but NP=$NP; rerun ./rebuild-mesh.sh with the same NP" >&2
    exit 1
}

require_file "$CASE/constant/polyMesh/boundary"
require_file "$CASE/processor0/constant/polyMesh/boundary"

grep -q "body" "$CASE/constant/polyMesh/boundary" || {
    echo "error: $CASE/constant/polyMesh/boundary is missing the body patch" >&2
    exit 1
}
grep -q "body" "$CASE/processor0/constant/polyMesh/boundary" || {
    echo "error: $CASE/processor0/constant/polyMesh/boundary is missing the body patch" >&2
    exit 1
}

for field in U p T k omega nut alphat; do
    require_file "$CASE/processor0/0/$field"
done

pushd "$CASE" >/dev/null

# Clean prior run outputs while preserving the final mesh and initial fields.
foamListTimes -rm -processor >/dev/null 2>&1 || true
foamListTimes -rm            >/dev/null 2>&1 || true
rm -rf postProcessing log.rhoCentralFoam log.rhoCentralFoam.dryRun log.reconstructPar

if [ "$DRY_RUN" -eq 1 ]; then
    if ! mpirun -np "$NP" rhoCentralFoam -parallel -dry-run 2>&1 | tee log.rhoCentralFoam.dryRun; then
        if ! grep -q "MPI_ERR_TRUNCATE" log.rhoCentralFoam.dryRun; then
            exit 1
        fi
        {
            echo
            echo "warning: parallel rhoCentralFoam -dry-run hit MPI_ERR_TRUNCATE; retrying serial dry-run on reconstructed mesh"
        } | tee -a log.rhoCentralFoam.dryRun
        rhoCentralFoam -dry-run 2>&1 | tee -a log.rhoCentralFoam.dryRun
    fi
    popd >/dev/null
    echo "dry-run OK : $CASE"
    exit 0
fi

mpirun -np "$NP" rhoCentralFoam -parallel 2>&1 | tee log.rhoCentralFoam
reconstructPar -latestTime 2>&1 | tee log.reconstructPar
: > case.foam

popd >/dev/null

python3 "$ROOT/scripts/post_process.py" "$CASE"
