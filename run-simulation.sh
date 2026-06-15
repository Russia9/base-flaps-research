#!/usr/bin/env bash
set -euo pipefail

# Run hisa on a decomposed case, reconstruct the final time step, and
# emit coefficient CSV output.
#
# Usage:
#   ./run-simulation.sh [--dry-run] [case-dir]
#
# --dry-run validates the decomposed mesh in parallel (checkMesh); hisa itself
# has no -dry-run option.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# ── configuration (env overridable) ──────────────────────────────────────────
CASE_ARG=openfoam/test
DRY_RUN=0
NP=${NP:-12}
DEPENDENCIES="mpirun hisa checkMesh foamListTimes reconstructPar"
FIELDS="U p T k omega nut alphat"

# ── helpers ───────────────────────────────────────────────────────────────────
usage() {
    sed -n '1,10p' "$0" >&2
}

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

require_body_patch() {
    grep -q "type.*wall" "$1" || {
        echo "error: $1 has no wall-type patch" >&2
        exit 1
    }
}

# ── pipeline stages ───────────────────────────────────────────────────────────
parse_args() {
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
}

resolve_paths() {
    case "$CASE_ARG" in
        /*) CASE=$CASE_ARG ;;
        *)  CASE="$ROOT/$CASE_ARG" ;;
    esac
}

validate_config() {
    local exe

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

    for exe in $DEPENDENCIES; do
        need_command "$exe"
    done
}

preflight() {
    local proc_count field

    [ -d "$CASE/processor0" ] || {
        echo "error: $CASE is not decomposed; run ./rebuild-mesh.sh $CASE_ARG first" >&2
        exit 1
    }

    proc_count=$(processor_count)
    [ "$proc_count" -eq "$NP" ] || {
        echo "error: $CASE has $proc_count processor directories but NP=$NP; rerun ./rebuild-mesh.sh with the same NP" >&2
        exit 1
    }

    require_file "$CASE/constant/polyMesh/boundary"
    require_file "$CASE/processor0/constant/polyMesh/boundary"
    require_body_patch "$CASE/constant/polyMesh/boundary"
    require_body_patch "$CASE/processor0/constant/polyMesh/boundary"

    for field in $FIELDS; do
        require_file "$CASE/processor0/0/$field"
    done
}

clean_outputs() {
    # Clean prior run outputs while preserving the final mesh and initial fields.
    foamListTimes -rm -processor >/dev/null 2>&1 || true
    foamListTimes -rm            >/dev/null 2>&1 || true
    rm -rf postProcessing log.hisa log.checkMesh.dryRun log.reconstructPar
}

run_dry_run() {
    # hisa has no -dry-run; validate the decomposed mesh/patches instead, which
    # is the usual startup failure mode. checkMesh exits non-zero on a failed
    # check, which (set -e) aborts the dry-run as intended.
    mpirun -np "$NP" checkMesh -parallel 2>&1 | tee log.checkMesh.dryRun
}

run_solver() {
    mpirun -np "$NP" hisa -parallel 2>&1 | tee log.hisa
    reconstructPar -latestTime 2>&1 | tee log.reconstructPar
    : > case.foam
}

# ── orchestration ─────────────────────────────────────────────────────────────
main() {
    parse_args "$@"
    resolve_paths
    validate_config
    preflight

    pushd "$CASE" >/dev/null
    clean_outputs

    if [ "$DRY_RUN" -eq 1 ]; then
        run_dry_run
        popd >/dev/null
        echo "dry-run OK : $CASE (mesh/decomposition checked)"
        exit 0
    fi

    run_solver
    popd >/dev/null

    python3 "$ROOT/scripts/post_process.py" "$CASE"
}

main "$@"
