#!/usr/bin/env bash
set -euo pipefail

# Rebuild mesh artifacts for an existing case directory.
#
# Usage: ./rebuild-mesh.sh [--geometry path/to/model.scad] [case-dir]
#
# If the case directory does not exist, it is initialized from openfoam/template.
# Geometry parameters are read from constant/caseProperties.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# ── configuration (env overridable) ──────────────────────────────────────────
CASE_ARG=openfoam/cases/test
NP=${NP:-12}
MAX_CELLS=${MAX_CELLS:-11000000}   # 11M cell budget; MAX_CELLS=0 disables the guard
TEMPLATE="$ROOT/openfoam/template"
GEOMETRY=${GEOMETRY:-"$ROOT/geometry/model.scad"}

# Tools resolved during validation.
OPENSCAD=
DEPENDENCIES="surfaceCheck surfaceClean surfaceFeatureExtract blockMesh \
decomposePar foamDictionary mpirun snappyHexMesh reconstructParMesh checkMesh"

# ── helpers ───────────────────────────────────────────────────────────────────
usage_error() {
    echo "error: $*" >&2
    exit 2
}

need_command() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "error: '$1' is not on PATH; source OpenFOAM v2512 and install dependencies" >&2
        exit 127
    }
}

# Resolve the OpenSCAD binary, preferring openscad-nightly and falling back to
# plain openscad. Both accept identical -o/-D flags, so callers use "$OPENSCAD".
detect_openscad() {
    local exe
    for exe in openscad-nightly openscad; do
        if command -v "$exe" >/dev/null 2>&1; then
            printf '%s\n' "$exe"
            return 0
        fi
    done
    echo "error: neither 'openscad-nightly' nor 'openscad' is on PATH" >&2
    exit 127
}

foam_scalar() {
    local file=$1
    local key=$2
    local default=$3
    local value
    value=$(awk -v key="$key" '$1 == key { v=$2; gsub(/;/, "", v); print v; exit }' "$file" 2>/dev/null || true)
    if [ -n "$value" ]; then
        printf '%s\n' "$value"
    else
        printf '%s\n' "$default"
    fi
}

strip_frozen_points_zone() {
    local mesh_dir=$1
    local zone_file="$mesh_dir/pointZones"

    [ -f "$zone_file" ] || return 0
    grep -q "frozenPoints" "$zone_file" || return 0
    rm -f "$zone_file"
}

strip_frozen_points_zones() {
    local mesh_dir

    strip_frozen_points_zone constant/polyMesh
    for mesh_dir in processor*/constant/polyMesh; do
        [ -d "$mesh_dir" ] || continue
        strip_frozen_points_zone "$mesh_dir"
    done
}

# ── pipeline stages ───────────────────────────────────────────────────────────
parse_args() {
    while [ "$#" -gt 0 ]; do
        case "$1" in
            --geometry)
                [ "$#" -ge 2 ] || usage_error "--geometry requires a path"
                GEOMETRY=$2
                shift
                ;;
            -h|--help)
                sed -n '1,8p' "$0"
                exit 0
                ;;
            -*)
                usage_error "unknown option: $1"
                ;;
            *)
                CASE_ARG=$1
                ;;
        esac
        shift
    done
}

resolve_paths() {
    case "$GEOMETRY" in
        /*) ;;
        *)  GEOMETRY="$ROOT/$GEOMETRY" ;;
    esac

    case "$CASE_ARG" in
        /*) CASE=$CASE_ARG ;;
        *)  CASE="$ROOT/$CASE_ARG" ;;
    esac
}

validate_config() {
    local exe

    [ -d "$TEMPLATE" ] || usage_error "missing template directory: $TEMPLATE"
    [ -f "$GEOMETRY" ] || usage_error "missing geometry file: $GEOMETRY"

    case "$NP" in
        ''|*[!0-9]*) usage_error "NP must be a positive integer" ;;
    esac
    [ "$NP" -gt 0 ] || usage_error "NP must be a positive integer"

    case "$MAX_CELLS" in
        ''|*[!0-9]*) usage_error "MAX_CELLS must be a non-negative integer" ;;
    esac

    for exe in $DEPENDENCIES; do
        need_command "$exe"
    done
    OPENSCAD=$(detect_openscad)
}

init_case() {
    if [ ! -d "$CASE" ]; then
        mkdir -p "$(dirname "$CASE")"
        cp -R "$TEMPLATE" "$CASE"
    fi

    local params="$CASE/constant/caseProperties"
    [ -f "$params" ] || usage_error "missing $params; create the case with scripts/create_case.py"

    D=$(foam_scalar "$params" D 80.0)
    N=$(foam_scalar "$params" N 2)
    XI=$(foam_scalar "$params" xi 45)
    LD=$(foam_scalar "$params" LD 1.0)
    TD=$(foam_scalar "$params" TD 0.02)

    echo "case      : $CASE"
    echo "geometry  : $GEOMETRY"
    echo "openscad  : $OPENSCAD"
    echo "params    : D=${D}mm N=$N xi=$XI LD=$LD TD=$TD"
    echo "parallel  : $NP ranks"
}

clean_artifacts() {
    local path time_dir

    for path in \
        "$CASE"/processor* \
        "$CASE"/postProcessing \
        "$CASE"/dynamicCode \
        "$CASE"/log.* \
        "$CASE"/constant/polyMesh \
        "$CASE"/constant/triSurface \
        "$CASE"/constant/extendedFeatureEdgeMesh
    do
        [ -e "$path" ] || continue
        rm -rf "$path"
    done

    for time_dir in "$CASE"/[1-9]* "$CASE"/0.*; do
        [ -d "$time_dir" ] || continue
        rm -rf "$time_dir"
    done
}

generate_surface() {
    mkdir -p constant/triSurface
    "$OPENSCAD" \
        -o constant/triSurface/body.stl \
        -D "D=$D; N=$N; xi=$XI; LD=$LD; TD=$TD;" \
        "$GEOMETRY"

    # surfaceClean repairs non-watertight STLs: model.scad emits a non-closed
    # surface with duplicate-vertex "illegal" triangles that it strips. But on an
    # already-clean closed surface its collapseBase pass mangles benign sub-micron
    # CGAL union slivers (and aborts on the long cylinder slivers
    # arc_stabilizers.scad produces). So only clean when surfaceCheck reports the
    # STL is not already watertight.
    surfaceCheck constant/triSurface/body.stl 2>&1 | tee log.surfaceCheck >/dev/null || true
    if grep -q "Surface has no illegal triangles" log.surfaceCheck \
       && grep -q "Surface is closed" log.surfaceCheck; then
        echo "surfaceClean: skipped (body.stl already closed with no illegal triangles)"
    else
        echo "surfaceClean: repairing body.stl"
        surfaceClean constant/triSurface/body.stl 5e-05 1e-4 constant/triSurface/body.stl
    fi

    surfaceFeatureExtract
}

build_mesh() {
    blockMesh
    decomposePar -force
    mpirun -np "$NP" snappyHexMesh -parallel 2>&1 | tee log.snappyHexMesh
    reconstructParMesh -constant 2>&1 | tee log.reconstructParMesh
    strip_frozen_points_zones

    grep -q "body" constant/polyMesh/boundary || {
        echo "error: reconstructed constant/polyMesh is missing the body patch" >&2
        exit 1
    }

    rm -rf processor*
    decomposePar -force
    strip_frozen_points_zones

    grep -q "body" processor0/constant/polyMesh/boundary || {
        echo "error: decomposed mesh is missing the body patch" >&2
        exit 1
    }
}

verify_mesh() {
    local cell_count

    checkMesh -constant -noZero 2>&1 | tee log.checkMesh
    grep -q "Mesh OK" log.checkMesh || {
        echo "error: checkMesh did not report Mesh OK" >&2
        exit 1
    }

    cell_count=$(awk '$1 == "cells:" { print $2; exit }' log.checkMesh)
    [ -n "$cell_count" ] || {
        echo "error: could not read final cell count from log.checkMesh" >&2
        exit 1
    }
    echo "mesh cells : $cell_count"
    if [ "$MAX_CELLS" -gt 0 ] && [ "$cell_count" -gt "$MAX_CELLS" ]; then
        echo "error: final mesh has $cell_count cells, exceeding MAX_CELLS=$MAX_CELLS" >&2
        exit 1
    fi

    : > case.foam
}

# ── orchestration ─────────────────────────────────────────────────────────────
main() {
    parse_args "$@"
    resolve_paths
    validate_config
    init_case
    clean_artifacts

    pushd "$CASE" >/dev/null
    generate_surface
    build_mesh
    verify_mesh
    popd >/dev/null

    echo "mesh ready : $CASE"
}

main "$@"
