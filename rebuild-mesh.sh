#!/usr/bin/env bash
set -euo pipefail

# Rebuild mesh artifacts for an existing case directory.
#
# Usage: ./rebuild-mesh.sh [--geometry path/to/model.scad] [case-dir]
#
# If the case directory does not exist, it is initialized from openfoam/template.
# Geometry parameters are read from constant/caseProperties.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CASE_ARG=openfoam/cases/test
NP=${NP:-6}
MAX_CELLS=${MAX_CELLS:-2000000}

TEMPLATE="$ROOT/openfoam/template"
GEOMETRY=${GEOMETRY:-"$ROOT/geometry/model.scad"}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --geometry)
            [ "$#" -ge 2 ] || {
                echo "error: --geometry requires a path" >&2
                exit 2
            }
            GEOMETRY=$2
            shift
            ;;
        -h|--help)
            sed -n '1,8p' "$0"
            exit 0
            ;;
        -*)
            echo "error: unknown option: $1" >&2
            exit 2
            ;;
        *)
            CASE_ARG=$1
            ;;
    esac
    shift
done

case "$GEOMETRY" in
    /*) ;;
    *)  GEOMETRY="$ROOT/$GEOMETRY" ;;
esac

case "$CASE_ARG" in
    /*) CASE=$CASE_ARG ;;
    *)  CASE="$ROOT/$CASE_ARG" ;;
esac

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

[ -d "$TEMPLATE" ] || usage_error "missing template directory: $TEMPLATE"
[ -f "$GEOMETRY" ] || usage_error "missing geometry file: $GEOMETRY"

case "$NP" in
    ''|*[!0-9]*) usage_error "NP must be a positive integer" ;;
esac
[ "$NP" -gt 0 ] || usage_error "NP must be a positive integer"

case "$MAX_CELLS" in
    ''|*[!0-9]*) usage_error "MAX_CELLS must be a non-negative integer" ;;
esac

for exe in openscad surfaceClean surfaceFeatureExtract blockMesh decomposePar foamDictionary mpirun snappyHexMesh reconstructParMesh checkMesh; do
    need_command "$exe"
done

if [ ! -d "$CASE" ]; then
    mkdir -p "$(dirname "$CASE")"
    cp -R "$TEMPLATE" "$CASE"
fi

PARAMS="$CASE/constant/caseProperties"
[ -f "$PARAMS" ] || usage_error "missing $PARAMS; create the case with scripts/create_case.py"

D=$(foam_scalar "$PARAMS" D 80.0)
N=$(foam_scalar "$PARAMS" N 2)
XI=$(foam_scalar "$PARAMS" xi 45)
LD=$(foam_scalar "$PARAMS" LD 1.0)
TD=$(foam_scalar "$PARAMS" TD 0.02)

echo "case      : $CASE"
echo "geometry  : $GEOMETRY"
echo "params    : D=${D}mm N=$N xi=$XI LD=$LD TD=$TD"
echo "parallel  : $NP ranks"

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

mkdir -p "$CASE/constant/triSurface"
openscad \
    -o "$CASE/constant/triSurface/body.stl" \
    -D "D=$D; N=$N; xi=$XI; LD=$LD; TD=$TD;" \
    "$GEOMETRY"

pushd "$CASE" >/dev/null

foamDictionary system/decomposeParDict -entry numberOfSubdomains -set "$NP" >/dev/null

# surfaceClean repairs non-watertight STLs: model.scad emits a non-closed surface
# with duplicate-vertex "illegal" triangles that it strips. But on an already-clean
# closed surface its collapseBase pass mangles benign sub-micron CGAL union slivers
# (and aborts on the long cylinder slivers arc_stabilizers.scad produces). So only
# clean when surfaceCheck reports the STL is not already watertight.
surfaceCheck constant/triSurface/body.stl 2>&1 | tee log.surfaceCheck >/dev/null || true
if grep -q "Surface has no illegal triangles" log.surfaceCheck \
   && grep -q "Surface is closed" log.surfaceCheck; then
    echo "surfaceClean: skipped (body.stl already closed with no illegal triangles)"
else
    echo "surfaceClean: repairing body.stl"
    surfaceClean constant/triSurface/body.stl 5e-05 1e-4 constant/triSurface/body.stl
fi

surfaceFeatureExtract
blockMesh
decomposePar -force
mpirun -np "$NP" snappyHexMesh -parallel -overwrite 2>&1 | tee log.snappyHexMesh
reconstructParMesh -constant 2>&1 | tee log.reconstructParMesh
strip_frozen_points_zones

if ! grep -q "body" constant/polyMesh/boundary; then
    echo "error: reconstructed constant/polyMesh is missing the body patch" >&2
    exit 1
fi

rm -rf processor*
decomposePar -force
strip_frozen_points_zones

if ! grep -q "body" processor0/constant/polyMesh/boundary; then
    echo "error: decomposed mesh is missing the body patch" >&2
    exit 1
fi

checkMesh -constant -noZero 2>&1 | tee log.checkMesh
grep -q "Mesh OK" log.checkMesh || {
    echo "error: checkMesh did not report Mesh OK" >&2
    exit 1
}

CELL_COUNT=$(awk '$1 == "cells:" { print $2; exit }' log.checkMesh)
[ -n "$CELL_COUNT" ] || {
    echo "error: could not read final cell count from log.checkMesh" >&2
    exit 1
}
echo "mesh cells : $CELL_COUNT"
if [ "$MAX_CELLS" -gt 0 ] && [ "$CELL_COUNT" -gt "$MAX_CELLS" ]; then
    echo "error: final mesh has $CELL_COUNT cells, exceeding MAX_CELLS=$MAX_CELLS" >&2
    exit 1
fi

: > case.foam

popd >/dev/null

echo "mesh ready : $CASE"
