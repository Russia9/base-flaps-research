#!/usr/bin/env bash
set -euo pipefail

# Rebuild mesh artifacts for an existing case directory.
#
# Usage: ./rebuild-mesh.sh [case-dir]
#
# If the case directory does not exist, it is initialized from openfoam/template.
# Geometry parameters are read from constant/caseProperties.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CASE_ARG=${1:-openfoam/test}
NP=${NP:-6}
MAX_CELLS=${MAX_CELLS:-2000000}

case "$CASE_ARG" in
    /*) CASE=$CASE_ARG ;;
    *)  CASE="$ROOT/$CASE_ARG" ;;
esac

TEMPLATE="$ROOT/openfoam/template"
GEOMETRY="$ROOT/geometry/model.scad"

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

strip_empty_frozen_points_zone() {
    local mesh_dir=$1
    local zone_file="$mesh_dir/pointZones"

    [ -f "$zone_file" ] || return 0
    grep -q "names[[:space:]]*( frozenPoints );" "$zone_file" || return 0
    grep -q "pointLabels[[:space:]]*List<label>[[:space:]]*0;" "$zone_file" || return 0
    rm -f "$zone_file"
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
echo "geometry  : D=${D}mm N=$N xi=$XI LD=$LD TD=$TD"
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
surfaceClean constant/triSurface/body.stl 1e-4 1e-4 constant/triSurface/body.stl
surfaceFeatureExtract
blockMesh
decomposePar -force
mpirun -np "$NP" snappyHexMesh -parallel -overwrite 2>&1 | tee log.snappyHexMesh
reconstructParMesh -constant 2>&1 | tee log.reconstructParMesh
strip_empty_frozen_points_zone constant/polyMesh

if ! grep -q "body" constant/polyMesh/boundary; then
    echo "error: reconstructed constant/polyMesh is missing the body patch" >&2
    exit 1
fi

rm -rf processor*
decomposePar -force

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
