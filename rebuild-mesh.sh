#!/usr/bin/env bash
set -euo pipefail

# Rebuild mesh artifacts for an existing case directory.
#
# Usage: ./rebuild-mesh.sh [--geometry path/to/geometry.scad] <case-dir>
#
# Geometry parameters are read from constant/caseProperties. The case must
# already exist; create one with scripts/create_case.py. A case with
# system/gmshParams is meshed by scripts/gen_gmsh.py (zoned structured Gmsh
# mesh, stitched in OpenFOAM; no OpenSCAD or snappyHexMesh). Otherwise the
# snappy pipeline runs.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# ── configuration (env overridable) ──────────────────────────────────────────
CASE_ARG=openfoam/arc
NP=${NP:-12}
MAX_CELLS=${MAX_CELLS:-11000000}   # 11M cell budget; MAX_CELLS=0 disables the guard
GEOMETRY=${GEOMETRY:-"$ROOT/geometry/arc_stabilizers.scad"}

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
                sed -n '1,9p' "$0"
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

    [ -f "$GEOMETRY" ] || usage_error "missing geometry file: $GEOMETRY"

    case "$NP" in
        ''|*[!0-9]*) usage_error "NP must be a positive integer" ;;
    esac
    [ "$NP" -gt 0 ] || usage_error "NP must be a positive integer"

    case "$MAX_CELLS" in
        ''|*[!0-9]*) usage_error "MAX_CELLS must be a non-negative integer" ;;
    esac

    MESHER=snappy
    if [ -f "$CASE/system/gmshParams" ]; then
        MESHER=gmsh
        DEPENDENCIES="${GMSH_PYTHON:-uv} gmshToFoam transformPoints foamDictionary stitchMesh createPatch decomposePar checkMesh"
    fi

    for exe in $DEPENDENCIES; do
        need_command "$exe"
    done
    [ "$MESHER" != snappy ] || OPENSCAD=$(detect_openscad)
}

init_case() {
    [ -d "$CASE" ] || usage_error "missing case directory: $CASE; create it with scripts/create_case.py"

    local params="$CASE/constant/caseProperties"
    [ -f "$params" ] || usage_error "missing $params; create the case with scripts/create_case.py"

    D=$(foam_scalar "$params" D 80.0)
    N=$(foam_scalar "$params" N 4)
    XI=$(foam_scalar "$params" xi 90)
    L=$(foam_scalar "$params" L 140.0)

    echo "case      : $CASE"
    if [ "$MESHER" != snappy ]; then
        echo "mesher    : $MESHER"
    else
        echo "geometry  : $GEOMETRY"
        echo "openscad  : $OPENSCAD"
    fi
    echo "params    : D=${D}mm N=$N xi=$XI L=${L}mm"
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

# Export one EXPORT part of the geometry to constant/triSurface/<part>.stl and
# repair it if it is not already watertight.
#
# surfaceClean strips duplicate-vertex "illegal" triangles from a non-closed
# surface. But on an already-clean closed surface its collapseBase pass mangles
# benign sub-micron CGAL union slivers (and aborts on the long cylinder slivers
# arc_stabilizers.scad produces). So only clean when surfaceCheck reports the
# STL is not already watertight.
export_part() {
    local part=$1
    local stl="constant/triSurface/$part.stl"
    local log="log.surfaceCheck.$part"

    "$OPENSCAD" \
        -o "$stl" \
        -D "D=$D; N=$N; xi=$XI; L=$L; EXPORT=\"$part\";" \
        "$GEOMETRY"

    surfaceCheck "$stl" 2>&1 | tee "$log" >/dev/null || true
    if grep -q "Surface has no illegal triangles" "$log" \
       && grep -q "Surface is closed" "$log"; then
        echo "surfaceClean: skipped ($part.stl already closed with no illegal triangles)"
    else
        echo "surfaceClean: repairing $part.stl"
        surfaceClean "$stl" 5e-05 1e-4 "$stl"
    fi
}

generate_surface() {
    mkdir -p constant/triSurface
    export_part fuselage
    if [ "${N%%.*}" -gt 0 ]; then
        export_part stabilizers
    else
        echo "stabilizers: skipped (N=$N, clean-body case)"
    fi

    surfaceFeatureExtract
}

build_mesh() {
    blockMesh
    decomposePar -force
    mpirun -np "$NP" snappyHexMesh -parallel -overwrite 2>&1 | tee log.snappyHexMesh
    reconstructParMesh -constant 2>&1 | tee log.reconstructParMesh
    strip_frozen_points_zones

    grep -q "fuselage" constant/polyMesh/boundary || {
        echo "error: reconstructed constant/polyMesh is missing the fuselage patch" >&2
        exit 1
    }

    rm -rf processor*
    decomposePar -force
    strip_frozen_points_zones

    grep -q "fuselage" processor0/constant/polyMesh/boundary || {
        echo "error: decomposed mesh is missing the fuselage patch" >&2
        exit 1
    }
}

build_gmsh_mesh() {
    # GMSH_PYTHON: a Python with the gmsh package (e.g. a venv on the CFD
    # server, which has no uv); otherwise uv provides it.
    if [ -n "${GMSH_PYTHON:-}" ]; then
        "$GMSH_PYTHON" "$ROOT/scripts/gen_mesh.py" . --volume 2>&1 | tee log.gmsh
    else
        uv run "$ROOT/scripts/gen_mesh.py" . --volume 2>&1 | tee log.gmsh
    fi
    # -keepOrientation: gmshToFoam's own orientation fix flips sound 3 um
    # wall cells on curved walls into negative volumes
    gmshToFoam -keepOrientation mesh.msh 2>&1 | tee log.gmshToFoam
    transformPoints -scale 0.001 2>&1 | tee log.transformPoints   # Gmsh works in mm
    foamDictionary constant/polyMesh/boundary -entry entry0/fuselage/type -set wall >/dev/null
    if grep -q "^ *stabilizers$" constant/polyMesh/boundary; then
        foamDictionary constant/polyMesh/boundary -entry entry0/stabilizers/type -set wall >/dev/null
    fi
    rm -f mesh.msh
    stitch_gmsh_seams
    decomposePar -force
}

# Join the zone groups (non-conformal, ESI stitchMesh integral mode, one pair
# per call). gen_mesh.py names one pair per master group and plane,
# seam_<master>_<i>_{master,slave}; every seam is flat with straight edges and
# off the walls (see scripts/gmesh/layout.py). The inner groups' seams (outer
# and wake masters) go first, so the slaves of the far zone's seams, which
# span several groups, are one patch when their turn comes. stitchMesh and
# createPatch read the fields, which have no seam entries, so 0/ is moved
# aside meanwhile.
stitch_gmsh_seams() {
    local master pair n
    mv 0 0.fields
    for master in outer wake far inner; do
        for pair in $(grep -oE "seam_${master}_[0-9]+_master" constant/polyMesh/boundary \
                | sed 's/_master$//' | sort -u -t_ -k3 -n); do
            if stitchMesh -overwrite "${pair}_master" "${pair}_slave" > "log.stitchMesh.$pair" 2>&1; then
                echo "stitched $pair"
            else
                rm -rf 0
                mv 0.fields 0
                echo "error: stitchMesh failed for $pair; see log.stitchMesh.$pair" >&2
                exit 1
            fi
            rm -rf 0   # stitchMesh writes 0/meshPhi, which the next stitch would read at the old size
        done
    done
    printf 'FoamFile { version 2.0; format ascii; class dictionary; object createPatchDict; }\npointSync false;\npatches ();\n' \
        > createPatchDict.removeEmpty
    createPatch -overwrite -dict createPatchDict.removeEmpty > log.createPatch 2>&1   # drops the emptied seams
    rm -f constant/polyMesh/meshModifiers   # stitchMesh's sliding-interface definition
    # createPatch writes meshPhi into a fresh 0/; replace it with the fields.
    rm -rf 0
    mv 0.fields 0
    n=$(grep -c "seam_" constant/polyMesh/boundary || true)
    [ "$n" -eq 0 ] || {
        echo "error: $n seam patches still have faces after stitching" >&2
        exit 1
    }
}

verify_mesh() {
    local cell_count

    checkMesh -constant -noZero 2>&1 | tee log.checkMesh
    # Hard errors stop the run; quality flags (aspect ratio, skewness,
    # non-orthogonality) are reported only: the wall-resolved fin meshes
    # carry them at the fin-root junctions and the LE apex.
    if ! grep -q "Mesh OK" log.checkMesh; then
        if grep -qE '[*]{3}(Zero or negative cell volume|Open cells|Error in face pyramids|Boundary openness)|incorrectly oriented|Number of regions: [2-9]' log.checkMesh; then
            echo "error: checkMesh found broken cells or faces; see log.checkMesh" >&2
            exit 1
        fi
        echo "warning: checkMesh quality flags (no broken cells):" >&2
        grep -E "^ \*\*\*" log.checkMesh >&2 || true
    fi

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
    if [ "$MESHER" = gmsh ]; then
        build_gmsh_mesh
    else
        generate_surface
        build_mesh
    fi
    verify_mesh
    popd >/dev/null

    echo "mesh ready : $CASE"
}

main "$@"
