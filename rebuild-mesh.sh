#!/bin/bash

# Prepare the template
rm -rf openfoam/test
cp -r openfoam/template openfoam/test
mkdir -p openfoam/test/constant/triSurface

# Prepare the model
openscad -o openfoam/test/constant/triSurface/body.stl geometry/model.scad
cd openfoam/test
surfaceClean constant/triSurface/body.stl 1e-4 1e-4 constant/triSurface/body.stl

# Create the mesh
surfaceFeatureExtract
blockMesh
decomposePar
mpirun -np 8 snappyHexMesh -parallel
reconstructParMesh -constant

# Output the cell count
checkMesh 2>&1 | grep "cells:"
