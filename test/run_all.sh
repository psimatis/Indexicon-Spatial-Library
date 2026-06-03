#!/bin/bash

set -e
cd "$(dirname "$0")"

TESTS=(
    rtree/rtree_point_2d.cpp
    rtree/rtree_point_3d.cpp
    rtree/rtree_mbr_2d.cpp
    rtree/rtree_mbr_3d.cpp
    quadtree/quadtree_point_2d.cpp
    quadtree/mxcif_quadtree_mbr_2d.cpp
    quadtree/octtree_point_3d.cpp
    kdtree/kdtree_point_2d.cpp
    kdtree/kdtree_point_3d.cpp
)

mkdir -p build

for t in "${TESTS[@]}"; do
    echo "=== Compiling $t ==="
    exe="build/${t%.cpp}.exe"
    mkdir -p "$(dirname "$exe")"
    g++ -std=c++17 -O3 -o "$exe" "$t"
    echo "=== Running $t ==="
    ./"$exe"
    echo ""
done
