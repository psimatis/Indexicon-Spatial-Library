#!/bin/bash

set -e
cd "$(dirname "$0")"

TESTS=(
    test_rtree_point
    test_rtree_mbr
    test_quadtree
    test_mxcif_quadtree
    test_octtree
    test_kdtree
)

for t in "${TESTS[@]}"; do
    echo "=== Compiling $t ==="
    g++ -std=c++17 -O3 -o "$t.exe" "$t.cpp"
    echo "=== Running $t ==="
    ./"$t.exe"
    echo ""
done
