cmake_minimum_required(VERSION 3.15...3.30)
project(indexicon VERSION "${SKBUILD_PROJECT_VERSION}" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_core
    bindings/module.cpp
    bindings/rtree_point.cpp
    bindings/rtree_mbr.cpp
    bindings/quadtree.cpp
    bindings/mxcif_quadtree.cpp
    bindings/octtree.cpp
    bindings/kdtree.cpp
)

target_include_directories(_core PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../indexes")
target_compile_definitions(_core PRIVATE INDEXICON_CPP_VERSION="${SKBUILD_PROJECT_VERSION_FULL}")
install(TARGETS _core LIBRARY DESTINATION indexicon)
