"""Indexicon spatial index bindings."""

from __future__ import annotations
from importlib.metadata import version
from indexicon._core import (
    KDTree2D as _KDTree2D,
    KDTree3D as _KDTree3D,
    MXCIFQuadTree,
    OctTree,
    QuadTree,
    RTreeBox2D as _RTreeBox2D,
    RTreeBox3D as _RTreeBox3D,
    RTreePoint2D as _RTreePoint2D,
    RTreePoint3D as _RTreePoint3D,
    __indexicon_cpp_version__,
)

__version__ = version("indexicon")

__all__ = [
    "KDTree",
    "MXCIFQuadTree",
    "OctTree",
    "QuadTree",
    "RTree",
    "__indexicon_cpp_version__",
]


def RTree(data: str = "point", dims: int = 2):
    """Create an R-tree index. data: point or box. dims: 2 or 3."""
    if data == "point":
        if dims == 2:
            return _RTreePoint2D()
        if dims == 3:
            return _RTreePoint3D()
    elif data == "box":
        if dims == 2:
            return _RTreeBox2D()
        if dims == 3:
            return _RTreeBox3D()
    raise ValueError('RTree: data must be "point" or "box", dims must be 2 or 3')


def KDTree(dims: int = 2, split: str = "longest_axis"):
    """Create a KD-tree. split: adaptive, round_robin, longest_axis."""
    if dims == 2:
        return _KDTree2D(split)
    if dims == 3:
        return _KDTree3D(split)
    raise ValueError("KDTree supports dims 2 or 3")
