import numpy as np
import indexicon as idx


def test_rtree_point_2d_basic():
    tree = idx.RTree(data="point", dims=2)
    coords = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]], dtype=np.float64)
    ids = np.array([10, 20, 30], dtype=np.int64)

    tree.bulk_load(coords, ids)
    assert tree.statistics()["count"] == 3

    tree.insert(np.array([3.0, 3.0]), 40)
    assert tree.statistics()["count"] == 4

    hits = tree.query(np.array([0.5, 0.5]), np.array([2.5, 2.5]))
    assert hits.dtype == np.int64
    assert set(hits) == {20, 30}

    knn = tree.knn(np.array([0.1, 0.1]), 2)
    assert knn.dtype == np.int64
    assert set(knn) == {10, 20}

    assert tree.remove(np.array([1.0, 1.0]), 20)
    assert tree.statistics()["count"] == 3

    stats = tree.statistics()
    assert stats["count"] == 3
    assert stats["height"] >= 1


def test_rtree_box_2d_basic():
    tree = idx.RTree(data="box", dims=2)
    mins = np.array([[0.0, 0.0], [2.0, 2.0]], dtype=np.float64)
    maxs = np.array([[1.0, 1.0], [3.0, 3.0]], dtype=np.float64)
    ids = np.array([1, 2], dtype=np.int64)

    tree.bulk_load(mins, maxs, ids)
    hits = tree.query(np.array([0.5, 0.5]), np.array([2.5, 2.5]))
    assert hits.dtype == np.int64
    assert set(hits) == {1, 2}


def test_quadtree_basic():
    tree = idx.QuadTree(split="pr")
    coords = np.array([[0.0, 0.0], [1.0, 1.0], [5.0, 5.0]], dtype=np.float64)
    ids = np.array([1, 2, 3], dtype=np.int64)
    tree.bulk_load(coords, ids)

    hits = tree.query(np.array([0.0, 0.0]), np.array([2.0, 2.0]))
    assert hits.dtype == np.int64
    assert set(hits) == {1, 2}
    assert tree.split == "pr"


def test_mxcif_quadtree_basic():
    tree = idx.MXCIFQuadTree()
    mins = np.array([[0.0, 0.0], [2.0, 2.0]], dtype=np.float64)
    maxs = np.array([[1.0, 1.0], [3.0, 3.0]], dtype=np.float64)
    ids = np.array([1, 2], dtype=np.int64)
    tree.bulk_load(mins, maxs, ids)

    hits = tree.query(np.array([0.0, 0.0]), np.array([2.5, 2.5]))
    assert hits.dtype == np.int64
    assert set(hits) == {1, 2}


def test_octtree_basic():
    tree = idx.OctTree()
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=np.float64)
    ids = np.array([1, 2], dtype=np.int64)
    tree.bulk_load(coords, ids)

    hits = tree.query(np.array([0.0, 0.0, 0.0]), np.array([2.0, 2.0, 2.0]))
    assert hits.dtype == np.int64
    assert set(hits) == {1, 2}


def test_kdtree_basic():
    tree = idx.KDTree(dims=2, split="adaptive")
    coords = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]], dtype=np.float64)
    ids = np.array([1, 2, 3], dtype=np.int64)
    tree.bulk_load(coords, ids)

    hits = tree.query(np.array([0.5, 0.5]), np.array([2.5, 2.5]))
    assert hits.dtype == np.int64
    assert set(hits) == {2, 3}


def test_auto_ids_bulk_and_insert():
    tree = idx.RTree(data="point", dims=2)
    coords = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]], dtype=np.float64)
    tree.bulk_load(coords)
    tree.remove(np.array([1.0, 1.0]), 1)
    tree.insert(np.array([3.0, 3.0]))
    assert tree.statistics()["count"] == 3
    hits = tree.query(np.array([-1.0, -1.0]), np.array([10.0, 10.0]))
    assert set(hits) == {0, 2, 3}


def test_auto_ids_after_custom_bulk():
    tree = idx.RTree(data="point", dims=2)
    coords = np.array([[0.0, 0.0], [1.0, 1.0]], dtype=np.float64)
    ids = np.array([10, 20], dtype=np.int64)
    tree.bulk_load(coords, ids)
    tree.insert(np.array([2.0, 2.0]))
    hits = tree.query(np.array([-1.0, -1.0]), np.array([10.0, 10.0]))
    assert set(hits) == {10, 20, 21}


def test_auto_ids_insert_on_empty():
    tree = idx.KDTree(dims=2)
    tree.insert(np.array([0.0, 0.0]))
    tree.insert(np.array([1.0, 1.0]))
    hits = tree.query(np.array([-1.0, -1.0]), np.array([10.0, 10.0]))
    assert set(hits) == {0, 1}

