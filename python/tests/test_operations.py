import numpy as np

import indexicon as idx


def test_rtree_box_insert_remove_and_knn():
    tree = idx.RTree(data="box", dims=2)
    tree.bulk_load(
        np.array([[0.0, 0.0], [5.0, 5.0]]),
        np.array([[1.0, 1.0], [6.0, 6.0]]),
        np.array([1, 2], dtype=np.int64),
    )
    tree.insert(np.array([10.0, 10.0]), np.array([11.0, 11.0]), id=3)

    assert tree.knn(np.array([10.5, 10.5]), k=1).tolist() == [3]
    assert tree.remove(np.array([0.0, 0.0]), np.array([1.0, 1.0]), 1)
    assert tree.statistics()["count"] == 2


def test_quadtree_insert_remove_and_knn():
    tree = idx.QuadTree()
    tree.bulk_load(
        np.array([[0.0, 0.0], [5.0, 5.0]]),
        np.array([1, 2], dtype=np.int64),
    )
    tree.insert(10.0, 10.0, id=3)

    assert tree.knn(np.array([9.9, 9.9]), k=1).tolist() == [3]
    assert tree.remove(0.0, 0.0, 1)
    assert tree.statistics()["count"] == 2


def test_mxcif_quadtree_insert_remove_and_knn():
    tree = idx.MXCIFQuadTree()
    tree.bulk_load(
        np.array([[0.0, 0.0], [5.0, 5.0]]),
        np.array([[1.0, 1.0], [6.0, 6.0]]),
        np.array([1, 2], dtype=np.int64),
    )
    tree.insert(10.0, 10.0, 11.0, 11.0, id=3)

    assert tree.knn(np.array([10.5, 10.5]), k=1).tolist() == [3]
    assert tree.remove(0.0, 0.0, 1.0, 1.0, 1)
    assert tree.statistics()["count"] == 2


def test_octtree_insert_remove_and_knn():
    tree = idx.OctTree()
    tree.bulk_load(
        np.array([[0.0, 0.0, 0.0], [5.0, 5.0, 5.0]]),
        np.array([1, 2], dtype=np.int64),
    )
    tree.insert(10.0, 10.0, 10.0, id=3)

    assert tree.knn(np.array([9.9, 9.9, 9.9]), k=1).tolist() == [3]
    assert tree.remove(0.0, 0.0, 0.0, 1)
    assert tree.statistics()["count"] == 2


def test_kdtree_insert_remove_and_knn():
    tree = idx.KDTree(dims=2)
    tree.bulk_load(
        np.array([[0.0, 0.0], [5.0, 5.0]]),
        np.array([1, 2], dtype=np.int64),
    )
    tree.insert(np.array([10.0, 10.0]), id=3)

    assert tree.knn(np.array([9.9, 9.9]), k=1).tolist() == [3]
    assert tree.remove(np.array([0.0, 0.0]))
    assert tree.statistics()["count"] == 2


def test_three_dimensional_factories_query_data():
    coords = np.array([[0.0, 0.0, 0.0], [2.0, 2.0, 2.0]])
    ids = np.array([10, 20], dtype=np.int64)
    minimum = np.array([-1.0, -1.0, -1.0])
    maximum = np.array([3.0, 3.0, 3.0])

    for tree in (idx.RTree(data="point", dims=3), idx.KDTree(dims=3)):
        tree.bulk_load(coords, ids)
        assert set(tree.query(minimum, maximum).tolist()) == {10, 20}
