import numpy as np
import pytest

import indexicon as idx


def test_automatic_id_advances_past_explicit_insert_id():
    tree = idx.RTree(data="point", dims=2)
    tree.insert(np.array([0.0, 0.0]), id=0)
    tree.insert(np.array([1.0, 1.0]))

    hits = tree.query(np.array([-1.0, -1.0]), np.array([2.0, 2.0]))

    assert sorted(hits.tolist()) == [0, 1]


def test_kdtree_query_rejects_corner_with_wrong_dimension():
    tree = idx.KDTree(dims=2)
    tree.bulk_load(np.array([[0.0, 0.0], [1.0, 1.0]]))

    with pytest.raises(ValueError, match="expected length 2"):
        tree.query(np.array([0.0, 0.0, 0.0]), np.array([2.0, 2.0]))


@pytest.mark.parametrize(
    ("tree", "load"),
    [
        (
            idx.RTree(data="point", dims=2),
            lambda tree, offset: tree.bulk_load(
                np.array([[offset, offset], [offset + 1, offset + 1]])
            ),
        ),
        (
            idx.RTree(data="box", dims=2),
            lambda tree, offset: tree.bulk_load(
                np.array([[offset, offset]]),
                np.array([[offset + 1, offset + 1]]),
            ),
        ),
        (
            idx.KDTree(dims=2),
            lambda tree, offset: tree.bulk_load(
                np.array([[offset, offset], [offset + 1, offset + 1]])
            ),
        ),
    ],
)
def test_bulk_load_replaces_existing_tree(tree, load):
    load(tree, 0.0)
    load(tree, 10.0)

    hits = tree.query(np.array([9.0, 9.0]), np.array([12.0, 12.0]))

    assert tree.statistics()["count"] in {1, 2}
    assert set(hits.tolist()) == set(range(tree.statistics()["count"]))


def test_explicit_id_10_then_automatic_id_11():
    tree = idx.RTree(data="point", dims=2)
    tree.insert(np.array([0.0, 0.0]), id=10)
    tree.insert(np.array([1.0, 1.0]))

    hits = tree.query(np.array([-1.0, -1.0]), np.array([10.0, 10.0]))
    assert sorted(hits.tolist()) == [10, 11]


def test_all_negative_bulk_then_automatic_id_zero():
    tree = idx.RTree(data="point", dims=2)
    coords = np.array([[0.0, 0.0], [1.0, 1.0]], dtype=np.float64)
    ids = np.array([-1, -2], dtype=np.int64)
    tree.bulk_load(coords, ids)
    tree.insert(np.array([2.0, 2.0]))

    hits = tree.query(np.array([-1.0, -1.0]), np.array([10.0, 10.0]))
    assert sorted(hits.tolist()) == [-2, -1, 0]


def test_automatic_id_exhausted_after_explicit_int64_max():
    tree = idx.RTree(data="point", dims=2)
    tree.insert(np.array([0.0, 0.0]), id=np.int64(2**63 - 1))

    with pytest.raises(OverflowError, match="automatic id exhausted"):
        tree.insert(np.array([1.0, 1.0]))
