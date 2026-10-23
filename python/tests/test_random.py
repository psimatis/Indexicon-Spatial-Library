import numpy as np
import pytest

import indexicon as idx

N = 5000
rng = np.random.default_rng(7)


def point_cases():
    arr_insert = lambda t, p, i: t.insert(p, id=i)
    cases = []
    for d in (2, 3):
        cases.append((f"rtree_{d}d", lambda d=d: idx.RTree("point", d), d,
                      arr_insert, lambda t, p, i: t.remove(p, i)))
        for split in ("longest_axis", "adaptive", "round_robin"):
            cases.append((f"kdtree_{d}d_{split}", lambda d=d, s=split: idx.KDTree(d, s), d,
                          arr_insert, lambda t, p, i: t.remove(p)))
    for split in ("pr", "pseudo_median", "longest_axis"):
        cases.append((f"quadtree_{split}", lambda s=split: idx.QuadTree(split=s), 2,
                      lambda t, p, i: t.insert(p[0], p[1], id=i),
                      lambda t, p, i: t.remove(p[0], p[1], i)))
    cases.append(("octtree", idx.OctTree, 3,
                  lambda t, p, i: t.insert(p[0], p[1], p[2], id=i),
                  lambda t, p, i: t.remove(p[0], p[1], p[2], i)))
    return cases


def box_cases():
    cases = []
    for d in (2, 3):
        cases.append((f"rtree_box_{d}d", lambda d=d: idx.RTree("box", d), d,
                      lambda t, lo, hi, i: t.insert(lo, hi, id=i),
                      lambda t, lo, hi, i: t.remove(lo, hi, i)))
    cases.append(("mxcif_quadtree", idx.MXCIFQuadTree, 2,
                  lambda t, lo, hi, i: t.insert(lo[0], lo[1], hi[0], hi[1], id=i),
                  lambda t, lo, hi, i: t.remove(lo[0], lo[1], hi[0], hi[1], i)))
    return cases


def check_queries(tree, alive, lo, hi, dims):
    a = np.array(sorted(alive))
    assert tree.statistics()["count"] == len(a)
    for _ in range(100):
        c = rng.uniform(0, 1000, dims)
        h = rng.uniform(1, 80, dims)
        got = tree.query(c - h, c + h)
        hit = np.all((lo[a] <= c + h) & (hi[a] >= c - h), axis=1)
        assert sorted(got.tolist()) == a[hit].tolist()
    for k in (1, 10):
        for _ in range(30):
            q = rng.uniform(0, 1000, dims)
            dist = lambda ids: np.sqrt((np.maximum(0, np.maximum(lo[ids] - q, q - hi[ids])) ** 2).sum(1))
            got = tree.knn(q, k)
            assert len(got) == k
            assert np.allclose(np.sort(dist(got)), np.sort(dist(a))[:k])


def run_workload(make, dims, insert, remove, lo, hi, args):
    nb = int(N * 0.8)
    tree = make()
    tree.bulk_load(*[x[:nb] for x in args], np.arange(nb, dtype=np.int64))
    for i in range(nb, N):
        insert(tree, *[x[i] for x in args], i)
    alive = set(range(N))
    for i in rng.choice(N, N // 20, replace=False):
        assert remove(tree, *[x[i] for x in args], int(i))
        alive.discard(int(i))
    check_queries(tree, alive, lo, hi, dims)


def make_points(dims):
    pts = rng.uniform(0, 1000, (N, dims))
    pts[int(N * 0.8)::50, 0] += 600
    return pts


@pytest.mark.parametrize("name,make,dims,insert,remove", point_cases(), ids=lambda c: c if isinstance(c, str) else "")
def test_points_match_brute_force(name, make, dims, insert, remove):
    pts = make_points(dims)
    run_workload(make, dims, insert, remove, pts, pts, (pts,))


@pytest.mark.parametrize("name,make,dims,insert,remove", box_cases(), ids=lambda c: c if isinstance(c, str) else "")
def test_boxes_match_brute_force(name, make, dims, insert, remove):
    lo = make_points(dims)
    hi = lo + rng.uniform(0, 10, (N, dims))
    run_workload(make, dims, insert, remove, lo, hi, (lo, hi))
