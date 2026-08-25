"""2D KD-tree: bulk-load 80%, insert 20%, delete 5%, range query, kNN."""

import numpy as np
import indexicon as idx

# --- configuration ---
N = 10_000
BULK_FRAC = 0.80
DELETE_FRAC = 0.05
KNN_K = 10
SEED = 42
DIMS = 2

SPLIT = "longest_axis"
# SPLIT = "adaptive"
# SPLIT = "round_robin"

rng = np.random.default_rng(SEED)
n_bulk = int(N * BULK_FRAC)
n_delete = int(N * DELETE_FRAC)

# --- random data ---
coords = rng.uniform(0.0, 100.0, (N, DIMS))
ids = np.arange(N, dtype=np.int64)

# --- build index ---
tree = idx.KDTree(dims=DIMS, split=SPLIT)

# --- bulk load (80%) ---
tree.bulk_load(coords[:n_bulk], ids[:n_bulk])

# --- insert (20%) ---
for i in range(n_bulk, N):
    tree.insert(coords[i], int(ids[i]))

# --- delete (5%) --- remove matches by point position (no id)
for i in range(n_delete):
    tree.remove(coords[i])

# --- range query ---
min_corner = np.array([20.0, 20.0])
max_corner = np.array([80.0, 80.0])
hits = tree.query(min_corner, max_corner)

# --- kNN query ---
query_point = rng.uniform(0.0, 100.0, DIMS)
near = tree.knn(query_point, k=KNN_K)

# --- results ---
stats = tree.statistics()
print(f"count={stats['count']}, range_hits={len(hits)}, knn_hits={len(near)}")
