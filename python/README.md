# Indexicon for Python

Indexicon provides NumPy-native Python bindings for high-performance R-trees, Quad-trees, Oct-trees, and KD-trees backed by C++.

## Installation

Indexicon requires Python 3.10 or newer.

```bash
pip install indexicon
```

## Example

Every index supports `bulk_load`, `insert`, `remove`, `query`, `knn`, and `statistics`.

```python
import numpy as np
import indexicon as idx

tree = idx.RTree(data="point", dims=2)
coords = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]])
tree.bulk_load(coords)

tree.insert(np.array([3.0, 3.0]))
hits = tree.query(np.array([0.5, 0.5]), np.array([2.5, 2.5]))
near = tree.knn(np.array([0.1, 0.1]), k=2)

print(hits)
print(near)
print(tree.statistics())
```

More examples are available in the [`python/examples`](https://github.com/psimatis/Indexicon-Spatial-Library/tree/main/python/examples) directory.

## Supported Indexes

| Index | Data | Dimensions | Constructor |
| :--- | :--- | :---: | :--- |
| R-tree | Points | 2D / 3D | `RTree(data="point", dims=2)` |
| R-tree | Boxes | 2D / 3D | `RTree(data="box", dims=2)` |
| Quad-tree | Points | 2D | `QuadTree(split="pr")` |
| MX-CIF Quad-tree | Boxes | 2D | `MXCIFQuadTree()` |
| Oct-tree | Points | 3D | `OctTree()` |
| KD-tree | Points | 2D / 3D | `KDTree(dims=2)` |

Quad-tree splitting strategies are `"pr"`, `"pseudo_median"`, and `"longest_axis"`. KD-tree splitting strategies are `"longest_axis"`, `"adaptive"`, and `"round_robin"`.

Python bindings use `float64` coordinates, `int64` IDs, and a fixed node or bucket capacity of 128. Capacity and other template parameters can be customized through the C++ API.

## Object IDs

Use custom IDs by passing an `(N,)` NumPy array of `int64` values to `bulk_load(..., ids=ids)` or an integer to `insert(..., id=100)`. When IDs are omitted, Indexicon assigns them automatically. Automatic IDs advance past nonnegative custom IDs and are not recycled after removal.

## Citation

If you use Indexicon in a project, paper, benchmark, or product, please cite the [Indexicon paper](https://arxiv.org/abs/2606.04676).

## License

Indexicon is distributed under the MIT License.
