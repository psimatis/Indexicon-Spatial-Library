#include "common.hpp"
#include "mxcif_quadtree.hpp"
#include <memory>
#include <utility>

namespace {

struct MXCIFQuadTreeIndex {
    using Node = MXCIFQuadTreeNode<Coord, Value, 128>;
    using Box = MXCIFBox<Coord, Value>;

    std::unique_ptr<Node> tree_;
    NextId next_id_ = 0;

    void bulk_load(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& min_corners,
        const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& max_corners,
        const nb::object& ids = nb::none()
    ) {
        if (min_corners.ndim() != 2 || max_corners.ndim() != 2) {
            throw std::invalid_argument("MXCIFQuadTree.bulk_load: min_corners and max_corners must be 2D arrays");
        }
        check_coords_shape(min_corners, 2, "MXCIFQuadTree.min_corners");
        check_coords_shape(max_corners, 2, "MXCIFQuadTree.max_corners");
        const std::size_t n = static_cast<std::size_t>(min_corners.shape(0));
        if (static_cast<std::size_t>(max_corners.shape(0)) != n) {
            throw std::invalid_argument("MXCIFQuadTree.bulk_load: min_corners and max_corners must have the same number of rows");
        }
        const auto id_values = resolve_ids(ids, n, "MXCIFQuadTree");

        std::vector<Box> boxes;
        boxes.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            boxes.emplace_back(
                id_values[i],
                min_corners(i, 0),
                min_corners(i, 1),
                max_corners(i, 0),
                max_corners(i, 1)
            );
        }
        tree_ = std::make_unique<Node>(boxes.begin(), boxes.end());
        update_next_id_after_bulk(ids, id_values, next_id_);
    }

    void insert(Coord xmin, Coord ymin, Coord xmax, Coord ymax, const nb::object& id = nb::none()) {
        if (!tree_) {
            tree_ = std::make_unique<Node>();
        }
        tree_->insert(Box(resolve_insert_id(id, next_id_), xmin, ymin, xmax, ymax));
    }

    bool remove(Coord xmin, Coord ymin, Coord xmax, Coord ymax, Value id) {
        if (!tree_) {
            return false;
        }
        return tree_->remove(Box(id, xmin, ymin, xmax, ymax));
    }

    IdArray query(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& min_corner,
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& max_corner
    ) {
        check_corner_shape(min_corner, 2, "MXCIFQuadTree.query min_corner");
        check_corner_shape(max_corner, 2, "MXCIFQuadTree.query max_corner");
        if (!tree_) {
            return empty_array();
        }
        Box box(0, min_corner(0), min_corner(1), max_corner(0), max_corner(1));
        std::vector<Value> out;
        tree_->query(box, std::back_inserter(out));
        return to_array(std::move(out));
    }

    IdArray knn(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point, std::size_t k) {
        check_corner_shape(point, 2, "MXCIFQuadTree.knn point");
        if (!tree_) {
            return empty_array();
        }
        std::array<Coord, 2> q{point(0), point(1)};
        std::vector<Value> out;
        tree_->knnQuery(q, k, std::back_inserter(out));
        return to_array(std::move(out));
    }

    nb::dict statistics() const {
        if (!tree_) {
            return stats_dict(0, 0, 0, 0, 0);
        }
        const auto s = tree_->getStatistics();
        return stats_dict(s.height, s.numPoints, s.numLeaves, s.numInternalNodes, s.sizeBytes);
    }

};

}  // namespace

void bind_mxcif_quadtree(nb::module_& m) {
    nb::class_<MXCIFQuadTreeIndex>(m, "MXCIFQuadTree")
        .def(nb::init<>())
        .def("bulk_load", &MXCIFQuadTreeIndex::bulk_load, nb::arg("min_corners"), nb::arg("max_corners"), nb::arg("ids") = nb::none())
        .def("insert", &MXCIFQuadTreeIndex::insert, nb::arg("xmin"), nb::arg("ymin"), nb::arg("xmax"), nb::arg("ymax"), nb::arg("id") = nb::none())
        .def("remove", &MXCIFQuadTreeIndex::remove)
        .def("query", &MXCIFQuadTreeIndex::query, nb::arg("min_corner"), nb::arg("max_corner"))
        .def("knn", &MXCIFQuadTreeIndex::knn, nb::arg("point"), nb::arg("k"))
        .def("statistics", &MXCIFQuadTreeIndex::statistics);
}
