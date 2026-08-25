#include "common.hpp"
#include "octtree.hpp"
#include <memory>
#include <utility>

namespace {

struct OctTreeIndex {
    using Node = OctreeNode<Coord, Value, 128>;
    using Point = OPoint<Coord, Value>;
    using Box = OBox<Coord, Value>;

    std::unique_ptr<Node> tree_;
    NextId next_id_ = 0;

    void bulk_load(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& coords,
        const nb::object& ids = nb::none()
    ) {
        check_coords_shape(coords, 3, "OctTree");
        const std::size_t n = static_cast<std::size_t>(coords.shape(0));
        const auto id_values = resolve_ids(ids, n, "OctTree");

        std::vector<Point> points;
        points.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            points.emplace_back(id_values[i], coords(i, 0), coords(i, 1), coords(i, 2));
        }
        tree_ = std::make_unique<Node>(points.begin(), points.end());
        update_next_id_after_bulk(ids, id_values, next_id_);
    }

    void insert(Coord x, Coord y, Coord z, const nb::object& id = nb::none()) {
        if (!tree_) {
            tree_ = std::make_unique<Node>();
        }
        tree_->insert(Point(resolve_insert_id(id, next_id_), x, y, z));
    }

    bool remove(Coord x, Coord y, Coord z, Value id) {
        if (!tree_) {
            return false;
        }
        return tree_->remove(Point(id, x, y, z));
    }

    IdArray query(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& min_corner,
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& max_corner
    ) {
        check_corner_shape(min_corner, 3, "OctTree.query min_corner");
        check_corner_shape(max_corner, 3, "OctTree.query max_corner");
        if (!tree_) {
            return empty_array();
        }
        Box box(min_corner(0), min_corner(1), min_corner(2), max_corner(0), max_corner(1), max_corner(2));
        std::vector<Value> out;
        tree_->query(box, std::back_inserter(out));
        return to_array(std::move(out));
    }

    IdArray knn(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point, std::size_t k) {
        check_corner_shape(point, 3, "OctTree.knn point");
        if (!tree_) {
            return empty_array();
        }
        std::array<Coord, 3> q{point(0), point(1), point(2)};
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

void bind_octtree(nb::module_& m) {
    nb::class_<OctTreeIndex>(m, "OctTree")
        .def(nb::init<>())
        .def("bulk_load", &OctTreeIndex::bulk_load, nb::arg("coords"), nb::arg("ids") = nb::none())
        .def("insert", &OctTreeIndex::insert, nb::arg("x"), nb::arg("y"), nb::arg("z"), nb::arg("id") = nb::none())
        .def("remove", &OctTreeIndex::remove)
        .def("query", &OctTreeIndex::query, nb::arg("min_corner"), nb::arg("max_corner"))
        .def("knn", &OctTreeIndex::knn, nb::arg("point"), nb::arg("k"))
        .def("statistics", &OctTreeIndex::statistics);
}
