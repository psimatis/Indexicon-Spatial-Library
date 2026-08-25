#include "common.hpp"
#include <optional>
#include <utility>

#define RTree RTreePointImpl
#include "rtree_point.hpp"
#undef RTree

namespace {

template <std::size_t Dims>
struct RTreePointIndex {
    using Tree = RTreePointImpl<Value, double, Dims>;
    using Point = typename Tree::PointType;
    using Box = typename Tree::BoxType;

    std::optional<Tree> tree_{std::in_place};
    NextId next_id_ = 0;

    void bulk_load(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& coords,
        const nb::object& ids = nb::none()
    ) {
        check_coords_shape(coords, Dims, "RTreePoint");
        const std::size_t n = static_cast<std::size_t>(coords.shape(0));
        const auto id_values = resolve_ids(ids, n, "RTreePoint");

        std::vector<std::pair<Point, Value>> entries;
        entries.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            Point p;
            for (std::size_t d = 0; d < Dims; ++d) {
                p[d] = coords(i, static_cast<int>(d));
            }
            entries.emplace_back(p, id_values[i]);
        }
        tree_.emplace(entries.begin(), entries.end());
        update_next_id_after_bulk(ids, id_values, next_id_);
    }

    void insert(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point,
        const nb::object& id = nb::none()
    ) {
        if (static_cast<std::size_t>(point.shape(0)) != Dims) {
            throw std::invalid_argument("RTreePoint.insert: point must have length equal to dims");
        }
        Point p;
        for (std::size_t d = 0; d < Dims; ++d) {
            p[d] = point(d);
        }
        tree_.value().insert(p, resolve_insert_id(id, next_id_));
    }

    bool remove(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point,
        Value id
    ) {
        if (static_cast<std::size_t>(point.shape(0)) != Dims) {
            throw std::invalid_argument("RTreePoint.remove: point must have length equal to dims");
        }
        Point p;
        for (std::size_t d = 0; d < Dims; ++d) {
            p[d] = point(d);
        }
        return tree_.value().remove(p, id);
    }

    IdArray query(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& min_corner,
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& max_corner
    ) const {
        if (static_cast<std::size_t>(min_corner.shape(0)) != Dims ||
            static_cast<std::size_t>(max_corner.shape(0)) != Dims) {
            throw std::invalid_argument("RTreePoint.query: box corners must have length equal to dims");
        }
        Box box;
        for (std::size_t d = 0; d < Dims; ++d) {
            box.minCorner[d] = min_corner(d);
            box.maxCorner[d] = max_corner(d);
        }
        std::vector<Value> out;
        tree_.value().query(box, std::back_inserter(out));
        return to_array(std::move(out));
    }

    IdArray knn(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point,
        std::size_t k
    ) const {
        if (static_cast<std::size_t>(point.shape(0)) != Dims) {
            throw std::invalid_argument("RTreePoint.knn: point must have length equal to dims");
        }
        Point p;
        for (std::size_t d = 0; d < Dims; ++d) {
            p[d] = point(d);
        }
        std::vector<Value> out;
        tree_.value().knnQuery(p, k, std::back_inserter(out));
        return to_array(std::move(out));
    }

    nb::dict statistics() const {
        const auto s = tree_.value().getStatistics();
        return stats_dict(s.height, s.numPoints, s.numLeaves, s.numInternalNodes, s.sizeBytes);
    }
};

}  // namespace

void bind_rtree_point(nb::module_& m) {
    nb::class_<RTreePointIndex<2>>(m, "RTreePoint2D")
        .def(nb::init<>())
        .def("bulk_load", &RTreePointIndex<2>::bulk_load, nb::arg("coords"), nb::arg("ids") = nb::none())
        .def("insert", &RTreePointIndex<2>::insert, nb::arg("point"), nb::arg("id") = nb::none())
        .def("remove", &RTreePointIndex<2>::remove)
        .def("query", &RTreePointIndex<2>::query, nb::arg("min_corner"), nb::arg("max_corner"))
        .def("knn", &RTreePointIndex<2>::knn, nb::arg("point"), nb::arg("k"))
        .def("statistics", &RTreePointIndex<2>::statistics);

    nb::class_<RTreePointIndex<3>>(m, "RTreePoint3D")
        .def(nb::init<>())
        .def("bulk_load", &RTreePointIndex<3>::bulk_load, nb::arg("coords"), nb::arg("ids") = nb::none())
        .def("insert", &RTreePointIndex<3>::insert, nb::arg("point"), nb::arg("id") = nb::none())
        .def("remove", &RTreePointIndex<3>::remove)
        .def("query", &RTreePointIndex<3>::query, nb::arg("min_corner"), nb::arg("max_corner"))
        .def("knn", &RTreePointIndex<3>::knn, nb::arg("point"), nb::arg("k"))
        .def("statistics", &RTreePointIndex<3>::statistics);
}
