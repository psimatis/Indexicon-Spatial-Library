#include "common.hpp"
#include <optional>
#include <utility>

#define RTree RTreeBoxImpl
#include "rtree_mbr.hpp"
#undef RTree

namespace {

template <std::size_t Dims>
struct RTreeBoxIndex {
    using Tree = RTreeBoxImpl<Value, double, Dims>;
    using Point = typename Tree::PointType;
    using Box = typename Tree::BoxType;

    std::optional<Tree> tree_{std::in_place};
    NextId next_id_ = 0;

    void bulk_load(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& min_corners,
        const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& max_corners,
        const nb::object& ids = nb::none()
    ) {
        if (min_corners.ndim() != 2 || max_corners.ndim() != 2) {
            throw std::invalid_argument("RTreeBox.bulk_load: min_corners and max_corners must be 2D arrays");
        }
        check_coords_shape(min_corners, Dims, "RTreeBox.min_corners");
        check_coords_shape(max_corners, Dims, "RTreeBox.max_corners");
        const std::size_t n = static_cast<std::size_t>(min_corners.shape(0));
        if (static_cast<std::size_t>(max_corners.shape(0)) != n) {
            throw std::invalid_argument("RTreeBox.bulk_load: min_corners and max_corners must have the same number of rows");
        }
        const auto id_values = resolve_ids(ids, n, "RTreeBox");

        std::vector<std::pair<Box, Value>> entries;
        entries.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            Box box;
            for (std::size_t d = 0; d < Dims; ++d) {
                box.minCorner[d] = min_corners(i, static_cast<int>(d));
                box.maxCorner[d] = max_corners(i, static_cast<int>(d));
            }
            entries.emplace_back(box, id_values[i]);
        }
        tree_.emplace(entries.begin(), entries.end());
        update_next_id_after_bulk(ids, id_values, next_id_);
    }

    void insert(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& min_corner,
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& max_corner,
        const nb::object& id = nb::none()
    ) {
        if (static_cast<std::size_t>(min_corner.shape(0)) != Dims ||
            static_cast<std::size_t>(max_corner.shape(0)) != Dims) {
            throw std::invalid_argument("RTreeBox.insert: box corners must have length equal to dims");
        }
        Box box;
        for (std::size_t d = 0; d < Dims; ++d) {
            box.minCorner[d] = min_corner(d);
            box.maxCorner[d] = max_corner(d);
        }
        tree_.value().insert(box, resolve_insert_id(id, next_id_));
    }

    bool remove(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& min_corner,
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& max_corner,
        Value id
    ) {
        if (static_cast<std::size_t>(min_corner.shape(0)) != Dims ||
            static_cast<std::size_t>(max_corner.shape(0)) != Dims) {
            throw std::invalid_argument("RTreeBox.remove: box corners must have length equal to dims");
        }
        Box box;
        for (std::size_t d = 0; d < Dims; ++d) {
            box.minCorner[d] = min_corner(d);
            box.maxCorner[d] = max_corner(d);
        }
        return tree_.value().remove(box, id);
    }

    IdArray query(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& min_corner,
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& max_corner
    ) const {
        if (static_cast<std::size_t>(min_corner.shape(0)) != Dims ||
            static_cast<std::size_t>(max_corner.shape(0)) != Dims) {
            throw std::invalid_argument("RTreeBox.query: box corners must have length equal to dims");
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
            throw std::invalid_argument("RTreeBox.knn: point must have length equal to dims");
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

void bind_rtree_mbr(nb::module_& m) {
    nb::class_<RTreeBoxIndex<2>>(m, "RTreeBox2D")
        .def(nb::init<>())
        .def("bulk_load", &RTreeBoxIndex<2>::bulk_load, nb::arg("min_corners"), nb::arg("max_corners"), nb::arg("ids") = nb::none())
        .def("insert", &RTreeBoxIndex<2>::insert, nb::arg("min_corner"), nb::arg("max_corner"), nb::arg("id") = nb::none())
        .def("remove", &RTreeBoxIndex<2>::remove)
        .def("query", &RTreeBoxIndex<2>::query, nb::arg("min_corner"), nb::arg("max_corner"))
        .def("knn", &RTreeBoxIndex<2>::knn, nb::arg("point"), nb::arg("k"))
        .def("statistics", &RTreeBoxIndex<2>::statistics);

    nb::class_<RTreeBoxIndex<3>>(m, "RTreeBox3D")
        .def(nb::init<>())
        .def("bulk_load", &RTreeBoxIndex<3>::bulk_load, nb::arg("min_corners"), nb::arg("max_corners"), nb::arg("ids") = nb::none())
        .def("insert", &RTreeBoxIndex<3>::insert, nb::arg("min_corner"), nb::arg("max_corner"), nb::arg("id") = nb::none())
        .def("remove", &RTreeBoxIndex<3>::remove)
        .def("query", &RTreeBoxIndex<3>::query, nb::arg("min_corner"), nb::arg("max_corner"))
        .def("knn", &RTreeBoxIndex<3>::knn, nb::arg("point"), nb::arg("k"))
        .def("statistics", &RTreeBoxIndex<3>::statistics);
}
