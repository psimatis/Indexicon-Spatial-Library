#include "common.hpp"
#include "kdtree.hpp"
#include <memory>
#include <optional>
#include <string>
#include <utility>

namespace {

KDSplit parse_kd_split(const std::string& split) {
    if (split == "adaptive" || split == "ADAPTIVE")
        return KDSplit::ADAPTIVE;
    if (split == "round_robin" || split == "ROUND_ROBIN")
        return KDSplit::ROUND_ROBIN;
    if (split == "longest_axis" || split == "LONGEST_AXIS")
        return KDSplit::LONGEST_AXIS;
    throw invalid_argument("KDTree: split must be one of 'adaptive', 'round_robin', 'longest_axis'");
}

template <std::size_t Dims, typename Point>
Point read_kd_point(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point) {
    if (static_cast<std::size_t>(point.shape(0)) != Dims)
        throw invalid_argument("KDTree: point must have length equal to dims");
    Point p;
    for (std::size_t d = 0; d < Dims; ++d)
        p[d] = point(d);
    return p;
}

struct KDTreeBase {
    virtual ~KDTreeBase() = default;
    virtual void bulk_load(const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& coords, const nb::object& ids = nb::none()) = 0;
    virtual void insert(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point, const nb::object& id = nb::none()) = 0;
    virtual bool remove(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point) = 0;
    virtual IdArray query(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& min_corner, const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& max_corner) = 0;
    virtual IdArray knn(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point, std::size_t k) = 0;
    virtual nb::dict statistics() const = 0;
};

template <std::size_t Dims, KDSplit SplitMethod>
struct KDTreeImpl final : KDTreeBase {
    using Tree = KDTree<Coord, Value, Dims, 128, SplitMethod>;
    using Point = typename Tree::Point;
    using Region = typename Tree::Region;

    std::optional<Tree> tree_{std::in_place};
    NextId next_id_ = 0;

    void bulk_load(const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& coords, const nb::object& ids = nb::none()) override {
        check_coords_shape(coords, Dims, "KDTree");
        const std::size_t n = static_cast<std::size_t>(coords.shape(0));
        const auto id_values = resolve_ids(ids, n, "KDTree");

        std::vector<typename Tree::PointPair> entries;
        entries.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            Point p;
            for (std::size_t d = 0; d < Dims; ++d)
                p[d] = coords(i, static_cast<int>(d));
            entries.emplace_back(p, id_values[i]);
        }
        tree_.emplace(entries.begin(), entries.end());
        update_next_id_after_bulk(ids, id_values, next_id_);
    }

    void insert(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point, const nb::object& id = nb::none()) override {
        tree_.value().insert(read_kd_point<Dims, Point>(point), resolve_insert_id(id, next_id_));
    }

    bool remove(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point) override {
        return tree_.value().remove(read_kd_point<Dims, Point>(point));
    }

    IdArray query(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& min_corner, const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& max_corner) override {
        check_corner_shape(min_corner, Dims, "KDTree.query min_corner");
        check_corner_shape(max_corner, Dims, "KDTree.query max_corner");
        Region region;
        for (std::size_t d = 0; d < Dims; ++d) {
            region.low[d] = min_corner(d);
            region.high[d] = max_corner(d);
        }
        std::vector<Value> out;
        tree_.value().query(region, std::back_inserter(out));
        return to_array(std::move(out));
    }

    IdArray knn(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point, std::size_t k) override {
        std::vector<Value> out;
        tree_.value().knnQuery(read_kd_point<Dims, Point>(point), k, std::back_inserter(out));
        return to_array(std::move(out));
    }

    nb::dict statistics() const override {
        const auto s = tree_.value().getStatistics();
        return stats_dict(s.height, s.numPoints, s.numLeaves, s.numInternalNodes, s.sizeBytes);
    }
};

template <std::size_t Dims>
std::unique_ptr<KDTreeBase> make_kdtree_impl(const std::string& split) {
    switch (parse_kd_split(split)) {
        case KDSplit::ADAPTIVE:
            return std::make_unique<KDTreeImpl<Dims, KDSplit::ADAPTIVE>>();
        case KDSplit::ROUND_ROBIN:
            return std::make_unique<KDTreeImpl<Dims, KDSplit::ROUND_ROBIN>>();
        case KDSplit::LONGEST_AXIS:
            return std::make_unique<KDTreeImpl<Dims, KDSplit::LONGEST_AXIS>>();
    }
    throw invalid_argument("KDTree: invalid split");
}

template <std::size_t Dims>
struct KDTreeIndex {
    std::string split_;
    std::unique_ptr<KDTreeBase> impl_;

    explicit KDTreeIndex(const std::string& split) : split_(split), impl_(make_kdtree_impl<Dims>(split)) {}

    const std::string& split() const { return split_; }

    void bulk_load(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& coords,
        const nb::object& ids = nb::none()
    ) {
        impl_->bulk_load(coords, ids);
    }

    void insert(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point, const nb::object& id = nb::none()) {
        impl_->insert(point, id);
    }

    bool remove(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point) {
        return impl_->remove(point);
    }

    IdArray query(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& min_corner,
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& max_corner
    ) {
        return impl_->query(min_corner, max_corner);
    }

    IdArray knn(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point,
        std::size_t k
    ) {
        return impl_->knn(point, k);
    }

    nb::dict statistics() const { return impl_->statistics(); }
};

}  // namespace

void bind_kdtree(nb::module_& m) {
    nb::class_<KDTreeIndex<2>>(m, "KDTree2D")
        .def(nb::init<const std::string&>(), nb::arg("split") = "longest_axis")
        .def_prop_ro("split", &KDTreeIndex<2>::split)
        .def("bulk_load", &KDTreeIndex<2>::bulk_load, nb::arg("coords"), nb::arg("ids") = nb::none())
        .def("insert", &KDTreeIndex<2>::insert, nb::arg("point"), nb::arg("id") = nb::none())
        .def("remove", &KDTreeIndex<2>::remove)
        .def("query", &KDTreeIndex<2>::query, nb::arg("min_corner"), nb::arg("max_corner"))
        .def("knn", &KDTreeIndex<2>::knn, nb::arg("point"), nb::arg("k"))
        .def("statistics", &KDTreeIndex<2>::statistics);

    nb::class_<KDTreeIndex<3>>(m, "KDTree3D")
        .def(nb::init<const std::string&>(), nb::arg("split") = "longest_axis")
        .def_prop_ro("split", &KDTreeIndex<3>::split)
        .def("bulk_load", &KDTreeIndex<3>::bulk_load, nb::arg("coords"), nb::arg("ids") = nb::none())
        .def("insert", &KDTreeIndex<3>::insert, nb::arg("point"), nb::arg("id") = nb::none())
        .def("remove", &KDTreeIndex<3>::remove)
        .def("query", &KDTreeIndex<3>::query, nb::arg("min_corner"), nb::arg("max_corner"))
        .def("knn", &KDTreeIndex<3>::knn, nb::arg("point"), nb::arg("k"))
        .def("statistics", &KDTreeIndex<3>::statistics);
}
