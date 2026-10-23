#include "common.hpp"
#include "quadtree.hpp"
#include <memory>
#include <string>
#include <utility>

namespace {

SplitStrategy parse_quad_split(const std::string& split) {
    if (split == "pr" || split == "PR")
        return SplitStrategy::PR;
    if (split == "pseudo_median" || split == "PSEUDO_MEDIAN")
        return SplitStrategy::PSEUDO_MEDIAN;
    if (split == "longest_axis" || split == "POINT_LONGEST_AXIS")
        return SplitStrategy::POINT_LONGEST_AXIS;
    throw std::invalid_argument("QuadTree: split must be one of 'pr', 'pseudo_median', 'longest_axis'");
}

struct QuadTreeBase {
    virtual ~QuadTreeBase() = default;
    virtual void bulk_load(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& coords,
        const nb::object& ids = nb::none()
    ) = 0;
    virtual void insert(Coord x, Coord y, const nb::object& id = nb::none()) = 0;
    virtual bool remove(Coord x, Coord y, Value id) = 0;
    virtual IdArray query(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& min_corner,
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& max_corner
    ) = 0;
    virtual IdArray knn(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point, std::size_t k) = 0;
    virtual nb::dict statistics() const = 0;
};

template <SplitStrategy Strategy>
struct QuadTreeImpl final : QuadTreeBase {
    using Node = QuadTreeNode<Coord, Value, 128, Strategy>;
    using Point = QPoint<Coord, Value>;
    using Box = QBox<Coord, Value>;

    std::unique_ptr<Node> tree_;
    NextId next_id_ = 0;

    void bulk_load(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& coords,
        const nb::object& ids = nb::none()
    ) override {
        check_coords_shape(coords, 2, "QuadTree");
        const std::size_t n = static_cast<std::size_t>(coords.shape(0));
        const auto id_values = resolve_ids(ids, n, "QuadTree");

        std::vector<Point> points;
        points.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            points.emplace_back(id_values[i], coords(i, 0), coords(i, 1));
        }
        tree_ = std::make_unique<Node>(points.begin(), points.end());
        update_next_id_after_bulk(ids, id_values, next_id_);
    }

    void insert(Coord x, Coord y, const nb::object& id = nb::none()) override {
        if (!tree_) {
            tree_ = std::make_unique<Node>();
        }
        tree_->insert(Point(resolve_insert_id(id, next_id_), x, y));
    }

    bool remove(Coord x, Coord y, Value id) override {
        if (!tree_) {
            return false;
        }
        return tree_->remove(Point(id, x, y));
    }

    IdArray query(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& min_corner,
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& max_corner
    ) override {
        check_corner_shape(min_corner, 2, "QuadTree.query min_corner");
        check_corner_shape(max_corner, 2, "QuadTree.query max_corner");
        if (!tree_) {
            return empty_array();
        }
        Box box(min_corner(0), min_corner(1), max_corner(0), max_corner(1), 0);
        std::vector<Value> out;
        tree_->query(box, std::back_inserter(out));
        return to_array(std::move(out));
    }

    IdArray knn(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point, std::size_t k) override {
        check_corner_shape(point, 2, "QuadTree.knn point");
        if (!tree_) {
            return empty_array();
        }
        std::array<Coord, 2> q{point(0), point(1)};
        std::vector<Value> out;
        tree_->knnQuery(q, k, std::back_inserter(out));
        return to_array(std::move(out));
    }

    nb::dict statistics() const override {
        if (!tree_) {
            return stats_dict(0, 0, 0, 0, 0, 0, "num_reroots");
        }
        const auto s = tree_->getStatistics();
        return stats_dict(s.height, s.numPoints, s.numLeaves, s.numInternalNodes, s.sizeBytes, s.numReroots, "num_reroots");
    }

};

struct QuadTreeIndex {
    std::string split_;
    std::unique_ptr<QuadTreeBase> impl_;

    explicit QuadTreeIndex(const std::string& split) : split_(split) {
        switch (parse_quad_split(split)) {
            case SplitStrategy::PR:
                impl_ = std::make_unique<QuadTreeImpl<SplitStrategy::PR>>();
                break;
            case SplitStrategy::PSEUDO_MEDIAN:
                impl_ = std::make_unique<QuadTreeImpl<SplitStrategy::PSEUDO_MEDIAN>>();
                break;
            case SplitStrategy::POINT_LONGEST_AXIS:
                impl_ = std::make_unique<QuadTreeImpl<SplitStrategy::POINT_LONGEST_AXIS>>();
                break;
        }
        split_ = split;
    }

    const std::string& split() const { return split_; }

    void bulk_load(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& coords,
        const nb::object& ids = nb::none()
    ) {
        impl_->bulk_load(coords, ids);
    }

    void insert(Coord x, Coord y, const nb::object& id = nb::none()) { impl_->insert(x, y, id); }
    bool remove(Coord x, Coord y, Value id) { return impl_->remove(x, y, id); }
    IdArray query(
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& min_corner,
        const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& max_corner
    ) {
        return impl_->query(min_corner, max_corner);
    }
    IdArray knn(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& point, std::size_t k) {
        return impl_->knn(point, k);
    }
    nb::dict statistics() const { return impl_->statistics(); }
};

}  // namespace

void bind_quadtree(nb::module_& m) {
    nb::class_<QuadTreeIndex>(m, "QuadTree")
        .def(nb::init<const std::string&>(), nb::arg("split") = "pr")
        .def_prop_ro("split", &QuadTreeIndex::split)
        .def("bulk_load", &QuadTreeIndex::bulk_load, nb::arg("coords"), nb::arg("ids") = nb::none())
        .def("insert", &QuadTreeIndex::insert, nb::arg("x"), nb::arg("y"), nb::arg("id") = nb::none())
        .def("remove", &QuadTreeIndex::remove)
        .def("query", &QuadTreeIndex::query, nb::arg("min_corner"), nb::arg("max_corner"))
        .def("knn", &QuadTreeIndex::knn, nb::arg("point"), nb::arg("k"))
        .def("statistics", &QuadTreeIndex::statistics);
}
