#ifndef QUADTREE_HPP
#define QUADTREE_HPP

#include <algorithm>
#include <array>
#include <cfloat>
#include <cstdint>
#include <queue>
#include <vector>

using namespace std;

#ifndef XLOW_DEFINED
#define XLOW_DEFINED
constexpr size_t XLOW = 0;
constexpr size_t YLOW = 1;
constexpr size_t XHIGH = 2;
constexpr size_t YHIGH = 3;

constexpr size_t NW = 0;
constexpr size_t NE = 1;
constexpr size_t SW = 2;
constexpr size_t SE = 3;
#endif

struct QuadTreeStats {
    size_t height;
    size_t numPoints;
    size_t numLeaves;
    size_t numInternalNodes;
    size_t numReroots;
    size_t sizeBytes;
};

template <typename CoordType, typename Value>
class QPoint {
public:
    Value id;
    array<CoordType, 2> pt;

    QPoint(Value _id, CoordType x, CoordType y) {
        id = _id;
        pt[0] = x;
        pt[1] = y;
    }
};

template <typename CoordType, typename Value>
class QBox {
public:
    Value id;
    QPoint<CoordType, Value> lo;
    QPoint<CoordType, Value> hi;

    QBox(CoordType xlow, CoordType ylow, CoordType xhigh, CoordType yhigh, Value _id)
        : id(_id), lo(_id, xlow, ylow), hi(_id, xhigh, yhigh) {}

    bool intersects(const QBox<CoordType, Value>& r) const {
        return !(lo.pt[0] > r.hi.pt[0] || hi.pt[0] < r.lo.pt[0] ||
                 lo.pt[1] > r.hi.pt[1] || hi.pt[1] < r.lo.pt[1]);
    }
};

enum class SplitStrategy {
    PR,
    PSEUDO_MEDIAN,
    POINT_LONGEST_AXIS
};

template <typename CoordType, typename Value, size_t Capacity = 128, SplitStrategy Strategy = SplitStrategy::PR>
class QuadTreeNode {
public:
    QuadTreeNode(const array<CoordType, 4>& _boundary) : box(_boundary) {}

    array<QuadTreeNode<CoordType, Value, Capacity, Strategy>*, 4> children = {nullptr, nullptr, nullptr, nullptr};
    vector<QPoint<CoordType, Value>> data;
    array<CoordType, 4> box;
    size_t rerootCount = 0;

    QuadTreeNode() : box{0, 0, 1, 1} {}

    template <typename Iter>
    QuadTreeNode(Iter first, Iter last, size_t maxDepth = SIZE_MAX) {
        if (first == last) {
            box = {0, 0, 0, 0};
            return;
        }
        CoordType minX = first->pt[0], minY = first->pt[1];
        CoordType maxX = first->pt[0], maxY = first->pt[1];
        for (auto it = first; it != last; ++it) {
            minX = min(minX, it->pt[0]);
            minY = min(minY, it->pt[1]);
            maxX = max(maxX, it->pt[0]);
            maxY = max(maxY, it->pt[1]);
        }
        box = {minX, minY, maxX, maxY};
        pack(first, last, 0, maxDepth);
    }

    ~QuadTreeNode() {
        for (auto c : children) {
            delete c;
        }
    }

    bool isLeaf() const { return children[NW] == nullptr; }

    template <typename Iter>
    void computeSplit(Iter first, Iter last, CoordType& xMid, CoordType& yMid) {
        if constexpr (Strategy == SplitStrategy::PR) {
            xMid = (box[XLOW] + box[XHIGH]) / 2;
            yMid = (box[YLOW] + box[YHIGH]) / 2;
        } else if constexpr (Strategy == SplitStrategy::PSEUDO_MEDIAN) {
            size_t count = distance(first, last);
            size_t midIdx = count / 2;

            nth_element(first, first + midIdx, last, [](const QPoint<CoordType, Value>& a, const QPoint<CoordType, Value>& b) {
                return a.pt[0] < b.pt[0];
            });
            xMid = (first + midIdx)->pt[0];

            nth_element(first, first + midIdx, last, [](const QPoint<CoordType, Value>& a, const QPoint<CoordType, Value>& b) {
                return a.pt[1] < b.pt[1];
            });
            yMid = (first + midIdx)->pt[1];
        } else if constexpr (Strategy == SplitStrategy::POINT_LONGEST_AXIS) {
            size_t count = distance(first, last);
            size_t midIdx = count / 2;
            
            CoordType width = box[XHIGH] - box[XLOW];
            CoordType height = box[YHIGH] - box[YLOW];
            
            if (width >= height) {
                nth_element(first, first + midIdx, last, [](const QPoint<CoordType, Value>& a, const QPoint<CoordType, Value>& b) {
                    return a.pt[0] < b.pt[0];
                });
                xMid = (first + midIdx)->pt[0];
                yMid = (first + midIdx)->pt[1];
            } else {
                nth_element(first, first + midIdx, last, [](const QPoint<CoordType, Value>& a, const QPoint<CoordType, Value>& b) {
                    return a.pt[1] < b.pt[1];
                });
                xMid = (first + midIdx)->pt[0];
                yMid = (first + midIdx)->pt[1];
            }
        }
    }

    size_t getQuadrant(const QPoint<CoordType, Value>& p, CoordType xMid, CoordType yMid) const {
        if (p.pt[0] < xMid) {
            if (p.pt[1] >= yMid) 
                return NW;
            else 
                return SW;
        } else {
            if (p.pt[1] >= yMid) 
                return NE;
            else 
                return SE;
        }
    }

    template <typename Iter>
    void pack(Iter first, Iter last, size_t depth = 0, size_t maxDepth = SIZE_MAX) {
        size_t count = distance(first, last);
        if (count <= Capacity || depth >= maxDepth) {
            data.assign(first, last);
            return;
        }

        CoordType xMid, yMid;
        computeSplit(first, last, xMid, yMid);

        children[NW] = new QuadTreeNode<CoordType, Value, Capacity, Strategy>(array<CoordType, 4>{box[XLOW], yMid, xMid, box[YHIGH]});
        children[NE] = new QuadTreeNode<CoordType, Value, Capacity, Strategy>(array<CoordType, 4>{xMid, yMid, box[XHIGH], box[YHIGH]});
        children[SW] = new QuadTreeNode<CoordType, Value, Capacity, Strategy>(array<CoordType, 4>{box[XLOW], box[YLOW], xMid, yMid});
        children[SE] = new QuadTreeNode<CoordType, Value, Capacity, Strategy>(array<CoordType, 4>{xMid, box[YLOW], box[XHIGH], yMid});

        array<vector<QPoint<CoordType, Value>>, 4> buckets;
        for (Iter it = first; it != last; ++it) {
            size_t q = getQuadrant(*it, xMid, yMid);
            buckets[q].push_back(*it);
        }

        bool allIdentical = true;
        for (Iter it = first + 1; it != last && allIdentical; ++it) {
            if (it->pt[0] != first->pt[0] || it->pt[1] != first->pt[1])
                allIdentical = false;
        }

        for (size_t i = 0; i < 4; ++i) {
            if (buckets[i].size() < count && !allIdentical)
                children[i]->pack(buckets[i].begin(), buckets[i].end(), depth + 1, maxDepth);
            else // Allow oversized leaf to handle duplicates
                children[i]->data = move(buckets[i]);
        }
    }

    void rerootExpand(const QPoint<CoordType, Value>& r) {
        ++rerootCount;
        CoordType w = box[XHIGH] - box[XLOW];
        CoordType h = box[YHIGH] - box[YLOW];
        if (w <= 0) 
            w = CoordType(1);
        if (h <= 0) 
            h = CoordType(1);

        CoordType newXLOW = (r.pt[0] < box[XLOW]) ? box[XLOW] - w : box[XLOW];
        CoordType newXHIGH = (r.pt[0] > box[XHIGH]) ? box[XHIGH] + w : box[XHIGH];
        CoordType newYLOW = (r.pt[1] < box[YLOW]) ? box[YLOW] - h : box[YLOW];
        CoordType newYHIGH = (r.pt[1] > box[YHIGH]) ? box[YHIGH] + h : box[YHIGH];

        CoordType oldCenterX = (box[XLOW] + box[XHIGH]) / 2;
        CoordType oldCenterY = (box[YLOW] + box[YHIGH]) / 2;

        auto* oldNode = new QuadTreeNode<CoordType, Value, Capacity, Strategy>(box);
        oldNode->data = move(data);
        oldNode->children = children;
        data.clear();

        box = {newXLOW, newYLOW, newXHIGH, newYHIGH};
        
        CoordType xMid, yMid;
        if constexpr (Strategy == SplitStrategy::PR) {
            xMid = (newXLOW + newXHIGH) / 2;
            yMid = (newYLOW + newYHIGH) / 2;
        } else {
            xMid = (box[XLOW] == oldNode->box[XLOW]) ? oldNode->box[XHIGH] : oldNode->box[XLOW];
            yMid = (box[YLOW] == oldNode->box[YLOW]) ? oldNode->box[YHIGH] : oldNode->box[YLOW];
            
            if (xMid == box[XLOW] || xMid == box[XHIGH]) 
                xMid = (newXLOW + newXHIGH) / 2;
            if (yMid == box[YLOW] || yMid == box[YHIGH]) 
                yMid = (newYLOW + newYHIGH) / 2;
        }

        size_t oldQ = (oldCenterX >= xMid ? 1 : 0) + (oldCenterY < yMid ? 2 : 0);

        children[NW] = new QuadTreeNode<CoordType, Value, Capacity, Strategy>({newXLOW, yMid, xMid, newYHIGH});
        children[NE] = new QuadTreeNode<CoordType, Value, Capacity, Strategy>({xMid, yMid, newXHIGH, newYHIGH});
        children[SW] = new QuadTreeNode<CoordType, Value, Capacity, Strategy>({newXLOW, newYLOW, xMid, yMid});
        children[SE] = new QuadTreeNode<CoordType, Value, Capacity, Strategy>({xMid, newYLOW, newXHIGH, yMid});

        delete children[oldQ];
        children[oldQ] = oldNode;
    }

    void insert(const QPoint<CoordType, Value>& r, size_t depth = 0, size_t maxDepth = SIZE_MAX) {
        if (isLeaf()) {
            data.push_back(r);
            if (data.size() > Capacity && depth < maxDepth) {
                bool allIdentical = true;
                for (size_t i = 1; i < data.size() && allIdentical; ++i) {
                    if (data[i].pt[0] != data[0].pt[0] || data[i].pt[1] != data[0].pt[1])
                        allIdentical = false;
                }
                if (!allIdentical) {
                    for (const auto& rec : data) {
                        box[XLOW] = min(box[XLOW], rec.pt[0]);
                        box[YLOW] = min(box[YLOW], rec.pt[1]);
                        box[XHIGH] = max(box[XHIGH], rec.pt[0]);
                        box[YHIGH] = max(box[YHIGH], rec.pt[1]);
                    }
                    divide();
                    CoordType xMid = children[NE]->box[XLOW];
                    CoordType yMid = children[NE]->box[YLOW];

                    for (const auto& rec : data) {
                        size_t q = getQuadrant(rec, xMid, yMid);
                        children[q]->data.push_back(rec);
                    }
                    data.clear();
                }
                // else: keep as oversized leaf
            }
        } else {
            while (r.pt[0] < box[XLOW] || r.pt[0] > box[XHIGH] ||
                   r.pt[1] < box[YLOW] || r.pt[1] > box[YHIGH])
                rerootExpand(r);
            
            CoordType xMid = children[NE]->box[XLOW];
            CoordType yMid = children[NE]->box[YLOW];
            size_t q = getQuadrant(r, xMid, yMid);
            children[q]->insert(r, depth + 1, maxDepth);
        }
    }

    void divide() {
        CoordType xMid, yMid;
        computeSplit(data.begin(), data.end(), xMid, yMid);

        children[NW] = new QuadTreeNode<CoordType, Value, Capacity, Strategy>(array<CoordType, 4>{box[XLOW], yMid, xMid, box[YHIGH]});
        children[NE] = new QuadTreeNode<CoordType, Value, Capacity, Strategy>(array<CoordType, 4>{xMid, yMid, box[XHIGH], box[YHIGH]});
        children[SW] = new QuadTreeNode<CoordType, Value, Capacity, Strategy>(array<CoordType, 4>{box[XLOW], box[YLOW], xMid, yMid});
        children[SE] = new QuadTreeNode<CoordType, Value, Capacity, Strategy>(array<CoordType, 4>{xMid, box[YLOW], box[XHIGH], yMid});
    }

    bool intersects(const QBox<CoordType, Value>& q) {
        return !(box[XLOW] > q.hi.pt[0] || q.lo.pt[0] > box[XHIGH] || box[YLOW] > q.hi.pt[1] || q.lo.pt[1] > box[YHIGH]);
    }

    template <typename OutIter>
    void queryDumpAll(OutIter& out) const {
        if (isLeaf()) {
            for (const auto& r : data)
                *out++ = r.id;
        } else {
            for (auto c : children)
                c->queryDumpAll(out);
        }
    }

    template <typename OutIter>
    void query(const QBox<CoordType, Value>& q, OutIter out) {
        queryImpl(q, out);
    }

    template <typename OutIter>
    void queryImpl(const QBox<CoordType, Value>& q, OutIter& out) {
        if (q.lo.pt[0] <= box[XLOW] && q.hi.pt[0] >= box[XHIGH] &&
            q.lo.pt[1] <= box[YLOW] && q.hi.pt[1] >= box[YHIGH]) {
            queryDumpAll(out);
            return;
        }

        if (isLeaf()) {
            for (const auto& r : data) {
                if (r.pt[0] >= q.lo.pt[0] && r.pt[0] <= q.hi.pt[0] &&
                    r.pt[1] >= q.lo.pt[1] && r.pt[1] <= q.hi.pt[1])
                    *out++ = r.id;
            }
        } else {
            for (auto c : children) {
                if (c->intersects(q))
                    c->queryImpl(q, out);
            }
        }
    }

    double minSqrDist(const array<CoordType, 4>& r) const {
        CoordType dx = max({CoordType(0), box[XLOW] - r[2], r[0] - box[XHIGH]});
        CoordType dy = max({CoordType(0), box[YLOW] - r[3], r[1] - box[YHIGH]});
        return dx * dx + dy * dy;
    }

    template <typename OutIter>
    void knnQuery(array<CoordType, 2> q, size_t k, OutIter out) {
        struct Neighbor {
            double dist;
            Value id;
        };
        auto neighbor_cmp = [](const Neighbor& a, const Neighbor& b) { return a.dist < b.dist; };
        vector<Neighbor> neighbors;
        neighbors.reserve(k);

        struct KnnNode {
            double dist;
            QuadTreeNode<CoordType, Value, Capacity, Strategy>* sn;
            bool operator>(const KnnNode& o) const { return dist > o.dist; }
        };
        priority_queue<KnnNode, vector<KnnNode>, greater<KnnNode>> unseenNodes;

        auto ignoreNode = [&](double d) { return neighbors.size() == k && neighbors.front().dist <= d; };

        auto calcSqrDist = [](const array<CoordType, 2>& x, const array<CoordType, 2>& y) {
            CoordType dx = x[0] - y[0], dy = x[1] - y[1];
            return (double)(dx * dx + dy * dy);
        };

        array<CoordType, 4> query{q[0], q[1], q[0], q[1]};
        KnnNode start; start.dist = minSqrDist(query); start.sn = this;
        unseenNodes.push(start);

        while (!unseenNodes.empty()) {
            double d = unseenNodes.top().dist;
            auto* node = unseenNodes.top().sn;
            unseenNodes.pop();

            if (ignoreNode(d))
                break;

            if (node->isLeaf()) {
                for (const auto& p : node->data) {
                    double dd = calcSqrDist(p.pt, q);
                    if (neighbors.size() < k) {
                        Neighbor nb; nb.dist = dd; nb.id = p.id;
                        neighbors.push_back(nb);
                        if (neighbors.size() == k)
                            make_heap(neighbors.begin(), neighbors.end(), neighbor_cmp);
                    } else if (dd < neighbors.front().dist) {
                        pop_heap(neighbors.begin(), neighbors.end(), neighbor_cmp);
                        neighbors.back().dist = dd;
                        neighbors.back().id = p.id;
                        push_heap(neighbors.begin(), neighbors.end(), neighbor_cmp);
                    }
                }
            } else {
                for (auto c : node->children) {
                    double cd = c->minSqrDist(query);
                    if (!ignoreNode(cd)) {
                        KnnNode kn; kn.dist = cd; kn.sn = c;
                        unseenNodes.push(kn);
                    }
                }
            }
        }

        for (const auto& n : neighbors)
            *out++ = n.id;
    }

    bool remove(const QPoint<CoordType, Value>& r) {
        if (isLeaf()) {
            auto it = find_if(data.begin(), data.end(), [&r](const QPoint<CoordType, Value>& e) { return e.id == r.id; });
            if (it == data.end()) 
                return false;
            data.erase(it);
            return true;
        } else {
            CoordType xMid = children[NE]->box[XLOW];
            CoordType yMid = children[NE]->box[YLOW];
            size_t q = getQuadrant(r, xMid, yMid);
            if (children[q]->remove(r)) {
                tryMerge();
                return true;
            }
            return false;
        }
    }

    void tryMerge() {
        for (auto c : children)
            if (!c->isLeaf()) 
                return;

        size_t total = 0;
        for (auto c : children)
            total += c->data.size();

        if (total <= Capacity) {
            for (auto c : children) {
                data.insert(data.end(), c->data.begin(), c->data.end());
                delete c;
            }
            children = {nullptr, nullptr, nullptr, nullptr};
        }
    }

    QuadTreeStats getStatistics() const {
        QuadTreeStats stats = {0, 0, 0, 0, 0, 0};
        stats.numReroots = rerootCount;
        stats.sizeBytes += sizeof(*this);
        stats.sizeBytes += data.capacity() * sizeof(QPoint<CoordType, Value>);
        if (isLeaf()) {
            stats.numLeaves++;
            stats.numPoints += data.size();
        } else {
            stats.numInternalNodes++;
            for (auto c : children) {
                QuadTreeStats childStats = c->getStatistics();
                stats.height = max(stats.height, childStats.height);
                stats.numPoints += childStats.numPoints;
                stats.numLeaves += childStats.numLeaves;
                stats.numInternalNodes += childStats.numInternalNodes;
                stats.numReroots += childStats.numReroots;
                stats.sizeBytes += childStats.sizeBytes;
            }
            stats.height += 1;
        }
        return stats;
    }

};

#endif // QUADTREE_HPP
