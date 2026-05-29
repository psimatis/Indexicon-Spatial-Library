#ifndef MXCIF_QUADTREE_HPP
#define MXCIF_QUADTREE_HPP

#include <algorithm>
#include <array>
#include <cfloat>
#include <limits>
#include <cstdint>
#include <queue>
#include <vector>

using namespace std;

constexpr size_t XLOW  = 0;
constexpr size_t YLOW  = 1;
constexpr size_t XHIGH = 2;
constexpr size_t YHIGH = 3;

constexpr size_t NW = 0;
constexpr size_t NE = 1;
constexpr size_t SW = 2;
constexpr size_t SE = 3;

struct MXCIFStats {
    size_t height;
    size_t numPoints;
    size_t numLeaves;
    size_t numInternalNodes;
    size_t sizeBytes;
};

template <typename CoordType, typename Value>
class MXCIFBox {
public:
    Value id;
    array<CoordType, 4> box;

    MXCIFBox(Value _id, CoordType xlow, CoordType ylow, CoordType xhigh, CoordType yhigh) {
        id = _id;
        box[XLOW] = xlow;
        box[YLOW] = ylow;
        box[XHIGH] = xhigh;
        box[YHIGH] = yhigh;
    }

    bool intersects(const MXCIFBox& r) const {
        return !(box[XLOW] > r.box[XHIGH] || r.box[XLOW] > box[XHIGH] ||
                 box[YLOW] > r.box[YHIGH] || r.box[YLOW] > box[YHIGH]);
    }

    array<CoordType, 2> toKNNPoint() const { return array<CoordType, 2>({box[XLOW], box[YLOW]}); }
};

template <typename Value>
struct MXCIFKnnPoint {
    double dist = FLT_MAX;
    Value id;
    bool operator<(const MXCIFKnnPoint<Value>& second) const { return dist < second.dist; }
};

template <typename CoordType, typename Value, size_t Capacity = 128>
class MXCIFQuadTreeNode {
public:
    array<MXCIFQuadTreeNode*, 4> children = {nullptr, nullptr, nullptr, nullptr};
    vector<MXCIFBox<CoordType, Value>> xStraddle; // crosses vertical midline, sorted by YLOW
    vector<MXCIFBox<CoordType, Value>> yStraddle; // crosses horizontal midline only, sorted by XLOW
    vector<MXCIFBox<CoordType, Value>> data; // leaf buffer for non-straddling items
    array<CoordType, 4> box;
    CoordType xMid, yMid;

    MXCIFQuadTreeNode() :
        box{numeric_limits<CoordType>::max(), numeric_limits<CoordType>::max(),
            numeric_limits<CoordType>::lowest(), numeric_limits<CoordType>::lowest()},
        xMid(0), yMid(0) {}

    MXCIFQuadTreeNode(const array<CoordType, 4>& _box) : 
        box(_box),
        xMid((_box[XLOW] + _box[XHIGH]) / 2),
        yMid((_box[YLOW] + _box[YHIGH]) / 2) {}

    template <typename Iter>
    MXCIFQuadTreeNode(Iter first, Iter last) {
        if (first == last) {
            box = {0, 0, 0, 0};
            xMid = yMid = 0;
            return;
        }
        box[XLOW] = first->box[XLOW];
        box[YLOW] = first->box[YLOW];
        box[XHIGH] = first->box[XHIGH];
        box[YHIGH] = first->box[YHIGH];
        for (auto it = first; it != last; ++it) {
            box[XLOW] = min(box[XLOW], it->box[XLOW]);
            box[YLOW] = min(box[YLOW], it->box[YLOW]);
            box[XHIGH] = max(box[XHIGH], it->box[XHIGH]);
            box[YHIGH] = max(box[YHIGH], it->box[YHIGH]);
        }
        updateMidpoints();
        pack(first, last);
    }

    ~MXCIFQuadTreeNode() {
        for (auto c : children)
            if (c) delete c;
    }

    bool isLeaf() const { 
        return children[NW] == nullptr; 
    }

    void updateMidpoints() {
        xMid = (box[XLOW] + box[XHIGH]) / 2;
        yMid = (box[YLOW] + box[YHIGH]) / 2;
    }

    bool crossesXMidline(const MXCIFBox<CoordType, Value>& r) const {
        return r.box[XLOW] < xMid && r.box[XHIGH] > xMid;
    }

    bool crossesYMidline(const MXCIFBox<CoordType, Value>& r) const {
        return r.box[YLOW] < yMid && r.box[YHIGH] > yMid;
    }

    bool storedHere(const MXCIFBox<CoordType, Value>& r) const {
        return crossesXMidline(r) || crossesYMidline(r);
    }

    void storeHere(const MXCIFBox<CoordType, Value>& r) {
        if (crossesXMidline(r))
            xStraddle.push_back(r);
        else
            yStraddle.push_back(r);
    }

    size_t getQuadrant(const MXCIFBox<CoordType, Value>& r) const {
        bool east = r.box[XLOW] >= xMid;
        bool north = r.box[YLOW] >= yMid;
        if (!east && north) 
            return NW;
        if (east && north) 
            return NE;
        if (!east && !north) 
            return SW;
        return SE;
    }

    bool intersectsBox(const MXCIFBox<CoordType, Value>& q) const {
        return !(box[XLOW] > q.box[XHIGH] || q.box[XLOW] > box[XHIGH] ||
                 box[YLOW] > q.box[YHIGH] || q.box[YLOW] > box[YHIGH]);
    }

    void rerootExpand(const MXCIFBox<CoordType, Value>& r) {
        CoordType w = box[XHIGH] - box[XLOW];
        CoordType h = box[YHIGH] - box[YLOW];
        if (w <= 0) 
            w = CoordType(1);
        if (h <= 0) 
            h = CoordType(1);

        CoordType newXLOW  = (r.box[XLOW]  < box[XLOW])  ? box[XLOW]  - w : box[XLOW];
        CoordType newXHIGH = (r.box[XHIGH] > box[XHIGH]) ? box[XHIGH] + w : box[XHIGH];
        CoordType newYLOW  = (r.box[YLOW]  < box[YLOW])  ? box[YLOW]  - h : box[YLOW];
        CoordType newYHIGH = (r.box[YHIGH] > box[YHIGH]) ? box[YHIGH] + h : box[YHIGH];

        CoordType newXMid = (newXLOW + newXHIGH) / 2;
        CoordType newYMid = (newYLOW + newYHIGH) / 2;

        CoordType oldCX = (box[XLOW] + box[XHIGH]) / 2;
        CoordType oldCY = (box[YLOW] + box[YHIGH]) / 2;
        size_t oldQ = (oldCX >= newXMid ? 1 : 0) + (oldCY < newYMid ? 2 : 0);

        auto* oldNode = new MXCIFQuadTreeNode<CoordType, Value, Capacity>(box);
        oldNode->xStraddle = move(xStraddle);
        oldNode->yStraddle = move(yStraddle);
        oldNode->data = move(data);
        oldNode->children = children;
        for (auto& c : children) c = nullptr;
        xStraddle.clear();
        yStraddle.clear();
        data.clear();

        box = {newXLOW, newYLOW, newXHIGH, newYHIGH};
        updateMidpoints();
        children[NW] = new MXCIFQuadTreeNode<CoordType, Value, Capacity>({newXLOW, newYMid, newXMid, newYHIGH});
        children[NE] = new MXCIFQuadTreeNode<CoordType, Value, Capacity>({newXMid, newYMid, newXHIGH, newYHIGH});
        children[SW] = new MXCIFQuadTreeNode<CoordType, Value, Capacity>({newXLOW, newYLOW, newXMid, newYMid});
        children[SE] = new MXCIFQuadTreeNode<CoordType, Value, Capacity>({newXMid, newYLOW, newXHIGH, newYMid});

        delete children[oldQ];
        children[oldQ] = oldNode;
    }

    void createChildren() {
        children[NW] = new MXCIFQuadTreeNode<CoordType, Value, Capacity>({box[XLOW], yMid, xMid, box[YHIGH]});
        children[NE] = new MXCIFQuadTreeNode<CoordType, Value, Capacity>({xMid, yMid, box[XHIGH], box[YHIGH]});
        children[SW] = new MXCIFQuadTreeNode<CoordType, Value, Capacity>({box[XLOW], box[YLOW], xMid, yMid});
        children[SE] = new MXCIFQuadTreeNode<CoordType, Value, Capacity>({xMid, box[YLOW], box[XHIGH], yMid});
    }

    template <typename Iter>
    void pack(Iter first, Iter last) {
        size_t count = distance(first, last);
        if (count == 0) 
            return;

        array<vector<MXCIFBox<CoordType, Value>>, 4> buckets;
        vector<MXCIFBox<CoordType, Value>> xs, ys;

        for (Iter it = first; it != last; ++it) {
            if (crossesXMidline(*it))       
                xs.push_back(*it);
            else if (crossesYMidline(*it))  
                ys.push_back(*it);
            else                            
                buckets[getQuadrant(*it)].push_back(*it);
        }

        sort(xs.begin(), xs.end(), [](const MXCIFBox<CoordType, Value>& a, const MXCIFBox<CoordType, Value>& b) {
            return a.box[YLOW] < b.box[YLOW];
        });
        sort(ys.begin(), ys.end(), [](const MXCIFBox<CoordType, Value>& a, const MXCIFBox<CoordType, Value>& b) {
            return a.box[XLOW] < b.box[XLOW];
        });
        xStraddle = move(xs);
        yStraddle = move(ys);

        size_t nonStraddle = 0;
        for (size_t i = 0; i < 4; ++i) 
            nonStraddle += buckets[i].size();

        if (nonStraddle == 0) 
            return;
        else if (nonStraddle <= Capacity) {
            for (size_t i = 0; i < 4; ++i) {
                for (auto& item : buckets[i])
                    data.push_back(item);
            }
        } else {
            createChildren();
            for (size_t i = 0; i < 4; ++i) {
                if (buckets[i].empty()) 
                    continue;
                if (buckets[i].size() < count)
                    children[i]->pack(buckets[i].begin(), buckets[i].end());
                else {
                    for (const auto& r : buckets[i]) 
                        children[i]->storeHere(r);
                }
            }
        }
    }

    void insert(const MXCIFBox<CoordType, Value>& r) {
        if (isLeaf()) {
            box[XLOW] = min(box[XLOW], r.box[XLOW]);
            box[YLOW] = min(box[YLOW], r.box[YLOW]);
            box[XHIGH] = max(box[XHIGH], r.box[XHIGH]);
            box[YHIGH] = max(box[YHIGH], r.box[YHIGH]);
            updateMidpoints();

            if (storedHere(r) || xMid == box[XLOW] || yMid == box[YLOW]) {
                storeHere(r);
                return;
            }

            data.push_back(r);
            if (data.size() <= Capacity)
                return;

            // Check if all non-straddling items land in the same quadrant
            bool canSplit = false;
            size_t firstQ = SIZE_MAX;
            for (const auto& item : data) {
                if (storedHere(item)) { 
                    canSplit = true; 
                    break; 
                }
                if (firstQ == SIZE_MAX) 
                    firstQ = getQuadrant(item);
                else if (getQuadrant(item) != firstQ) { 
                    canSplit = true; 
                    break; 
                }
            }
            if (!canSplit)
                return; // keep as oversized leaf

            createChildren();
            for (auto& item : data) {
                if (storedHere(item)) 
                    storeHere(item);
                else 
                    children[getQuadrant(item)]->insert(item);
            }
            data.clear();
            return;
        }

        while (r.box[XLOW] < box[XLOW] || r.box[XHIGH] > box[XHIGH] ||
               r.box[YLOW] < box[YLOW] || r.box[YHIGH] > box[YHIGH])
            rerootExpand(r);

        if (storedHere(r) || xMid == box[XLOW] || yMid == box[YLOW]) {
            storeHere(r);
            return;
        }
        children[getQuadrant(r)]->insert(r);
    }

    template <typename OutIter>
    void queryDumpAll(OutIter& out) {
        for (const auto& item : data) 
            *out++ = item.id;
        for (const auto& item : xStraddle) 
            *out++ = item.id;
        for (const auto& item : yStraddle) 
            *out++ = item.id;
        if (!isLeaf()) {
            for (auto c : children) 
                c->queryDumpAll(out);
        }
    }

    template <typename OutIter>
    void query(const MXCIFBox<CoordType, Value>& q, OutIter out) {
        queryImpl(q, out);
    }

    template <typename OutIter>
    void queryImpl(const MXCIFBox<CoordType, Value>& q, OutIter& out) {
        if (q.box[XLOW] <= box[XLOW] && q.box[XHIGH] >= box[XHIGH] &&
            q.box[YLOW] <= box[YLOW] && q.box[YHIGH] >= box[YHIGH]) {
            queryDumpAll(out);
            return;
        }

        for (const auto& item : xStraddle)
            if (q.intersects(item))
                *out++ = item.id;

        for (const auto& item : yStraddle)
            if (q.intersects(item))
                *out++ = item.id;

        if (isLeaf()) {
            for (const auto& item : data) {
                if (q.intersects(item))
                    *out++ = item.id;
            }
        } else {
            for (auto c : children) {
                if (c->intersectsBox(q))
                    c->queryImpl(q, out);
            }
        }
    }

    double minSqrDist(const array<CoordType, 4>& r) const {
        CoordType qx = r[XLOW], qy = r[YLOW];
        CoordType dx = max({CoordType(0), box[XLOW] - qx, qx - box[XHIGH]});
        CoordType dy = max({CoordType(0), box[YLOW] - qy, qy - box[YHIGH]});
        return dx * dx + dy * dy;
    }

    template <typename OutIter>
    void knnQuery(array<CoordType, 2> q, size_t k, OutIter out) {
        struct KnnNode {
            MXCIFQuadTreeNode<CoordType, Value, Capacity>* sn;
            double dist = FLT_MAX;
            bool operator<(const KnnNode& second) const { return dist > second.dist; }
        };

        auto mindistSqr = [&](const MXCIFBox<CoordType, Value>& box) {
            CoordType dx = max({CoordType(0), box.box[XLOW] - q[0], q[0] - box.box[XHIGH]});
            CoordType dy = max({CoordType(0), box.box[YLOW] - q[1], q[1] - box.box[YHIGH]});
            return dx * dx + dy * dy;
        };

        vector<MXCIFKnnPoint<Value>> tempPts(k);
        array<CoordType, 4> query{q[0], q[1], q[0], q[1]};
        priority_queue<MXCIFKnnPoint<Value>, vector<MXCIFKnnPoint<Value>>> knnPts(tempPts.begin(), tempPts.end());
        priority_queue<KnnNode, vector<KnnNode>> unseenNodes;
        unseenNodes.emplace(KnnNode{this, minSqrDist(query)});
        double d, minDist;
        MXCIFQuadTreeNode<CoordType, Value, Capacity>* node;

        while (!unseenNodes.empty()) {
            node = unseenNodes.top().sn;
            d = unseenNodes.top().dist;
            unseenNodes.pop();
            minDist = knnPts.top().dist;
            
            if (d >= minDist) 
                break;

            for (const auto& p : node->xStraddle) {
                minDist = knnPts.top().dist;
                d = mindistSqr(p);
                if (d < minDist) { 
                    knnPts.pop(); 
                    knnPts.push({d, p.id}); 
                }
            }
            for (const auto& p : node->yStraddle) {
                minDist = knnPts.top().dist;
                d = mindistSqr(p);
                if (d < minDist) { 
                    knnPts.pop(); 
                    knnPts.push({d, p.id}); 
                }
            }
            if (node->isLeaf()) {
                for (const auto& p : node->data) {
                    minDist = knnPts.top().dist;
                    d = mindistSqr(p);
                    if (d < minDist) { 
                        knnPts.pop(); 
                        knnPts.push({d, p.id}); 
                    }
                }
            } else {
                minDist = knnPts.top().dist;
                for (auto c : node->children) {
                    d = c->minSqrDist(query);
                    if (d < minDist) 
                        unseenNodes.push({c, d});
                }
            }
        }
        while (!knnPts.empty()) {
            *out++ = knnPts.top().id;
            knnPts.pop();
        }
    }

    bool remove(const MXCIFBox<CoordType, Value>& r) {
        auto& xs = xStraddle;
        auto it = find_if(xs.begin(), xs.end(), [&r](const MXCIFBox<CoordType, Value>& e) { return e.id == r.id; });
        if (it != xs.end()) { 
            xs.erase(it); 
            return true; 
        }

        auto& ys = yStraddle;
        it = find_if(ys.begin(), ys.end(), [&r](const MXCIFBox<CoordType, Value>& e) { return e.id == r.id; });
        if (it != ys.end()) { 
            ys.erase(it); 
            return true; 
        }

        if (isLeaf()) {
            auto dit = find_if(data.begin(), data.end(), [&r](const MXCIFBox<CoordType, Value>& e) { return e.id == r.id; });
            if (dit != data.end()) { 
                data.erase(dit); 
                return true; 
            }
            return false;
        }
        if (!storedHere(r)) {
            if (children[getQuadrant(r)]->remove(r)) {
                tryMerge();
                return true;
            }
        }
        return false;
    }

    void tryMerge() {
        for (auto c : children)
            if (!c->isLeaf()) 
                return;

        size_t total = 0;
        for (auto c : children)
            total += c->data.size() + c->xStraddle.size() + c->yStraddle.size();

        if (total <= Capacity) {
            for (auto c : children) {
                for (auto& item : c->data) data.push_back(move(item));
                for (auto& item : c->xStraddle) data.push_back(move(item));
                for (auto& item : c->yStraddle) data.push_back(move(item));
                delete c;
            }
            children = {nullptr, nullptr, nullptr, nullptr};
        }
    }

    MXCIFStats getStatistics() const {
        MXCIFStats stats = {0, 0, 0, 0, 0};
        stats.sizeBytes += sizeof(*this);
        stats.numPoints += xStraddle.size() + yStraddle.size() + data.size();
        stats.sizeBytes += (xStraddle.capacity() + yStraddle.capacity() + data.capacity()) * sizeof(MXCIFBox<CoordType, Value>);

        if (isLeaf()) 
            stats.numLeaves++;
        else {
            stats.numInternalNodes++;
            for (auto c : children) {
                MXCIFStats cs = c->getStatistics();
                stats.height = max(stats.height, cs.height + 1);
                stats.numPoints += cs.numPoints;
                stats.numLeaves += cs.numLeaves;
                stats.numInternalNodes += cs.numInternalNodes;
                stats.sizeBytes += cs.sizeBytes;
            }
        }
        return stats;
    }
};

#endif // MXCIF_QUADTREE_HPP
