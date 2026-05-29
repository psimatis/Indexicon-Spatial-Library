#ifndef OCTREE_HPP
#define OCTREE_HPP

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <queue>
#include <vector>

using namespace std;

constexpr size_t XLOW = 0, YLOW = 1, ZLOW = 2, XHIGH = 3, YHIGH = 4, ZHIGH = 5;

enum Octant { LLL=0, HLL=1, LHL=2, HHL=3, LLH=4, HLH=5, LHH=6, HHH=7 };

struct OctreeStats {
    size_t height;
    size_t numPoints;
    size_t numLeaves;
    size_t numInternalNodes;
    size_t sizeBytes;
};

template <typename CoordType, typename Value>
struct OPoint {
    Value id;
    CoordType x, y, z;

    OPoint(Value _id, CoordType _x, CoordType _y, CoordType _z) : id(_id), x(_x), y(_y), z(_z) {}

    bool intersects(const array<CoordType, 6>& b) const {
        return x >= b[XLOW] && x <= b[XHIGH] &&
               y >= b[YLOW] && y <= b[YHIGH] &&
               z >= b[ZLOW] && z <= b[ZHIGH];
    }

    array<CoordType, 3> toKNNPoint() const { return array<CoordType, 3>({x, y, z}); }
};

template <typename CoordType, typename Value>
struct OBox {
    Value id;
    array<CoordType, 6> box; 

    OBox(Value _id, CoordType x, CoordType y, CoordType z) {
        id = _id;
        box[XLOW] = box[XHIGH] = x; box[YLOW] = box[YHIGH] = y; box[ZLOW] = box[ZHIGH] = z;
    }

    OBox(CoordType xlow, CoordType ylow, CoordType zlow, CoordType xhigh, CoordType yhigh, CoordType zhigh)
        : id(Value{}), box{xlow, ylow, zlow, xhigh, yhigh, zhigh} {}

    bool intersects(const array<CoordType, 6>& b) const {
        return !(box[XLOW] > b[XHIGH] || box[XHIGH] < b[XLOW] ||
                 box[YLOW] > b[YHIGH] || box[YHIGH] < b[YLOW] ||
                 box[ZLOW] > b[ZHIGH] || box[ZHIGH] < b[ZLOW]);
    }

    array<CoordType, 3> toKNNPoint() const { return array<CoordType, 3>({box[XLOW], box[YLOW], box[ZLOW]}); }
};

template <typename CoordType, typename Value>
struct OKnnPoint {
    double dist = FLT_MAX;
    Value id;
    bool operator<(const OKnnPoint<CoordType, Value>& second) const { return dist < second.dist; }
};

template <typename CoordType, typename Value, size_t Capacity = 128>
class OctreeNode {
public:
    array<OctreeNode*, 8> children = {nullptr};
    vector<OPoint<CoordType, Value>> data;
    array<CoordType, 6> bounds;

    OctreeNode(array<CoordType, 6> _b) : bounds(_b) {}

    OctreeNode() {
        bounds = {0, 0, 0, 1, 1, 1};
    }

    template <typename Iter>
    OctreeNode(Iter first, Iter last) {
        if (first == last) {
            bounds = {0, 0, 0, 0, 0, 0};
            return;
        }
        
        CoordType minX = first->x, minY = first->y, minZ = first->z;
        CoordType maxX = first->x, maxY = first->y, maxZ = first->z;
        
        for (auto it = first; it != last; ++it) {
            minX = min(minX, it->x);
            minY = min(minY, it->y);
            minZ = min(minZ, it->z);
            maxX = max(maxX, it->x);
            maxY = max(maxY, it->y);
            maxZ = max(maxZ, it->z);
        }
        
        bounds = {minX, minY, minZ, maxX, maxY, maxZ};
        pack(first, last);
    }

    ~OctreeNode() { 
        for (auto c : children) 
            delete c;
    }

    bool isLeaf() const { return children[0] == nullptr; }

    void divide() {
        CoordType xm = (bounds[XLOW] + bounds[XHIGH]) / 2;
        CoordType ym = (bounds[YLOW] + bounds[YHIGH]) / 2;
        CoordType zm = (bounds[ZLOW] + bounds[ZHIGH]) / 2;

        children[LLL] = new OctreeNode({bounds[XLOW], bounds[YLOW], bounds[ZLOW], xm, ym, zm});
        children[HLL] = new OctreeNode({xm, bounds[YLOW], bounds[ZLOW], bounds[XHIGH], ym, zm});
        children[LHL] = new OctreeNode({bounds[XLOW], ym, bounds[ZLOW], xm, bounds[YHIGH], zm});
        children[HHL] = new OctreeNode({xm, ym, bounds[ZLOW], bounds[XHIGH], bounds[YHIGH], zm});
        children[LLH] = new OctreeNode({bounds[XLOW], bounds[YLOW], zm, xm, ym, bounds[ZHIGH]});
        children[HLH] = new OctreeNode({xm, bounds[YLOW], zm, bounds[XHIGH], ym, bounds[ZHIGH]});
        children[LHH] = new OctreeNode({bounds[XLOW], ym, zm, xm, bounds[YHIGH], bounds[ZHIGH]});
        children[HHH] = new OctreeNode({xm, ym, zm, bounds[XHIGH], bounds[YHIGH], bounds[ZHIGH]});
    }

    template <typename Iter>
    void pack(Iter first, Iter last) {
        size_t count = distance(first, last);
        if (count <= Capacity) {
            data.assign(first, last);
            return;
        }

        divide();
        array<vector<OPoint<CoordType, Value>>, 8> buckets;
        for (auto it = first; it != last; ++it) {
            for (size_t i = 0; i < 8; ++i) {
                if (it->intersects(children[i]->bounds)) {
                    buckets[i].push_back(*it);
                    break;
                }
            }
        }

        for (size_t i = 0; i < 8; ++i) {
            if (!buckets[i].empty()) {
                if (buckets[i].size() < count) 
                    children[i]->pack(buckets[i].begin(), buckets[i].end());
                else 
                    children[i]->data = move(buckets[i]);
            }
        }
    }

    void insert(const OPoint<CoordType, Value>& p) {
        if (isLeaf()) {
            data.push_back(p);
            if (data.size() > Capacity) {
                bool allIdentical = true;
                for(const auto& it : data) {
                    if(it.x != data[0].x || it.y != data[0].y || it.z != data[0].z) {
                        allIdentical = false; 
                        break;
                    }
                }
                if (!allIdentical) {
                    vector<OPoint<CoordType, Value>> temp = move(data);
                    divide();
                    for (const auto& item : temp) 
                        insert(item);
                }
            }
        } else {
            while (p.x < bounds[XLOW] || p.x > bounds[XHIGH] ||
                   p.y < bounds[YLOW] || p.y > bounds[YHIGH] ||
                   p.z < bounds[ZLOW] || p.z > bounds[ZHIGH])
                rerootExpand(p);
            for (auto c : children) {
                if (p.intersects(c->bounds)) { 
                    c->insert(p); 
                    return; 
                }
            }
        }
    }

    void rerootExpand(const OPoint<CoordType, Value>& p) {
        CoordType w = bounds[XHIGH] - bounds[XLOW];
        CoordType h = bounds[YHIGH] - bounds[YLOW];
        CoordType d = bounds[ZHIGH] - bounds[ZLOW];
        if (w <= 0) 
            w = CoordType(1);
        if (h <= 0) 
            h = CoordType(1);
        if (d <= 0) 
            d = CoordType(1);

        CoordType newXLOW  = (p.x < bounds[XLOW])  ? bounds[XLOW]  - w : bounds[XLOW];
        CoordType newXHIGH = (p.x > bounds[XHIGH]) ? bounds[XHIGH] + w : bounds[XHIGH];
        CoordType newYLOW  = (p.y < bounds[YLOW])  ? bounds[YLOW]  - h : bounds[YLOW];
        CoordType newYHIGH = (p.y > bounds[YHIGH]) ? bounds[YHIGH] + h : bounds[YHIGH];
        CoordType newZLOW  = (p.z < bounds[ZLOW])  ? bounds[ZLOW]  - d : bounds[ZLOW];
        CoordType newZHIGH = (p.z > bounds[ZHIGH]) ? bounds[ZHIGH] + d : bounds[ZHIGH];

        CoordType xm = (newXLOW + newXHIGH) / 2;
        CoordType ym = (newYLOW + newYHIGH) / 2;
        CoordType zm = (newZLOW + newZHIGH) / 2;

        CoordType oldCX = (bounds[XLOW] + bounds[XHIGH]) / 2;
        CoordType oldCY = (bounds[YLOW] + bounds[YHIGH]) / 2;
        CoordType oldCZ = (bounds[ZLOW] + bounds[ZHIGH]) / 2;
        
        size_t oldOctant = 0;
        if (oldCX >= xm) 
            oldOctant |= 1; 
        if (oldCY >= ym) 
            oldOctant |= 2; 
        if (oldCZ >= zm) 
            oldOctant |= 4; 

        auto* oldNode = new OctreeNode(bounds);
        oldNode->data = move(data);
        oldNode->children = children;
        for (auto& c : children) 
            c = nullptr;
        data.clear();

        bounds = {newXLOW, newYLOW, newZLOW, newXHIGH, newYHIGH, newZHIGH};
        children[LLL] = new OctreeNode({newXLOW, newYLOW, newZLOW, xm, ym, zm});
        children[HLL] = new OctreeNode({xm, newYLOW, newZLOW, newXHIGH, ym, zm});
        children[LHL] = new OctreeNode({newXLOW, ym, newZLOW, xm, newYHIGH, zm});
        children[HHL] = new OctreeNode({xm, ym, newZLOW, newXHIGH, newYHIGH, zm});
        children[LLH] = new OctreeNode({newXLOW, newYLOW, zm, xm, ym, newZHIGH});
        children[HLH] = new OctreeNode({xm, newYLOW, zm, newXHIGH, ym, newZHIGH});
        children[LHH] = new OctreeNode({newXLOW, ym, zm, xm, newYHIGH, newZHIGH});
        children[HHH] = new OctreeNode({xm, ym, zm, newXHIGH, newYHIGH, newZHIGH});

        delete children[oldOctant];
        children[oldOctant] = oldNode;
    }

    double minSqrDist(const array<CoordType, 6>& query) const {
        double dist = 0;
        if (query[XLOW] < bounds[XLOW]) {
            CoordType d = bounds[XLOW] - query[XHIGH];
            dist += d * d;
        } else if (query[XLOW] > bounds[XHIGH]) {
            CoordType d = query[XLOW] - bounds[XHIGH];
            dist += d * d;
        }
        if (query[YLOW] < bounds[YLOW]) {
            CoordType d = bounds[YLOW] - query[YHIGH];
            dist += d * d;
        } else if (query[YLOW] > bounds[YHIGH]) {
            CoordType d = query[YLOW] - bounds[YHIGH];
            dist += d * d;
        }
        if (query[ZLOW] < bounds[ZLOW]) {
            CoordType d = bounds[ZLOW] - query[ZHIGH];
            dist += d * d;
        } else if (query[ZLOW] > bounds[ZHIGH]) {
            CoordType d = query[ZLOW] - bounds[ZHIGH];
            dist += d * d;
        }
        return dist;
    }

    template <typename OutIter>
    void queryDumpAll(OutIter& out) const {
        if (isLeaf()) {
            for (const auto& item : data)
                *out++ = item.id;
        } else {
            for (const auto& child : children)
                if (child) child->queryDumpAll(out);
        }
    }

    template <typename OutIter>
    void query(const OBox<CoordType, Value>& qbox, OutIter out) const {
        queryImpl(qbox, out);
    }

    template <typename OutIter>
    void queryImpl(const OBox<CoordType, Value>& qbox, OutIter& out) const {
        const auto& queryBox = qbox.box;
        if (queryBox[XLOW] <= bounds[XLOW] && queryBox[XHIGH] >= bounds[XHIGH] &&
            queryBox[YLOW] <= bounds[YLOW] && queryBox[YHIGH] >= bounds[YHIGH] &&
            queryBox[ZLOW] <= bounds[ZLOW] && queryBox[ZHIGH] >= bounds[ZHIGH]) {
            queryDumpAll(out);
            return;
        }

        if (isLeaf()) {
            for (const auto& item : data) {
                if (item.intersects(queryBox)) 
                    *out++ = item.id;
            }
        } else {
            for (const auto& child : children) {
                if (child && !(child->bounds[XLOW] > queryBox[XHIGH] || child->bounds[XHIGH] < queryBox[XLOW] ||
                               child->bounds[YLOW] > queryBox[YHIGH] || child->bounds[YHIGH] < queryBox[YLOW] ||
                               child->bounds[ZLOW] > queryBox[ZHIGH] || child->bounds[ZHIGH] < queryBox[ZLOW])) {
                    child->queryImpl(qbox, out);
                }
            }
        }
    }

    bool remove(const OPoint<CoordType, Value>& p) {
        if (isLeaf()) {
            for (auto it = data.begin(); it != data.end(); ++it) {
                if (it->id == p.id && it->x == p.x && it->y == p.y && it->z == p.z) {
                    data.erase(it);
                    return true;
                }
            }
            return false;
        } else {
            for (auto c : children) {
                if (p.intersects(c->bounds)) {
                    if (c->remove(p)) {
                        consolidate();
                        return true;
                    }
                }
            }
            return false;
        }
    }

    void consolidate() {
        if (isLeaf()) return;
        
        size_t totalPoints = 0;
        bool allLeaves = true;
        
        for (auto c : children) {
            if (!c->isLeaf()) {
                allLeaves = false;
                break;
            }
            totalPoints += c->data.size();
        }
        
        if (allLeaves && totalPoints <= Capacity) {
            for (auto c : children) {
                data.insert(data.end(), c->data.begin(), c->data.end());
                delete c;
            }
            for (auto& c : children) c = nullptr;
        }
    }

    template <typename OutIter>
    void knnQuery(array<CoordType, 3> q, size_t k, OutIter out) {
        struct Neighbor {
            double dist;
            Value id;
        };
        auto neighbor_cmp = [](const Neighbor& a, const Neighbor& b) { return a.dist < b.dist; };
        vector<Neighbor> neighbors;
        neighbors.reserve(k);

        struct KnnNode {
            double dist;
            OctreeNode<CoordType, Value, Capacity>* sn;
            bool operator>(const KnnNode& o) const { return dist > o.dist; }
        };
        priority_queue<KnnNode, vector<KnnNode>, greater<KnnNode>> unseenNodes;

        auto ignoreNode = [&](double d) { return neighbors.size() == k && neighbors.front().dist <= d; };

        auto calcSqrDist = [](const array<CoordType, 3>& x, const array<CoordType, 3>& y) {
            CoordType dx = x[0] - y[0], dy = x[1] - y[1], dz = x[2] - y[2];
            return (double)(dx * dx + dy * dy + dz * dz);
        };

        array<CoordType, 6> query{q[0], q[1], q[2], q[0], q[1], q[2]};
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
                    double dd = calcSqrDist(p.toKNNPoint(), q);
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

    OctreeStats getStatistics() const {
        OctreeStats stats = {0, 0, 0, 0, sizeof(*this)};
        if (isLeaf()) {
            stats.numLeaves = 1;
            stats.numPoints = data.size();
            stats.sizeBytes += data.capacity() * sizeof(OPoint<CoordType, Value>);
        } else {
            stats.numInternalNodes = 1;
            for (auto c : children) {
                OctreeStats cs = c->getStatistics();
                stats.height = max(stats.height, cs.height);
                stats.numPoints += cs.numPoints;
                stats.numLeaves += cs.numLeaves;
                stats.numInternalNodes += cs.numInternalNodes;
                stats.sizeBytes += cs.sizeBytes;
            }
            stats.height += 1;
        }
        return stats;
    }
};

#endif
