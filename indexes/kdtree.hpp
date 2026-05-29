#ifndef KDTREE_HPP
#define KDTREE_HPP

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <queue>
#include <vector>

using namespace std;

#ifndef FIXED_VECTOR_DEFINED
#define FIXED_VECTOR_DEFINED
template <typename T, size_t Capacity>
struct FixedVector {
    array<T, Capacity> data;
    size_t count = 0;

    FixedVector() = default;

    size_t size() const { return count; }
    bool empty() const { return count == 0; }

    void reserve(size_t) {}
    void clear() { count = 0; }

    void push_back(const T& val) { data[count++] = val; }
    void push_back(T&& val) { data[count++] = move(val); }

    T& operator[](size_t i) { return data[i]; }
    const T& operator[](size_t i) const { return data[i]; }

    auto begin() { return data.begin(); }
    auto end() { return data.begin() + count; }
    auto begin() const { return data.begin(); }
    auto end() const { return data.begin() + count; }

    auto rbegin() { return reverse_iterator(end()); }
    auto rend() { return reverse_iterator(begin()); }

    template <typename Iter>
    void erase(Iter first, Iter last) {
        size_t offset = distance(data.begin(), first);
        size_t n = distance(first, last);
        for (size_t i = offset; i + n < count; ++i)
            data[i] = move(data[i + n]);
        count -= n;
    }

    template <typename Iter>
    void erase(Iter pos) { erase(pos, pos + 1); }

    template <typename Iter>
    void assign(Iter first, Iter last) {
        count = 0;
        for (Iter it = first; it != last; ++it)
            push_back(*it);
    }

    void resize(size_t new_size) { count = new_size; }
};
#endif

template <typename CoordType, size_t Dims>
struct KDPoint {
    array<CoordType, Dims> coords;
    KDPoint() { coords.fill(CoordType{}); }
    template <typename... Args>
    KDPoint(Args... args) : coords{{static_cast<CoordType>(args)...}} {}
    CoordType& operator[](size_t i) { return coords[i]; }
    const CoordType& operator[](size_t i) const { return coords[i]; }
};

template <typename CoordType, size_t Dims>
struct KDRegion {
    KDPoint<CoordType, Dims> low, high;

    KDRegion() {
        for (size_t i = 0; i < Dims; ++i) {
            low[i] = numeric_limits<CoordType>::lowest();
            high[i] = numeric_limits<CoordType>::max();
        }
    }

    bool contains(const KDPoint<CoordType, Dims>& pt) const {
        for (size_t i = 0; i < Dims; ++i) {
            if (pt[i] < low[i] || pt[i] > high[i]) 
                return false;
        }
        return true;
    }

    bool intersects(const KDRegion& other) const {
        for (size_t i = 0; i < Dims; ++i) {
            if (high[i] < other.low[i] || low[i] > other.high[i]) 
                return false;
        }
        return true;
    }

    bool containsRegion(const KDRegion& other) const {
        for (size_t i = 0; i < Dims; ++i) {
            if (other.low[i] < low[i] || other.high[i] > high[i]) 
                return false;
        }
        return true;
    }
};

template <typename CoordType, size_t Dims>
CoordType kdPointDistSq(const KDPoint<CoordType, Dims>& a, const KDPoint<CoordType, Dims>& b) {
    CoordType res = 0;
    for (size_t i = 0; i < Dims; ++i) { 
        CoordType d = a[i] - b[i]; 
        res += d * d; 
    }
    return res;
}

template <typename CoordType, size_t Dims>
CoordType kdMindistToRegion(const KDPoint<CoordType, Dims>& p, const KDRegion<CoordType, Dims>& reg) {
    CoordType dist = 0;
    for (size_t i = 0; i < Dims; ++i) {
        if (p[i] < reg.low[i]) {
            CoordType d = reg.low[i] - p[i];
            dist += d * d;
        }
        else if (p[i] > reg.high[i]) {
            CoordType d = p[i] - reg.high[i];
            dist += d * d;
        }
    }
    return dist;
}

struct KDTreeStats {
    size_t height;
    size_t numPoints;
    size_t numLeaves;
    size_t numInternalNodes;
    size_t sizeBytes;
};

enum class KDSplit { ADAPTIVE, ROUND_ROBIN, LONGEST_AXIS };

template <typename CoordType, typename Value, size_t Dims, size_t BucketSize = 128, KDSplit SplitMethod = KDSplit::ADAPTIVE>
class KDTree {
public:
    using PointPair = pair<KDPoint<CoordType, Dims>, Value>;
    using Region = KDRegion<CoordType, Dims>;
    using Point = KDPoint<CoordType, Dims>;

    struct Node {
        bool isLeaf;
        Region region;
        Node(bool leaf) : isLeaf(leaf) {}
    };

    struct LeafNode : Node {
        FixedVector<PointPair, BucketSize + 1> data;
        LeafNode() : Node(true) {}
    };

    struct InternalNode : Node {
        size_t splitDim;
        CoordType splitVal;
        unique_ptr<Node> left, right;
        InternalNode() : Node(false), splitDim(0), splitVal(0) {}

        ~InternalNode() {
            destroyChild(left);
            destroyChild(right);
        }

        static void destroyChild(unique_ptr<Node>& child) {
            if (!child) return;
            if (child->isLeaf)
                delete static_cast<LeafNode*>(child.release());
            else
                delete static_cast<InternalNode*>(child.release());
        }
    };

    unique_ptr<Node> root;
    size_t treeSize = 0;

    static LeafNode* asLeaf(Node* n) { return static_cast<LeafNode*>(n); }
    static const LeafNode* asLeaf(const Node* n) { return static_cast<const LeafNode*>(n); }
    static InternalNode* asInternal(Node* n) { return static_cast<InternalNode*>(n); }
    static const InternalNode* asInternal(const Node* n) { return static_cast<const InternalNode*>(n); }

    static unique_ptr<Node> makeLeaf() { return unique_ptr<Node>(new LeafNode()); }
    static unique_ptr<Node> makeInternal() { return unique_ptr<Node>(new InternalNode()); }

    KDTree() : root(makeLeaf()), treeSize(0) {}

    template <typename Iter>
    KDTree(Iter first, Iter last) : treeSize(distance(first, last)) {
        vector<PointPair> pts(first, last);
        root = buildRecursive(pts, 0, pts.size(), 0, Region());
    }

    void insert(const Point& pt, const Value& val) {
        insertImpl(root, pt, val, 0);
        treeSize++;
    }

    void insertImpl(unique_ptr<Node>& node, const Point& pt, const Value& val, size_t depth) {
        if (node->isLeaf) {
            auto* leaf = asLeaf(node.get());
            leaf->data.push_back({pt, val});
            if (leaf->data.size() > BucketSize)
                splitLeaf(node, depth);
            return;
        }
        auto* internal = asInternal(node.get());
        if (pt[internal->splitDim] < internal->splitVal)
            insertImpl(internal->left, pt, val, depth + 1);
        else
            insertImpl(internal->right, pt, val, depth + 1);
    }

    void remove(const Point& pt) {
        if (removeFrom(root, pt))
            treeSize--;
    }

    template <typename OutIter>
    void query(const Region& queryBox, OutIter out) const {
        queryRecursive(root, queryBox, out);
    }

    template <typename OutIter>
    void queryDumpAll(const unique_ptr<Node>& node, OutIter& out) const {
        if (!node)
            return;
        if (node->isLeaf) {
            const auto* leaf = asLeaf(node.get());
            for (const auto& p : leaf->data)
                *out++ = p.second;
        } else {
            const auto* internal = asInternal(node.get());
            queryDumpAll(internal->left, out);
            queryDumpAll(internal->right, out);
        }
    }

    template <typename OutIter>
    void queryRecursive(const unique_ptr<Node>& node, const Region& q, OutIter& out) const {
        if (!node)
            return;
        if (q.containsRegion(node->region)) {
            queryDumpAll(node, out);
            return;
        }
        if (node->isLeaf) {
            const auto* leaf = asLeaf(node.get());
            for (const auto& p : leaf->data)
                if (q.contains(p.first))
                    *out++ = p.second;
            return;
        }
        const auto* internal = asInternal(node.get());
        if (internal->left && q.low[internal->splitDim] <= internal->splitVal)
            queryRecursive(internal->left, q, out);
        if (internal->right && q.high[internal->splitDim] >= internal->splitVal)
            queryRecursive(internal->right, q, out);
    }

    template <typename OutIter>
    void knnQuery(Point pt, size_t k, OutIter out) const {
        struct Neighbor {
            CoordType dist;
            Value id;
        };
        auto neighbor_cmp = [](const Neighbor& a, const Neighbor& b) { return a.dist < b.dist; };
        vector<Neighbor> neighbors;
        neighbors.reserve(k);

        struct NodeEntry {
            CoordType dist;
            const Node* ptr;
            bool operator>(const NodeEntry& o) const { return dist > o.dist; }
        };
        priority_queue<NodeEntry, vector<NodeEntry>, greater<NodeEntry>> unseenNodes;

        auto ignoreNode = [&](CoordType d) { return neighbors.size() == k && neighbors.front().dist <= d; };

        NodeEntry start; start.dist = kdMindistToRegion(pt, root->region); start.ptr = root.get();
        unseenNodes.push(start);

        while (!unseenNodes.empty()) {
            CoordType dist = unseenNodes.top().dist;
            const Node* node = unseenNodes.top().ptr;
            unseenNodes.pop();

            if (ignoreNode(dist))
                break;

            if (node->isLeaf) {
                const auto* leaf = asLeaf(node);
                for (const auto& item : leaf->data) {
                    CoordType d = kdPointDistSq(pt, item.first);
                    if (neighbors.size() < k) {
                        Neighbor nb; nb.dist = d; nb.id = item.second;
                        neighbors.push_back(nb);
                        if (neighbors.size() == k)
                            make_heap(neighbors.begin(), neighbors.end(), neighbor_cmp);
                    } else if (d < neighbors.front().dist) {
                        pop_heap(neighbors.begin(), neighbors.end(), neighbor_cmp);
                        neighbors.back().dist = d;
                        neighbors.back().id = item.second;
                        push_heap(neighbors.begin(), neighbors.end(), neighbor_cmp);
                    }
                }
            } else {
                const auto* internal = asInternal(node);
                if (internal->left) {
                    CoordType dL = kdMindistToRegion(pt, internal->left->region);
                    if (!ignoreNode(dL)) {
                        NodeEntry ne; ne.dist = dL; ne.ptr = internal->left.get();
                        unseenNodes.push(ne);
                    }
                }
                if (internal->right) {
                    CoordType dR = kdMindistToRegion(pt, internal->right->region);
                    if (!ignoreNode(dR)) {
                        NodeEntry ne; ne.dist = dR; ne.ptr = internal->right.get();
                        unseenNodes.push(ne);
                    }
                }
            }
        }

        for (const auto& n : neighbors)
            *out++ = n.id;
    }

    size_t selectSplitDim(PointPair* pts, size_t n, size_t depth, const Region& reg) {
        if constexpr (SplitMethod == KDSplit::ROUND_ROBIN)
            return depth % Dims;
        if constexpr (SplitMethod == KDSplit::LONGEST_AXIS) {
            size_t bestDim = 0; CoordType maxSpread = -1;
            for (size_t i = 0; i < Dims; ++i) {
                CoordType spread = reg.high[i] - reg.low[i];
                if (spread > maxSpread) {
                    maxSpread = spread;
                    bestDim = i;
                }
            }
            return bestDim;
        }
        size_t bestDim = 0; CoordType maxSpread = -1;
        for (size_t i = 0; i < Dims; ++i) {
            CoordType minV = numeric_limits<CoordType>::max(), maxV = numeric_limits<CoordType>::lowest();
            for (size_t j = 0; j < n; ++j) {
                minV = min(minV, pts[j].first[i]);
                maxV = max(maxV, pts[j].first[i]);
            }
            if ((maxV - minV) > maxSpread) {
                maxSpread = maxV - minV;
                bestDim = i;
            }
        }
        return bestDim;
    }

    unique_ptr<Node> buildRecursive(vector<PointPair>& pts, size_t lo, size_t hi, size_t depth, Region reg) {
        size_t n = hi - lo;

        if (n <= BucketSize) {
            auto node = makeLeaf();
            node->region = reg;
            asLeaf(node.get())->data.assign(pts.begin() + lo, pts.begin() + hi);
            return node;
        }

        auto node = makeInternal();
        node->region = reg;
        auto* internal = asInternal(node.get());

        size_t dim = selectSplitDim(pts.data() + lo, n, depth, reg);
        internal->splitDim = dim;

        size_t mid = lo + n / 2;
        nth_element(pts.begin() + lo, pts.begin() + mid, pts.begin() + hi,
            [dim](const auto& a, const auto& b) { return a.first[dim] < b.first[dim]; });
        internal->splitVal = pts[mid].first[dim];

        Region leftReg = reg;
        leftReg.high[dim] = internal->splitVal;
        Region rightReg = reg;
        rightReg.low[dim] = internal->splitVal;

        internal->left = buildRecursive(pts, lo, mid, depth + 1, leftReg);
        internal->right = buildRecursive(pts, mid, hi, depth + 1, rightReg);
        return node;
    }

    void splitLeaf(unique_ptr<Node>& node, size_t depth) {
        auto* oldLeaf = asLeaf(node.get());
        size_t dim = selectSplitDim(&oldLeaf->data[0], oldLeaf->data.size(), depth, oldLeaf->region);
        size_t mid = oldLeaf->data.size() / 2;
        nth_element(oldLeaf->data.begin(), oldLeaf->data.begin() + mid, oldLeaf->data.end(),
            [dim](const auto& a, const auto& b) { return a.first[dim] < b.first[dim]; });
        CoordType splitVal = oldLeaf->data[mid].first[dim];

        auto left = makeLeaf();
        auto right = makeLeaf();
        left->region = oldLeaf->region;
        left->region.high[dim] = splitVal;
        right->region = oldLeaf->region;
        right->region.low[dim] = splitVal;

        auto* leftLeaf = asLeaf(left.get());
        auto* rightLeaf = asLeaf(right.get());
        for (size_t i = 0; i < mid; ++i)
            leftLeaf->data.push_back(move(oldLeaf->data[i]));
        for (size_t i = mid; i < oldLeaf->data.size(); ++i)
            rightLeaf->data.push_back(move(oldLeaf->data[i]));

        Region savedRegion = oldLeaf->region;

        auto newInternal = makeInternal();
        newInternal->region = savedRegion;
        auto* internalPtr = asInternal(newInternal.get());
        internalPtr->splitDim = dim;
        internalPtr->splitVal = splitVal;
        internalPtr->left = move(left);
        internalPtr->right = move(right);

        delete asLeaf(node.release());
        node = move(newInternal);
    }

    bool removeFrom(unique_ptr<Node>& node, const Point& pt) {
        if (!node)
            return false;
        if (node->isLeaf) {
            auto* leaf = asLeaf(node.get());
            auto it = find_if(leaf->data.begin(), leaf->data.end(), [&](const auto& p) {
                for (size_t i = 0; i < Dims; ++i) {
                    if (p.first[i] != pt[i])
                        return false;
                }
                return true;
            });
            if (it != leaf->data.end()) {
                leaf->data.erase(it);
                return true;
            }
            return false;
        }

        auto* internal = asInternal(node.get());
        bool removed = false;
        if (pt[internal->splitDim] <= internal->splitVal)
            removed = removeFrom(internal->left, pt);
        if (!removed && pt[internal->splitDim] >= internal->splitVal)
            removed = removeFrom(internal->right, pt);

        if (removed && internal->left && internal->right
            && internal->left->isLeaf && internal->right->isLeaf) {
            auto* leftLeaf = asLeaf(internal->left.get());
            auto* rightLeaf = asLeaf(internal->right.get());
            size_t total = leftLeaf->data.size() + rightLeaf->data.size();
            if (total <= BucketSize) {
                auto merged = makeLeaf();
                merged->region = node->region;
                auto* mergedLeaf = asLeaf(merged.get());
                for (size_t i = 0; i < leftLeaf->data.size(); ++i)
                    mergedLeaf->data.push_back(move(leftLeaf->data[i]));
                for (size_t i = 0; i < rightLeaf->data.size(); ++i)
                    mergedLeaf->data.push_back(move(rightLeaf->data[i]));

                delete asInternal(node.release());
                node = move(merged);
            }
        }

        return removed;
    }

    KDTreeStats getStatistics() const {
        KDTreeStats stats = {0, treeSize, 0, 0, 0};
        if (root)
            getStatsRecursive(root, 1, stats);
        stats.sizeBytes = stats.numLeaves * sizeof(LeafNode)
                        + stats.numInternalNodes * sizeof(InternalNode);
        return stats;
    }

    void getStatsRecursive(const unique_ptr<Node>& n, size_t h, KDTreeStats& s) const {
        s.height = max(s.height, h);
        if (n->isLeaf)
            s.numLeaves++;
        else {
            s.numInternalNodes++;
            const auto* internal = asInternal(n.get());
            if (internal->left)
                getStatsRecursive(internal->left, h + 1, s);
            if (internal->right)
                getStatsRecursive(internal->right, h + 1, s);
        }
    }
};

#endif
