#ifndef RTREE_POINT_HPP
#define RTREE_POINT_HPP

#include <algorithm>
#include <array>
#include <limits>
#include <queue>
#include <tuple>
#include <vector>

using namespace std;

template <typename Type, size_t Dims>
struct RPoint {
    array<Type, Dims> coords;

    RPoint() { coords.fill(Type{}); }

    template <typename... Args>
    RPoint(Args... args) : coords{{static_cast<Type>(args)...}} {}

    Type& operator[](size_t i) { return coords[i]; }
    const Type& operator[](size_t i) const { return coords[i]; }
};

template <typename Type, size_t Dims>
struct RBox {
    using PointType = RPoint<Type, Dims>;
    RPoint<Type, Dims> minCorner;
    RPoint<Type, Dims> maxCorner;

    RBox() {
        for (size_t i = 0; i < Dims; ++i) {
            minCorner[i] = numeric_limits<Type>::max();
            maxCorner[i] = numeric_limits<Type>::lowest();
        }
    }

    RBox(const RPoint<Type, Dims>& minPt, const RPoint<Type, Dims>& maxPt)
        : minCorner(minPt), maxCorner(maxPt) {}

    void expand(const RBox<Type, Dims>& other) {
        for (size_t i = 0; i < Dims; ++i) {
            minCorner[i] = min(minCorner[i], other.minCorner[i]);
            maxCorner[i] = max(maxCorner[i], other.maxCorner[i]);
        }
    }

    void expand(const RPoint<Type, Dims>& pt) {
        for (size_t i = 0; i < Dims; ++i) {
            minCorner[i] = min(minCorner[i], pt[i]);
            maxCorner[i] = max(maxCorner[i], pt[i]);
        }
    }

    Type area() const {
        Type result = Type{1};
        for (size_t i = 0; i < Dims; ++i) 
            result *= (maxCorner[i] - minCorner[i]);
        return result;
    }

    Type margin() const {
        Type result = Type{0};
        for (size_t i = 0; i < Dims; ++i) {
            Type face = Type{1};
            for (size_t j = 0; j < Dims; ++j) {
                if (i != j) 
                    face *= (maxCorner[j] - minCorner[j]);
            }
            result += face;
        }
        return result * Type{2};
    }

    RPoint<Type, Dims> center() const {
        RPoint<Type, Dims> c;
        for (size_t i = 0; i < Dims; ++i)
            c[i] = (minCorner[i] + maxCorner[i]) / Type{2};
        return c;
    }
};

template <typename Type, size_t Dims>
bool contains(const RBox<Type, Dims>& outer, const RBox<Type, Dims>& inner) {
    for (size_t i = 0; i < Dims; ++i) {
        if (outer.minCorner[i] > inner.minCorner[i] || outer.maxCorner[i] < inner.maxCorner[i])
            return false;
    }
    return true;
}

template <typename Type, size_t Dims>
bool contains(const RBox<Type, Dims>& box, const RPoint<Type, Dims>& point) {
    for (size_t i = 0; i < Dims; ++i) {
        if (point[i] < box.minCorner[i] || point[i] > box.maxCorner[i])
            return false;
    }
    return true;
}

template <typename Type, size_t Dims>
bool intersects(const RBox<Type, Dims>& a, const RBox<Type, Dims>& b) {
    for (size_t i = 0; i < Dims; ++i) {
        if (a.maxCorner[i] < b.minCorner[i] || a.minCorner[i] > b.maxCorner[i]) 
            return false;
    }
    return true;
}

template <typename Type, size_t Dims>
Type squaredDistance(const RPoint<Type, Dims>& a, const RPoint<Type, Dims>& b) {
    Type result = Type{0};
    for (size_t i = 0; i < Dims; ++i) {
        Type diff = a[i] - b[i];
        result += diff * diff;
    }
    return result;
}

template <typename Type, size_t Dims>
Type enlargement(const RBox<Type, Dims>& box, const RBox<Type, Dims>& add) {
    RBox<Type, Dims> expanded = box;
    expanded.expand(add);
    return expanded.area() - box.area();
}

template <typename Type, size_t Dims>
Type enlargement(const RBox<Type, Dims>& box, const RPoint<Type, Dims>& add) {
    RBox<Type, Dims> expanded = box;
    expanded.expand(add);
    return expanded.area() - box.area();
}

template <typename Type, size_t Dims>
Type mindistPointToBox(const RPoint<Type, Dims>& p, const RBox<Type, Dims>& box) {
    Type dist = Type{0};
    for (size_t i = 0; i < Dims; ++i) {
        if (p[i] < box.minCorner[i]) {
            Type diff = box.minCorner[i] - p[i];
            dist += diff * diff;
        } else if (p[i] > box.maxCorner[i]) {
            Type diff = p[i] - box.maxCorner[i];
            dist += diff * diff;
        }
    }
    return dist;
}

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

template <typename Value, typename BoxType>
struct Node {
    BoxType bounds;
    bool isLeaf;

    Node(bool leaf) : isLeaf(leaf) {}
    ~Node() = default;
};

template <typename Value, typename BoxType, typename PointType, size_t MaxCap>
struct LeafNode : Node<Value, BoxType> {
    FixedVector<pair<PointType, Value>, MaxCap> elements;

    LeafNode() : Node<Value, BoxType>(true) {}
};

template <typename Value, typename BoxType, typename PointType, size_t MaxLeafCap, size_t MaxInternalCap>
struct InternalNode : Node<Value, BoxType> {
    FixedVector<pair<BoxType, Node<Value, BoxType>*>, MaxInternalCap> children;

    InternalNode() : Node<Value, BoxType>(false) {}
    
    ~InternalNode() {
        for (auto& child : children) {
            if (child.second->isLeaf)
                delete static_cast<LeafNode<Value, BoxType, PointType, MaxLeafCap>*>(child.second);
            else
                delete static_cast<InternalNode<Value, BoxType, PointType, MaxLeafCap, MaxInternalCap>*>(child.second);
        }
    }
};

template <typename Value, typename Type, size_t Dims, size_t MaxElements = 128, size_t MinElements = (MaxElements * 3) / 10, size_t ReinsertCount = (MaxElements * 3) / 10, size_t OverlapCostThreshold = 32, size_t MaxFanout = MaxElements, size_t MinFanout = (MaxFanout * 3) / 10>
class RTree {
public:
    using PointType = RPoint<Type, Dims>;
    using BoxType = RBox<Type, Dims>;
    using NodeType = Node<Value, BoxType>;
    using LeafType = LeafNode<Value, BoxType, PointType, MaxElements + 1>;
    using InternalType = InternalNode<Value, BoxType, PointType, MaxElements + 1, MaxFanout + 1>;

    NodeType* root;
    size_t size;
    size_t height;
    static constexpr size_t maxElements = MaxElements;
    static constexpr size_t minElements = MinElements;
    static constexpr size_t reinsertCount = ReinsertCount;
    static constexpr size_t overlapCostThreshold = OverlapCostThreshold;
    static constexpr size_t maxFanout = MaxFanout;
    static constexpr size_t minFanout = MinFanout;

    RTree() : root(new LeafType()), size(0), height(1) {}

    template <typename Iter>
    RTree(Iter first, Iter last) : root(nullptr), size(0), height(0) {
        vector<pair<PointType, Iter>> entries;
        size = distance(first, last);
        entries.reserve(size);
        BoxType hintBox;
        for (Iter it = first; it != last; ++it) {
            hintBox.expand(it->first);
            entries.push_back({it->first, it});
        }

        size_t subtreeMax = 1;
        height = 1;
        if (size > maxElements) {
            subtreeMax = maxElements;
            height = 1;
            for (; subtreeMax < size; subtreeMax *= maxFanout)
                ++height;
            ++height;
        }

        root = packNode(entries.begin(), entries.end(), hintBox, size, subtreeMax);
    }

    static size_t biggestEdgeDim(const BoxType& box) {
        size_t dim = 0;
        Type maxLen = box.maxCorner[0] - box.minCorner[0];
        for (size_t i = 1; i < Dims; ++i) {
            Type len = box.maxCorner[i] - box.minCorner[i];
            if (len > maxLen) {
                maxLen = len;
                dim = i;
            }
        }
        return dim;
    }

    static size_t calculateMedianCount(size_t count, size_t subtreeMax, size_t subtreeMin) {
        size_t n = count / subtreeMax; 
        size_t r = count % subtreeMax; 

        if (r == 0)
            return (n / 2) * subtreeMax;

        if (subtreeMin <= r) 
            return ((n + 1) / 2) * subtreeMax;
        
        size_t countMinusMin = count - subtreeMin;
        r = countMinusMin % subtreeMax;
        if (r == 0)
            return ((n + 1) / 2) * subtreeMax;
        n = countMinusMin / subtreeMax;
        if (n == 0)
            return r;
        
        return ((n + 2) / 2) * subtreeMax;
    }    

    template <typename EIt>
    NodeType* packNode(EIt first, EIt last, const BoxType& hintBox, size_t count, size_t subtreeMax) {
        // Create leaf node
        if (subtreeMax <= 1) {
            auto* leaf = new LeafType();
            leaf->elements.reserve(count);
            for (EIt it = first; it != last; ++it) {
                leaf->elements.push_back(*(it->second));
                leaf->bounds.expand(it->second->first);
            }
            return leaf;
        }

        // Create internal node
        size_t nextSubtreeMax = (subtreeMax <= maxElements) ? 1 : subtreeMax / maxFanout;
        size_t subtreeMin = (subtreeMax <= maxElements ? minElements : minFanout) * nextSubtreeMax;
        auto* internal = new InternalType();
        pack(first, last, hintBox, count, subtreeMax, subtreeMin, nextSubtreeMax, internal);
        return internal;
    }

    template <typename EIt>
    void pack(EIt first, EIt last, const BoxType& hintBox, size_t count, size_t subtreeMax, size_t subtreeMin, size_t nextSubtreeMax, InternalType* internal) {
        if (count <= subtreeMax) {
            NodeType* child = packNode(first, last, hintBox, count, nextSubtreeMax);
            internal->children.push_back({child->bounds, child});
            internal->bounds.expand(child->bounds);
            return;
        }

        size_t medianCount = calculateMedianCount(count, subtreeMax, subtreeMin);
        EIt median = first + medianCount;

        size_t dim = biggestEdgeDim(hintBox);
        nth_element(first, median, last, [dim](const auto& a, const auto& b) {return a.first[dim] < b.first[dim];});

        BoxType leftBox = hintBox;
        BoxType rightBox = hintBox;
        Type mid = hintBox.center()[dim];
        leftBox.maxCorner[dim] = mid;
        rightBox.minCorner[dim] = mid;

        pack(first, median, leftBox, medianCount, subtreeMax, subtreeMin, nextSubtreeMax, internal);
        pack(median, last, rightBox, count - medianCount, subtreeMax, subtreeMin, nextSubtreeMax, internal);
    }

    ~RTree() {
        delete root;
    }

    void insert(const PointType& pt, const Value& value) {
        vector<pair<PointType, Value>> leafReinserts;
        vector<tuple<BoxType, NodeType*, size_t>> internalReinserts;

        NodeType* splitNode = insertImpl(root, pt, value, height, leafReinserts, internalReinserts);
        ++size;
        handleRootSplit(splitNode);

        processInternalReinserts(internalReinserts);

        for (auto it = leafReinserts.rbegin(); it != leafReinserts.rend(); ++it) {
            internalReinserts.clear();
            NodeType* s = insertNoLeafReinsert(root, it->first, it->second, height, internalReinserts);
            handleRootSplit(s);
            processInternalReinserts(internalReinserts);
        }
    }
        
    void handleRootSplit(NodeType* newNode) {
        if (newNode) {
            auto* newRoot = new InternalType();
            newRoot->children.push_back({root->bounds, root});
            newRoot->children.push_back({newNode->bounds, newNode});
            updateBounds(newRoot);
            root = newRoot;
            ++height;
        }
    }

    NodeType* insertImpl(NodeType* node, const PointType& pt, const Value& value, size_t level, vector<pair<PointType, Value>>& leafReinserts, vector<tuple<BoxType, NodeType*, size_t>>& internalReinserts) {
        if (node->isLeaf) {
            auto* leaf = static_cast<LeafType*>(node);
            leaf->elements.push_back({pt, value});
            leaf->bounds.expand(pt);

            if (leaf->elements.size() > maxElements) {
                if (reinsertCount > 0 && level < height) {
                    reinsertLeaf(leaf, leafReinserts);
                    return nullptr;
                }
                return splitLeaf(leaf);
            }
            return nullptr;
        } else {
            auto* internal = static_cast<InternalType*>(node);
            size_t bestIdx = chooseSubtree(internal, pt, level - 1);
            NodeType* childSplit = insertImpl(internal->children[bestIdx].second, pt, value, level - 1, leafReinserts, internalReinserts);

            internal->children[bestIdx].first = internal->children[bestIdx].second->bounds;

            if (childSplit) {
                internal->children.push_back({childSplit->bounds, childSplit});
                updateBounds(internal);
            } else 
                internal->bounds.expand(pt);

            if (internal->children.size() > maxFanout) {
                if (reinsertCount > 0 && level < height) {
                    reinsertInternal(internal, internalReinserts, level - 1);
                    return nullptr;
                }
                return splitInternal(internal);
            }
            return nullptr;
        }
    }

    void reinsertLeaf(LeafType* leaf, vector<pair<PointType, Value>>& reinsertList) {
        auto& elements = leaf->elements;
        PointType center = leaf->bounds.center();

        size_t removeCount = min(reinsertCount, elements.size() - minElements);

        typedef pair<Type, size_t> distPair;
        FixedVector<distPair, MaxElements + 1> sorted;
        sorted.resize(elements.size());
        for (size_t i = 0; i < elements.size(); ++i)
            sorted[i] = {squaredDistance(elements[i].first, center), i};

        partial_sort(sorted.begin(), sorted.begin() + removeCount, sorted.end(), [](const distPair& a, const distPair& b) { return a.first > b.first; });

        for (size_t i = 0; i < removeCount; ++i)
            reinsertList.push_back(elements[sorted[i].second]);

        FixedVector<pair<PointType, Value>, MaxElements + 1> remaining;
        remaining.resize(elements.size() - removeCount);
        for (size_t i = removeCount; i < sorted.size(); ++i)
            remaining[i - removeCount] = elements[sorted[i].second];
        elements.clear();
        for (const auto& r : remaining)
            elements.push_back(r);

        updateBounds(leaf);
    }

    void reinsertInternal(InternalType* internal, vector<tuple<BoxType, NodeType*, size_t>>& reinsertList, size_t child_level) {
        auto& children = internal->children;
        PointType center = internal->bounds.center();

        size_t removeCount = min(reinsertCount, children.size() - minFanout);

        typedef pair<Type, size_t> distPair;
        FixedVector<distPair, MaxFanout + 1> sorted;
        sorted.resize(children.size());
        for (size_t i = 0; i < children.size(); ++i)
            sorted[i] = {squaredDistance(children[i].first.center(), center), i};

        partial_sort(sorted.begin(), sorted.begin() + removeCount, sorted.end(),[](const distPair& a, const distPair& b) { return a.first > b.first; });

        for (size_t i = 0; i < removeCount; ++i) {
            size_t idx = sorted[i].second;
            reinsertList.push_back({children[idx].first, children[idx].second, child_level});
        }

        FixedVector<pair<BoxType, NodeType*>, MaxFanout + 1> remaining;
        remaining.resize(children.size() - removeCount);
        for (size_t i = removeCount; i < sorted.size(); ++i)
            remaining[i - removeCount] = children[sorted[i].second];
        children.clear();
        for (const auto& r : remaining)
            children.push_back(r);

        updateBounds(internal);
    }

    void processInternalReinserts(vector<tuple<BoxType, NodeType*, size_t>>& entries) {
        for (auto it = entries.rbegin(); it != entries.rend(); ++it) {
            auto& [b, n, lvl] = *it;
            NodeType* s = insertSubtreeNoReinsert(root, b, n, height, lvl);
            handleRootSplit(s);
        }
    }

    NodeType* insertNoLeafReinsert(NodeType* node, const PointType& pt, const Value& value, size_t level, vector<tuple<BoxType, NodeType*, size_t>>& internalReinserts) {
        if (node->isLeaf) {
            auto* leaf = static_cast<LeafType*>(node);
            leaf->elements.push_back({pt, value});
            leaf->bounds.expand(pt);

            if (leaf->elements.size() > maxElements)
                return splitLeaf(leaf);
            return nullptr;
        } else {
            auto* internal = static_cast<InternalType*>(node);
            size_t bestIdx = chooseSubtree(internal, pt, level - 1);
            NodeType* childSplit = insertNoLeafReinsert(internal->children[bestIdx].second, pt, value, level - 1, internalReinserts);

            internal->children[bestIdx].first = internal->children[bestIdx].second->bounds;

            if (childSplit) {
                internal->children.push_back({childSplit->bounds, childSplit});
                updateBounds(internal);
            } else
                internal->bounds.expand(pt);

            if (internal->children.size() > maxFanout) {
                if (reinsertCount > 0 && level < height) {
                    reinsertInternal(internal, internalReinserts, level - 1);
                    return nullptr;
                }
                return splitInternal(internal);
            }
            return nullptr;
        }
    }

    NodeType* insertSubtreeNoReinsert(NodeType* node, const BoxType& bounds, NodeType* subtree, size_t currentLevel, size_t subtree_level) {
        if (currentLevel == subtree_level + 1) {
            auto* internal = static_cast<InternalType*>(node);
            internal->children.push_back({bounds, subtree});
            internal->bounds.expand(bounds);
            if (internal->children.size() > maxFanout)
                return splitInternal(internal);
            return nullptr;
        } else {
            auto* internal = static_cast<InternalType*>(node);
            size_t bestIdx = chooseSubtree(internal, bounds, currentLevel - 1);
            NodeType* childSplit = insertSubtreeNoReinsert(
                internal->children[bestIdx].second, bounds, subtree, currentLevel - 1, subtree_level);
            internal->children[bestIdx].first = internal->children[bestIdx].second->bounds;
            if (childSplit) {
                internal->children.push_back({childSplit->bounds, childSplit});
                updateBounds(internal);
            } else 
                internal->bounds.expand(bounds);
            if (internal->children.size() > maxFanout)
                return splitInternal(internal);
            return nullptr;
        }
    }

    template <typename Container, typename BoxContainer>
    static void computePrefixSuffixPoints(Container& items, BoxContainer& prefix, BoxContainer& suffix, size_t axis) {
        sort(items.begin(), items.end(), [axis](const auto& a, const auto& b) { return a.first[axis] < b.first[axis]; });
        size_t n = items.size();
        prefix[0] = BoxType();
        for (size_t i = 0; i < n; ++i) { 
            prefix[i+1] = prefix[i]; 
            prefix[i+1].expand(items[i].first); 
        }
        suffix[n] = BoxType();
        for (size_t i = n; i > 0; --i) { 
            suffix[i-1] = suffix[i]; 
            suffix[i-1].expand(items[i-1].first); 
        }
    }

    template <typename Container, typename BoxContainer>
    static void computePrefixSuffixBoxes(Container& items, BoxContainer& prefix, BoxContainer& suffix, size_t axis, bool useMax) {
        if (useMax)
            sort(items.begin(), items.end(), [axis](const auto& a, const auto& b) { return a.first.maxCorner[axis] < b.first.maxCorner[axis]; });
        else
            sort(items.begin(), items.end(), [axis](const auto& a, const auto& b) { return a.first.minCorner[axis] < b.first.minCorner[axis]; });
        size_t n = items.size();
        prefix[0] = BoxType();
        for (size_t i = 0; i < n; ++i) { 
            prefix[i+1] = prefix[i]; 
            prefix[i+1].expand(items[i].first); 
        }
        suffix[n] = BoxType();
        for (size_t i = n; i > 0; --i) { 
            suffix[i-1] = suffix[i]; 
            suffix[i-1].expand(items[i-1].first); 
        }
    }

    size_t chooseSubtree(InternalType* node, const PointType& pt, size_t targetLevel) {
        const auto& children = node->children;
        
        if (targetLevel == 1) 
            return chooseSubtreeOverlap(node, pt);
        
        size_t bestIdx = 0;
        Type bestEnlargement = numeric_limits<Type>::max();
        Type bestContent = numeric_limits<Type>::max();

        for (size_t i = 0; i < children.size(); ++i) {
            Type enl = enlargement(children[i].first, pt);
            Type content = children[i].first.area();

            if (enl < bestEnlargement || (enl == bestEnlargement && content < bestContent)) {
                bestEnlargement = enl;
                bestContent = content;
                bestIdx = i;
            }
        }
        return bestIdx;
    }

    size_t chooseSubtree(InternalType* node, const BoxType& bounds, size_t targetLevel) {
        const auto& children = node->children;
        
        if (targetLevel == 1) 
            return chooseSubtreeOverlap(node, bounds);
        
        size_t bestIdx = 0;
        Type bestEnlargement = numeric_limits<Type>::max();
        Type bestContent = numeric_limits<Type>::max();

        for (size_t i = 0; i < children.size(); ++i) {
            Type enl = enlargement(children[i].first, bounds);
            Type content = children[i].first.area();

            if (enl < bestEnlargement || (enl == bestEnlargement && content < bestContent)) {
                bestEnlargement = enl;
                bestContent = content;
                bestIdx = i;
            }
        }
        return bestIdx;
    }

    size_t chooseSubtreeOverlap(InternalType* node, const PointType& pt) {
        const auto& children = node->children;
        const size_t n = children.size();
        
        struct ChildInfo {
            size_t idx;
            Type enlargement;
            Type content;
        };
        FixedVector<ChildInfo, MaxFanout + 1> infos;
        infos.resize(n);
        
        size_t bestIdx = 0;
        Type bestEnlargement = numeric_limits<Type>::max();
        Type bestContent = numeric_limits<Type>::max();
        for (size_t i = 0; i < n; ++i) {
            BoxType expanded = children[i].first;
            expanded.expand(pt);
            Type enl = expanded.area() - children[i].first.area();
            Type content = expanded.area();
            infos[i] = {i, enl, content};
            
            if (enl < bestEnlargement || (enl == bestEnlargement && content < bestContent)) {
                bestEnlargement = enl;
                bestContent = content;
                bestIdx = i;
            }
        }
        if (bestEnlargement <= numeric_limits<Type>::epsilon())
            return bestIdx;
        
        size_t check_count = min(n, overlapCostThreshold);
        nth_element(infos.begin(), infos.begin() + check_count, infos.end(),
            [](const ChildInfo& a, const ChildInfo& b) {
                return a.enlargement < b.enlargement || (a.enlargement == b.enlargement && a.content < b.content);
            });
        
        bestIdx = infos[0].idx;
        Type bestOverlapIncrease = numeric_limits<Type>::max();
        bestEnlargement = numeric_limits<Type>::max();
        bestContent = numeric_limits<Type>::max();
        
        for (size_t k = 0; k < check_count; ++k) {
            size_t i = infos[k].idx;
            BoxType expanded = children[i].first;
            expanded.expand(pt);
            
            Type overlapDiff = Type{0};
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    Type ov_after = intersectionContent(expanded, children[j].first);
                    if (ov_after)
                        overlapDiff += ov_after - intersectionContent(children[i].first, children[j].first);
                }
            }
            
            Type enl = infos[k].enlargement;
            Type content = infos[k].content;
            if (overlapDiff < bestOverlapIncrease ||
                (overlapDiff == bestOverlapIncrease && enl < bestEnlargement) ||
                (overlapDiff == bestOverlapIncrease && enl == bestEnlargement && content < bestContent)) {
                bestOverlapIncrease = overlapDiff;
                bestEnlargement = enl;
                bestContent = content;
                bestIdx = i;
            }
        }
        return bestIdx;
    }

    size_t chooseSubtreeOverlap(InternalType* node, const BoxType& bounds) {
        const auto& children = node->children;
        const size_t n = children.size();
        
        struct ChildInfo {
            size_t idx;
            Type enlargement;
            Type content;
        };
        FixedVector<ChildInfo, MaxFanout + 1> infos;
        infos.resize(n);
        
        size_t bestIdx = 0;
        Type bestEnlargement = numeric_limits<Type>::max();
        Type bestContent = numeric_limits<Type>::max();
        for (size_t i = 0; i < n; ++i) {
            BoxType expanded = children[i].first;
            expanded.expand(bounds);
            Type enl = expanded.area() - children[i].first.area();
            Type content = expanded.area();
            infos[i] = {i, enl, content};
            
            if (enl < bestEnlargement || (enl == bestEnlargement && content < bestContent)) {
                bestEnlargement = enl;
                bestContent = content;
                bestIdx = i;
            }
        }
        if (bestEnlargement <= numeric_limits<Type>::epsilon())
            return bestIdx;
        
        size_t check_count = min(n, overlapCostThreshold);
        nth_element(infos.begin(), infos.begin() + check_count, infos.end(),
            [](const ChildInfo& a, const ChildInfo& b) {
                return a.enlargement < b.enlargement || (a.enlargement == b.enlargement && a.content < b.content);
            });
        
        bestIdx = infos[0].idx;
        Type bestOverlapIncrease = numeric_limits<Type>::max();
        bestEnlargement = numeric_limits<Type>::max();
        bestContent = numeric_limits<Type>::max();
        
        for (size_t k = 0; k < check_count; ++k) {
            size_t i = infos[k].idx;
            BoxType expanded = children[i].first;
            expanded.expand(bounds);
            
            Type overlapDiff = Type{0};
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    Type ov_after = intersectionContent(expanded, children[j].first);
                    if (ov_after)
                        overlapDiff += ov_after - intersectionContent(children[i].first, children[j].first);
                }
            }
            
            Type enl = infos[k].enlargement;
            Type content = infos[k].content;
            if (overlapDiff < bestOverlapIncrease ||
                (overlapDiff == bestOverlapIncrease && enl < bestEnlargement) ||
                (overlapDiff == bestOverlapIncrease && enl == bestEnlargement && content < bestContent)) {
                bestOverlapIncrease = overlapDiff;
                bestEnlargement = enl;
                bestContent = content;
                bestIdx = i;
            }
        }
        return bestIdx;
    }

    static Type intersectionContent(const BoxType& a, const BoxType& b) {
        Type result = Type{1};
        for (size_t d = 0; d < Dims; ++d) {
            Type lo = max(a.minCorner[d], b.minCorner[d]);
            Type hi = min(a.maxCorner[d], b.maxCorner[d]);
            if (lo >= hi) 
                return Type{0};
            result *= (hi - lo);
        }
        return result;
    }

    NodeType* splitLeaf(LeafType* leaf) {
        auto& elements = leaf->elements;
        
        FixedVector<pair<PointType, Value>, MaxElements + 1> elements_copy;
        elements_copy.assign(elements.begin(), elements.end());
        const size_t n = elements_copy.size();

        FixedVector<BoxType, MaxElements + 2> prefix; prefix.resize(n + 1);
        FixedVector<BoxType, MaxElements + 2> suffix; suffix.resize(n + 1);
        size_t bestAxis = 0;
        Type bestMarginSum = numeric_limits<Type>::max();
        for (size_t axis = 0; axis < Dims; ++axis) {
            Type marginSum = 0;
            computePrefixSuffixPoints(elements_copy, prefix, suffix, axis);
            for (size_t split = minElements; split <= n - minElements; ++split)
                marginSum += prefix[split].margin() + suffix[split].margin();
            
            if (marginSum < bestMarginSum) {
                bestMarginSum = marginSum;
                bestAxis = axis;
            }
        }

        size_t bestSplit = minElements;
        Type bestOverlap = numeric_limits<Type>::max();
        Type bestContent = numeric_limits<Type>::max();

        computePrefixSuffixPoints(elements_copy, prefix, suffix, bestAxis);

        for (size_t split = minElements; split <= n - minElements; ++split) {
            Type overlap = intersectionContent(prefix[split], suffix[split]);
            Type content = prefix[split].area() + suffix[split].area();

            if (overlap < bestOverlap || (overlap == bestOverlap && content < bestContent)) {
                bestOverlap = overlap;
                bestContent = content;
                bestSplit = split;
            }
        }

        nth_element(elements.begin(), elements.begin() + bestSplit, elements.end(), [bestAxis](const auto& a, const auto& b) { return a.first[bestAxis] < b.first[bestAxis]; });

        auto* new_leaf = new LeafType();
        new_leaf->elements.assign(elements.begin() + bestSplit, elements.end());
        elements.resize(bestSplit);

        updateBounds(leaf);
        updateBounds(new_leaf);

        return new_leaf;
    }

    NodeType* splitInternal(InternalType* internal) {
        auto& children = internal->children;
        
        FixedVector<pair<BoxType, NodeType*>, MaxFanout + 1> children_copy;
        children_copy.assign(children.begin(), children.end());
        const size_t n = children_copy.size();

        FixedVector<BoxType, MaxFanout + 2> prefix; prefix.resize(n + 1);
        FixedVector<BoxType, MaxFanout + 2> suffix; suffix.resize(n + 1);

        size_t bestAxis = 0;
        Type bestMarginSum = numeric_limits<Type>::max();
        for (size_t axis = 0; axis < Dims; ++axis) {
            Type marginSum = 0;
            computePrefixSuffixBoxes(children_copy, prefix, suffix, axis, false);
            for (size_t split = minFanout; split <= n - minFanout; ++split)
                marginSum += prefix[split].margin() + suffix[split].margin();
            computePrefixSuffixBoxes(children_copy, prefix, suffix, axis, true);
            for (size_t split = minFanout; split <= n - minFanout; ++split)
                marginSum += prefix[split].margin() + suffix[split].margin();
            if (marginSum < bestMarginSum) {
                bestMarginSum = marginSum;
                bestAxis = axis;
            }
        }

        bool bestUseMax = false;
        size_t bestSplit = minFanout;
        Type bestOverlap = numeric_limits<Type>::max();
        Type bestContent = numeric_limits<Type>::max();

        computePrefixSuffixBoxes(children_copy, prefix, suffix, bestAxis, false);

        for (size_t split = minFanout; split <= n - minFanout; ++split) {
            Type overlap = intersectionContent(prefix[split], suffix[split]);
            Type content = prefix[split].area() + suffix[split].area();

            if (overlap < bestOverlap || (overlap == bestOverlap && content < bestContent)) {
                bestOverlap = overlap;
                bestContent = content;
                bestSplit = split;
                bestUseMax = false;
            }
        }

        computePrefixSuffixBoxes(children_copy, prefix, suffix, bestAxis, true);

        for (size_t split = minFanout; split <= n - minFanout; ++split) {
            Type overlap = intersectionContent(prefix[split], suffix[split]);
            Type content = prefix[split].area() + suffix[split].area();

            if (overlap < bestOverlap || (overlap == bestOverlap && content < bestContent)) {
                bestOverlap = overlap;
                bestContent = content;
                bestSplit = split;
                bestUseMax = true;
            }
        }

        if (bestUseMax)
            nth_element(children.begin(), children.begin() + bestSplit, children.end(), [bestAxis](const auto& a, const auto& b) { return a.first.maxCorner[bestAxis] < b.first.maxCorner[bestAxis]; });
        else
            nth_element(children.begin(), children.begin() + bestSplit, children.end(), [bestAxis](const auto& a, const auto& b) { return a.first.minCorner[bestAxis] < b.first.minCorner[bestAxis]; });

        auto* newInternal = new InternalType();
        for (size_t i = bestSplit; i < children.size(); ++i) 
            newInternal->children.push_back(move(children[i]));
        children.resize(bestSplit);

        updateBounds(internal);
        updateBounds(newInternal);

        return newInternal;
    }

    void updateBounds(LeafType* leaf) {
        leaf->bounds = BoxType();
        for (const auto& elem : leaf->elements) 
            leaf->bounds.expand(elem.first);
    }

    void updateBounds(InternalType* internal) {
        internal->bounds = BoxType();
        for (const auto& child : internal->children) 
            internal->bounds.expand(child.first);
    }

    bool remove(const PointType& pt, const Value& value) {
        vector<pair<PointType, Value>> orphans;
        bool found = removeImpl(root, pt, value, orphans);

        if (found) {
            --size;
            for (auto& orphan : orphans)
                insert(orphan.first, orphan.second);
            condenseRoot();
        }
        return found;
    }

    bool removeImpl(NodeType* node, const PointType& pt, const Value& value, vector<pair<PointType, Value>>& orphans) {
        if (node->isLeaf) {
            auto* leaf = static_cast<LeafType*>(node);
            for (auto it = leaf->elements.begin(); it != leaf->elements.end(); ++it) {
                if (it->second == value && it->first.coords == pt.coords) {
                    leaf->elements.erase(it);
                    updateBounds(leaf);
                    return true;
                }
            }
            return false;
        } else {
            auto* internal = static_cast<InternalType*>(node);
            for (size_t i = 0; i < internal->children.size(); ++i) {
                if (contains(internal->children[i].first, pt)) {
                    if (removeImpl(internal->children[i].second, pt, value, orphans)) {
                        NodeType* child = internal->children[i].second;
                        size_t child_size = child->isLeaf
                            ? static_cast<LeafType*>(child)->elements.size()
                            : static_cast<InternalType*>(child)->children.size();

                        if (child_size < (child->isLeaf ? minElements : minFanout) && internal->children.size() > 1) {
                            collectOrphans(child, orphans);
                            delete child;
                            internal->children.erase(internal->children.begin() + i);
                        } else 
                            internal->children[i].first = child->bounds;
                        updateBounds(internal);
                        return true;
                    }
                }
            }
            return false;
        }
    }

    void collectOrphans(NodeType* node, vector<pair<PointType, Value>>& orphans) {
        if (node->isLeaf) {
            auto* leaf = static_cast<LeafType*>(node);
            for (auto& elem : leaf->elements) 
                orphans.push_back(move(elem));
        } else {
            auto* internal = static_cast<InternalType*>(node);
            for (auto& child : internal->children) 
                collectOrphans(child.second, orphans);
            internal->children.clear();
        }
    }

    void condenseRoot() {
        while (!root->isLeaf) {
            auto* internal = static_cast<InternalType*>(root);
            if (internal->children.size() == 1) {
                NodeType* new_root = internal->children[0].second;
                internal->children.clear();
                delete root;
                root = new_root;
                --height;
            } else 
                break;
        }
    }

    struct RTreeStats {
        size_t height;
        size_t numPoints;
        size_t numLeaves;
        size_t numInternalNodes;
        size_t sizeBytes;
    };

    RTreeStats getStatistics() const {
        size_t leafCount = 0, internalCount = 0, internalChildren = 0;
        
        vector<const NodeType*> stack;
        stack.push_back(root);
        
        while (!stack.empty()) {
            const NodeType* node = stack.back();
            stack.pop_back();
            
            if (node->isLeaf)
                ++leafCount;
            else {
                const auto* internal = static_cast<const InternalType*>(node);
                ++internalCount;
                internalChildren += internal->children.size();
                for (const auto& child : internal->children)
                    stack.push_back(child.second);
            }
        }
        
        size_t bytes = leafCount * sizeof(LeafType);
        bytes += internalCount * sizeof(InternalType);
        
        return {height, size, leafCount, internalCount, bytes};
    }

    template <typename OutIter>
    void query(const BoxType& query_box, OutIter out) const {
        queryImpl(root, query_box, out);
    }

    template <typename OutIter>
    void queryPoint(const PointType& pt, OutIter out) const {
        queryPointImpl(root, pt, out);
    }

    template <typename OutIter>
    void queryDumpAll(const NodeType* node, OutIter& out) const {
        if (node->isLeaf) {
            const auto& elems = static_cast<const LeafType*>(node)->elements;
            for (const auto& elem : elems)
                *out++ = elem.second;
        } else {
            const auto& children = static_cast<const InternalType*>(node)->children;
            for (const auto& child : children)
                queryDumpAll(child.second, out);
        }
    }

    template <typename OutIter>
    void queryPointImpl(const NodeType* node, const PointType& pt, OutIter& out) const {
        if (node->isLeaf) {
            const auto& elems = static_cast<const LeafType*>(node)->elements;
            for (const auto& elem : elems) {
                if (elem.first.coords == pt.coords)
                    *out++ = elem.second;
            }
        } else {
            const auto& children = static_cast<const InternalType*>(node)->children;
            for (const auto& child : children) {
                if (contains(child.first, pt))
                    queryPointImpl(child.second, pt, out);
            }
        }
    }

    template <typename OutIter>
    void queryImpl(const NodeType* node, const BoxType& query_box, OutIter& out) const {
        if (node->isLeaf) {
            const auto& elems = static_cast<const LeafType*>(node)->elements;
            if (contains(query_box, node->bounds)) {
                for (const auto& elem : elems)
                    *out++ = elem.second;
            } else {
                for (const auto& elem : elems) {
                    if (contains(query_box, elem.first)) 
                        *out++ = elem.second;
                }
            }
        } else {
            const auto& children = static_cast<const InternalType*>(node)->children;
            for (const auto& child : children) {
                if (contains(query_box, child.first))
                    queryDumpAll(child.second, out);
                else if (intersects(child.first, query_box))
                    queryImpl(child.second, query_box, out);
            }
        }
    }

    template <typename OutIter>
    void knnQuery(const PointType& query, size_t k, OutIter out) const {
        struct neighbor {
            Type dist;
            const Value* ptr;
        };
        auto neighbor_cmp = [](const neighbor& a, const neighbor& b) {return a.dist < b.dist;};
        vector<neighbor> neighbors;
        neighbors.reserve(k);

        struct node {
            Type dist;
            size_t level;
            const NodeType* ptr;
        };
        auto node_cmp = [](const node& a, const node& b) {return a.dist > b.dist || (a.dist == b.dist && a.level > b.level);};
        priority_queue<node, vector<node>, decltype(node_cmp)> queue(node_cmp);
        auto ignoreNode = [&](Type d) {return neighbors.size() == k && neighbors.front().dist <= d;};

        const NodeType* ptr = root;
        size_t rl = height - 1;
        while (true) {
            if (rl > 0) {
                const auto* internal = static_cast<const InternalType*>(ptr);
                for (const auto& child : internal->children) {
                    Type d = mindistPointToBox(query, child.first);
                    if (!ignoreNode(d))
                        queue.push({d, rl - 1, child.second});
                }
            } else {
                const auto* leaf = static_cast<const LeafType*>(ptr);
                for (const auto& elem : leaf->elements) {
                    Type d = squaredDistance(query, elem.first);
                    if (neighbors.size() < k) {
                        neighbors.push_back({d, &elem.second});
                        if (neighbors.size() == k)
                            make_heap(neighbors.begin(), neighbors.end(), neighbor_cmp);
                    } else if (d < neighbors.front().dist) {
                        pop_heap(neighbors.begin(), neighbors.end(), neighbor_cmp);
                        neighbors.back() = {d, &elem.second};
                        push_heap(neighbors.begin(), neighbors.end(), neighbor_cmp);
                    }
                }
            }

            if (queue.empty() || ignoreNode(queue.top().dist))
                break;

            const auto& top = queue.top();
            ptr = top.ptr;
            rl = top.level;
            queue.pop();
        }

        for (const auto& n : neighbors)
            *out++ = *(n.ptr);
    }
};

#endif // RTREE_POINT_HPP
