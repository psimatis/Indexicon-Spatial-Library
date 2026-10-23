#include <array>
#include <fstream>
#include <iostream>
#include <vector>
#include "../../indexes/rtree_point.hpp"

using namespace std;

int main() {
    constexpr size_t Dims = 2;
    constexpr size_t Capacity = 128;

    // Load data
    ifstream fin("../data/input/sample_data.txt");
    if (!fin) { 
        cerr << "Cannot open sample_data.txt\n"; 
        return 1; 
    }

    // Each line: x y ...; read first 2 values, ignore the rest
    vector<array<double, Dims>> pts;
    double x, y, dummy;
    while (fin >> x >> y >> dummy >> dummy >> dummy >> dummy)
        pts.push_back({x, y});
    cout << "Loaded " << pts.size() << " points\n";

    using Tree = RTree<size_t, double, Dims, Capacity>;
    // Advanced configuration:
    // constexpr size_t MinElements = (Capacity * 3) / 10;
    // constexpr size_t ReinsertCount = (Capacity * 3) / 10;
    // constexpr size_t OverlapCostThreshold = 32;
    // constexpr size_t InternalNodeFanout = Capacity;
    // constexpr size_t MinInternalNodeFanout = (InternalNodeFanout * 3) / 10;
    // using Tree = RTree<size_t, double, Dims, Capacity, MinElements, ReinsertCount,
    //                    OverlapCostThreshold, InternalNodeFanout, MinInternalNodeFanout>;
    using Pt = Tree::PointType;
    using Box = Tree::BoxType;

    size_t split = pts.size() * 8 / 10;

    // Bulk load
    vector<pair<Pt, size_t>> bulk;
    bulk.reserve(split);
    for (size_t i = 0; i < split; ++i) {
        Pt p(pts[i][0], pts[i][1]);
        bulk.push_back({p, i});
    }

    Tree tree(bulk.begin(), bulk.end());
    cout << "Bulk loaded " << split << " points\n";

    // Insert
    for (size_t i = split; i < pts.size(); ++i) {
        Pt p(pts[i][0], pts[i][1]);
        tree.insert(p, i);
    }
    cout << "Inserted " << (pts.size() - split) << " more points\n";

    // Delete
    size_t del_count = 0;
    for (size_t i = 0; i < 100; ++i) {
        Pt p(pts[i][0], pts[i][1]);
        bool removed = tree.remove(p, i);
        if (removed)
            ++del_count;
    }
    cout << "Deleted " << del_count << " points\n";

    // Range query
    Box qbox;
    qbox.minCorner = Pt(-301.278286, 2883.648480);
    qbox.maxCorner = Pt(-221.070077, 2935.126180);
    vector<size_t> results;
    tree.query(qbox, back_inserter(results));
    cout << "Range query returned " << results.size() << " results\n";

    // kNN query
    Pt qpt(539.271901, 2915.970570);
    vector<size_t> knn;
    tree.knnQuery(qpt, 10, back_inserter(knn));
    cout << "10-NN query returned " << knn.size() << " results\n";

    // Statistics
    auto s = tree.getStatistics();
    cout << "\n=== R-tree (Point, 2D) statistics ===\n";
    cout << "  height          : " << s.height << "\n";
    cout << "  num points      : " << s.numPoints << "\n";
    cout << "  num leaves      : " << s.numLeaves << "\n";
    cout << "  num internals   : " << s.numInternalNodes << "\n";
    cout << "  memory (bytes)  : " << s.sizeBytes << "\n";
}
