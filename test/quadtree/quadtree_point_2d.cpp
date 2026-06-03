#include <array>
#include <fstream>
#include <iostream>
#include <vector>
#include "../../indexes/quadtree.hpp"

using namespace std;

int main() {
    constexpr size_t Dims = 2;
    constexpr size_t Capacity = 128;
    constexpr SplitStrategy Strategy = SplitStrategy::PR;
    // Other strategies: SplitStrategy::PSEUDO_MEDIAN, SplitStrategy::POINT_LONGEST_AXIS

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

    using Tree = QuadTreeNode<double, size_t, Capacity, Strategy>;
    using Pt = QPoint<double, size_t>;

    size_t split = pts.size() * 8 / 10;

    vector<Pt> bulk;
    bulk.reserve(split);
    for (size_t i = 0; i < split; ++i) {
        Pt p(i, pts[i][0], pts[i][1]);
        bulk.push_back(p);
    }

    // Bulk load
    Tree tree(bulk.begin(), bulk.end());
    cout << "Bulk loaded " << split << " points\n";

    // Insert
    for (size_t i = split; i < pts.size(); ++i) {
        Pt p(i, pts[i][0], pts[i][1]);
        tree.insert(p);
    }
    cout << "Inserted " << (pts.size() - split) << " more points\n";

    // Delete
    size_t del_count = 0;
    for (size_t i = 0; i < 100; ++i) {
        Pt p(i, pts[i][0], pts[i][1]);
        bool removed = tree.remove(p);
        if (removed)
            ++del_count;
    }
    cout << "Deleted " << del_count << " points\n";

    // Range query
    QBox<double, size_t> qbox(-301.278286, 2883.648480, -221.070077, 2935.126180, 0);
    vector<size_t> results;
    tree.query(qbox, back_inserter(results));
    cout << "Range query returned " << results.size() << " results\n";

    // kNN query
    array<double, Dims> qpt{539.271901, 2915.970570};
    vector<size_t> knn;
    tree.knnQuery(qpt, 10, back_inserter(knn));
    cout << "10-NN query returned " << knn.size() << " results\n";

    // Statistics
    auto s = tree.getStatistics();
    cout << "\n=== Quad-tree (Point, 2D) statistics ===\n";
    cout << "  height          : " << s.height << "\n";
    cout << "  num points      : " << s.numPoints << "\n";
    cout << "  num leaves      : " << s.numLeaves << "\n";
    cout << "  num internals   : " << s.numInternalNodes << "\n";
    cout << "  re-roots        : " << s.numReroots << "\n";
    cout << "  memory (bytes)  : " << s.sizeBytes << "\n";
}
