#include <array>
#include <fstream>
#include <iostream>
#include <vector>
#include "../../indexes/kdtree.hpp"

using namespace std;

int main() {
    constexpr size_t Dims = 3;
    constexpr size_t Capacity = 128;
    constexpr KDSplit SplitMethod = KDSplit::ADAPTIVE;
    // Other split methods: KDSplit::ROUND_ROBIN, KDSplit::LONGEST_AXIS

    // Load data
    ifstream fin("../data/input/sample_data.txt");
    if (!fin) {
        cerr << "Cannot open sample_data.txt\n";
        return 1;
    }

    // Each line: x y z ...; read first 3 values, ignore the rest
    vector<array<double, Dims>> pts;
    double x, y, z, dummy;
    while (fin >> x >> y >> z >> dummy >> dummy >> dummy)
        pts.push_back({x, y, z});
    cout << "Loaded " << pts.size() << " points\n";

    using Tree = KDTree<double, size_t, Dims, Capacity, SplitMethod>;
    using Pt = Tree::Point;
    using Reg = Tree::Region;

    size_t split = pts.size() * 8 / 10;

    vector<pair<Pt, size_t>> bulk;
    bulk.reserve(split);
    for (size_t i = 0; i < split; ++i) {
        Pt p(pts[i][0], pts[i][1], pts[i][2]);
        bulk.push_back({p, i});
    }

    // Bulk load
    Tree tree(bulk.begin(), bulk.end());
    cout << "Bulk loaded " << split << " points\n";

    // Insert
    for (size_t i = split; i < pts.size(); ++i) {
        Pt p(pts[i][0], pts[i][1], pts[i][2]);
        tree.insert(p, i);
    }
    cout << "Inserted " << (pts.size() - split) << " more points\n";

    // Delete
    size_t del_count = 0;
    for (size_t i = 0; i < 100; ++i) {
        Pt p(pts[i][0], pts[i][1], pts[i][2]);
        bool removed = tree.remove(p);
        if (removed)
            ++del_count;
    }
    cout << "Deleted " << del_count << " points\n";

    // Range query
    Reg q;
    q.low = Pt(-400.0, 2800.0, -22.0);
    q.high = Pt(0.0, 3100.0, -19.0);
    vector<size_t> results;
    tree.query(q, back_inserter(results));
    cout << "Range query returned " << results.size() << " results\n";

    // kNN query
    Pt qpt(539.271901, 2915.970570, -21.334010);
    vector<size_t> knn;
    tree.knnQuery(qpt, 10, back_inserter(knn));
    cout << "10-NN query returned " << knn.size() << " results\n";

    // Statistics
    auto s = tree.getStatistics();
    cout << "\n=== KD-tree (Point, 3D) statistics ===\n";
    cout << "  height          : " << s.height << "\n";
    cout << "  num points      : " << s.numPoints << "\n";
    cout << "  num leaves      : " << s.numLeaves << "\n";
    cout << "  num internals   : " << s.numInternalNodes << "\n";
    cout << "  memory (bytes)  : " << s.sizeBytes << "\n";
}
