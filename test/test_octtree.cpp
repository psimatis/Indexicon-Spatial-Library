#include <iostream>
#include <fstream>
#include <vector>
#include "../indexes/octtree.hpp"

using namespace std;

int main() {
    // Load data
    ifstream fin("../data/input/sample_data.txt");
    if (!fin) { cerr << "Cannot open sample_data.txt\n"; return 1; }

    // Each line: x y z ...  — read first 3 values, ignore the rest
    vector<array<double,3>> pts;
    double x, y, z, dummy;
    while (fin >> x >> y >> z >> dummy >> dummy >> dummy)
        pts.push_back({ x, y, z });
    cout << "Loaded " << pts.size() << " points\n";

    using Tree = OctreeNode<double, size_t>;
    using Pt   = OPoint<double, size_t>;

    size_t split = pts.size() * 8 / 10;

    vector<Pt> bulk;
    bulk.reserve(split);
    for (size_t i = 0; i < split; ++i)
        bulk.push_back(Pt(i, pts[i][0], pts[i][1], pts[i][2]));

    // Bulk load
    Tree tree(bulk.begin(), bulk.end());
    cout << "Bulk loaded " << split << " points\n";

    // Insert
    for (size_t i = split; i < pts.size(); ++i)
        tree.insert(Pt(i, pts[i][0], pts[i][1], pts[i][2]));
    cout << "Inserted " << (pts.size() - split) << " more points\n";

    // Delete
    size_t del_count = 0;
    for (size_t i = 0; i < 100; ++i)
        if (tree.remove(Pt(i, pts[i][0], pts[i][1], pts[i][2]))) ++del_count;
    cout << "Deleted " << del_count << " points\n";

    // Range query
    OBox<double,size_t> qbox(-400.0, 2800.0, -22.0, 0.0, 3100.0, -19.0);
    vector<size_t> results;
    tree.query(qbox, back_inserter(results));
    cout << "Range query returned " << results.size() << " results\n";

    // kNN query
    vector<size_t> knn;
    tree.knnQuery({539.271901, 2915.970570, -21.334010}, 10, back_inserter(knn));
    cout << "10-NN query returned " << knn.size() << " results\n";

    // Statistics
    auto s = tree.getStatistics();
    cout << "\n=== Oct-tree (3D) statistics ===\n";
    cout << "  height          : " << s.height           << "\n";
    cout << "  num points      : " << s.numPoints        << "\n";
    cout << "  num leaves      : " << s.numLeaves        << "\n";
    cout << "  num internals   : " << s.numInternalNodes << "\n";
    cout << "  memory (bytes)  : " << s.sizeBytes        << "\n";
}
