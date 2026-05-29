#include <iostream>
#include <fstream>
#include <vector>
#include "../indexes/rtree_mbr.hpp"

using namespace std;

int main() {
    // Load data
    ifstream fin("../data/input/sample_data.txt");
    if (!fin) { cerr << "Cannot open sample_data.txt\n"; return 1; }

    // Each line: x0 y0 z0 x1 y1 z1  (3-D MBB)
    vector<array<double,6>> recs;
    double x0, y0, z0, x1, y1, z1;
    while (fin >> x0 >> y0 >> z0 >> x1 >> y1 >> z1)
        recs.push_back({ x0, y0, z0, x1, y1, z1 });
    cout << "Loaded " << recs.size() << " MBRs\n";

    using Tree = RTree<size_t, double, 3>;
    using Box  = Tree::BoxType;
    using Pt   = Tree::PointType;

    size_t split = recs.size() * 8 / 10;

    // Bulk load
    vector<pair<Box,size_t>> bulk;
    bulk.reserve(split);
    for (size_t i = 0; i < split; ++i) {
        Box b;
        b.minCorner = Pt(recs[i][0], recs[i][1], recs[i][2]);
        b.maxCorner = Pt(recs[i][3], recs[i][4], recs[i][5]);
        bulk.push_back({ b, i });
    }

    Tree tree(bulk.begin(), bulk.end());
    cout << "Bulk loaded " << split << " MBRs\n";

    // Insert
    for (size_t i = split; i < recs.size(); ++i) {
        Box b;
        b.minCorner = Pt(recs[i][0], recs[i][1], recs[i][2]);
        b.maxCorner = Pt(recs[i][3], recs[i][4], recs[i][5]);
        tree.insert(b, i);
    }
    cout << "Inserted " << (recs.size() - split) << " more MBRs\n";

    // Delete
    size_t del_count = 0;
    for (size_t i = 0; i < 100; ++i) {
        Box b;
        b.minCorner = Pt(recs[i][0], recs[i][1], recs[i][2]);
        b.maxCorner = Pt(recs[i][3], recs[i][4], recs[i][5]);
        if (tree.remove(b, i)) ++del_count;
    }
    cout << "Deleted " << del_count << " MBRs\n";

    // Range query
    Box qbox;
    qbox.minCorner = Pt(-301.278286, 2883.648480, -20.416440);
    qbox.maxCorner = Pt(-221.070077, 2935.126180, -20.262942);
    vector<size_t> results;
    tree.query(qbox, back_inserter(results));
    cout << "Range query returned " << results.size() << " results\n";

    // kNN query
    Pt qpt(539.271901, 2915.970570, -21.334010);
    vector<size_t> knn;
    tree.knnQuery(qpt, 10, back_inserter(knn));
    cout << "10-NN query returned " << knn.size() << " results\n";

    // Statistics
    auto s = tree.getStatistics();
    cout << "\n=== R-tree (MBR, 3D) statistics ===\n";
    cout << "  height          : " << s.height           << "\n";
    cout << "  num MBRs        : " << s.numPoints        << "\n";
    cout << "  num leaves      : " << s.numLeaves        << "\n";
    cout << "  num internals   : " << s.numInternalNodes << "\n";
    cout << "  memory (bytes)  : " << s.sizeBytes        << "\n";
}
