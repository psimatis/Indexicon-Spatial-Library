#include <iostream>
#include <fstream>
#include <vector>
#include "../indexes/mxcif_quadtree.hpp"

using namespace std;

int main() {
    // Load data
    ifstream fin("../data/input/sample_data.txt");
    if (!fin) { 
        cerr << "Cannot open sample_data.txt\n"; 
        return 1; 
    }

    // Each line: x0 y0 z0 x1 y1 z1  — project to 2-D: (x0, y0, x1, y1)
    vector<array<double,4>> recs;
    double x0, y0, x1, y1, dummy;
    while (fin >> x0 >> y0 >> dummy >> x1 >> y1 >> dummy)
        recs.push_back({ x0, y0, x1, y1 });
    cout << "Loaded " << recs.size() << " MBRs\n";

    using Box  = MXCIFBox<double, size_t>;
    using Tree = MXCIFQuadTreeNode<double, size_t>;

    size_t split = recs.size() * 8 / 10;

    // Bulk load
    vector<Box> bulk;
    bulk.reserve(split);
    for (size_t i = 0; i < split; ++i)
        bulk.push_back(Box(i, recs[i][0], recs[i][1], recs[i][2], recs[i][3]));

    Tree tree(bulk.begin(), bulk.end());
    cout << "Bulk loaded " << split << " MBRs\n";

    // Insert
    for (size_t i = split; i < recs.size(); ++i)
        tree.insert(Box(i, recs[i][0], recs[i][1], recs[i][2], recs[i][3]));
    cout << "Inserted " << (recs.size() - split) << " more MBRs\n";

    // Delete
    size_t del_count = 0;
    for (size_t i = 0; i < 100; ++i)
        if (tree.remove(Box(i, recs[i][0], recs[i][1], recs[i][2], recs[i][3]))) ++del_count;
    cout << "Deleted " << del_count << " MBRs\n";

    // Range query
    Box qbox(0, -301.278286, 2883.648480, -221.070077, 2935.126180);
    vector<size_t> results;
    tree.query(qbox, back_inserter(results));
    cout << "Range query returned " << results.size() << " results\n";

    // kNN query
    vector<size_t> knn;
    tree.knnQuery({539.271901, 2915.970570}, 10, back_inserter(knn));
    cout << "10-NN query returned " << knn.size() << " results\n";

    // Statistics
    auto s = tree.getStatistics();
    cout << "\n=== MX-CIF Quad-tree (2D) statistics ===\n";
    cout << "  height          : " << s.height           << "\n";
    cout << "  num MBRs        : " << s.numPoints        << "\n";
    cout << "  num leaves      : " << s.numLeaves        << "\n";
    cout << "  num internals   : " << s.numInternalNodes << "\n";
    cout << "  memory (bytes)  : " << s.sizeBytes        << "\n";
}
