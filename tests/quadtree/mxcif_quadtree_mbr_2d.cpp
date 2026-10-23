#include <array>
#include <fstream>
#include <iostream>
#include <vector>
#include "../../indexes/mxcif_quadtree.hpp"

using namespace std;

int main() {
    constexpr size_t Dims = 2;
    constexpr size_t BoxCoordinates = Dims * 2;
    constexpr size_t Capacity = 128;

    // Load data
    ifstream fin("../data/input/sample_data.txt");
    if (!fin) { 
        cerr << "Cannot open sample_data.txt\n"; 
        return 1; 
    }

    // Each line: x0 y0 z0 x1 y1 z1; project to 2D MBRs
    vector<array<double, BoxCoordinates>> recs;
    double x0, y0, z0, x1, y1, z1;
    while (fin >> x0 >> y0 >> z0 >> x1 >> y1 >> z1)
        recs.push_back({x0, y0, x1, y1});
    cout << "Loaded " << recs.size() << " MBRs\n";

    using Box = MXCIFBox<double, size_t>;
    using Tree = MXCIFQuadTreeNode<double, size_t, Capacity>;

    size_t split = recs.size() * 8 / 10;

    // Bulk load
    vector<Box> bulk;
    bulk.reserve(split);
    for (size_t i = 0; i < split; ++i) {
        Box b(i, recs[i][0], recs[i][1], recs[i][2], recs[i][3]);
        bulk.push_back(b);
    }

    Tree tree(bulk.begin(), bulk.end());
    cout << "Bulk loaded " << split << " MBRs\n";

    // Insert
    for (size_t i = split; i < recs.size(); ++i) {
        Box b(i, recs[i][0], recs[i][1], recs[i][2], recs[i][3]);
        tree.insert(b);
    }
    cout << "Inserted " << (recs.size() - split) << " more MBRs\n";

    // Delete
    size_t del_count = 0;
    for (size_t i = 0; i < 100; ++i) {
        Box b(i, recs[i][0], recs[i][1], recs[i][2], recs[i][3]);
        bool removed = tree.remove(b);
        if (removed)
            ++del_count;
    }
    cout << "Deleted " << del_count << " MBRs\n";

    // Range query
    Box qbox(0, -301.278286, 2883.648480, -221.070077, 2935.126180);
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
    cout << "\n=== MX-CIF Quad-tree (MBR, 2D) statistics ===\n";
    cout << "  height          : " << s.height << "\n";
    cout << "  num MBRs        : " << s.numPoints << "\n";
    cout << "  num leaves      : " << s.numLeaves << "\n";
    cout << "  num internals   : " << s.numInternalNodes << "\n";
    cout << "  memory (bytes)  : " << s.sizeBytes << "\n";
}
