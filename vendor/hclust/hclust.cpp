#include "common.h"
#include "fastcluster.hpp"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <random>
#include <vector>

double distance(const point &a, const point &b) { return std::hypot(a.x - b.x, a.y - b.y); }

std::vector<double> create_distance_matrix(const point *data, size_t n) {
    std::vector<double> distmat((n * (n - 1)) / 2);

    int k, i, j;
    for (i = k = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            // compute distance between observables i and j
            distmat[k] = distance(data[i], data[j]);
            k++;
        }
    }

    return distmat;
}

// TODO move to another file, use a shared interface for both k-means and DBSCAN
int hclust(int myrank, const char *str, const point *pts, int np, int *res, int seed) {
    static std::default_random_engine eng;
    eng.seed(seed + 0xE4 * 2);
    auto min_k = 4;
    auto max_k = 10;
    static std::uniform_int_distribution<> dis(min_k, max_k); // range [min_k, max_k]
    size_t NUM_CLUSTERS = dis(eng);

    auto distmat = create_distance_matrix(pts, np);

    int *merge = new int[2 * (np - 1)];
    double *height = new double[np - 1];

    hclust_fast(np, distmat.data(), HCLUST_METHOD_SINGLE, merge, height);

    // partitioning into nclust clusters
    cutree_k(np, merge, NUM_CLUSTERS, res);
    //   // stop clustering at step with cluster distance >= cdist
    //   cutree_cdist(npoints, merge, height, cdist, labels);

    // Clean up allocated variables
    delete[] merge;
    delete[] height;

    bool LOG_TO_FILE = false;
    if (LOG_TO_FILE) {
        std::ofstream out_file("hclust_output.txt");
        out_file << "Feature 1\tFeature 2\tCluster" << std::endl;
        for (size_t i = 0; i < np; i++) {
            out_file << pts[i].x << "\t" << pts[i].y << "\t" << res[i] << std::endl;
        }
        out_file.close();
    }

    return 0;
}
