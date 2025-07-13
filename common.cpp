#include "common.h"
#include <algorithm>
#include <cmath>
#include <cstddef> // For offsetof
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <unordered_set>

void create_mpi_point_type(MPI_Datatype *mpi_point_type) {
    int block_lengths[2] = {1, 1}; // One double per field
    MPI_Aint offsets[2];           // Byte offsets of each field
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE}; // Field types

    // Compute offsets of struct members
    offsets[0] = offsetof(point, x);
    offsets[1] = offsetof(point, y);

    // Create the custom MPI struct datatype
    MPI_Type_create_struct(2, block_lengths, offsets, types, mpi_point_type);
    MPI_Type_commit(mpi_point_type);
}

std::vector<point> read_points(std::istream &file) {
    std::vector<point> points;
    std::string line;

    std::getline(file, line); // Skip first line

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        point p;
        char comma;

        // Read x and y, assuming CSV format
        if (ss >> p.x >> comma >> p.y) {
            points.emplace_back(p);
        }
    }

    return points;
}

std::vector<std::vector<size_t>> read_clusters(std::istream &file) {
    std::vector<std::vector<size_t>> clusters;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        size_t idx;
        char delimiter;
        std::vector<size_t> vals = std::vector<size_t>();

        while (ss >> idx) {
            ss >> delimiter;
            vals.emplace_back(idx);
        }
        clusters.emplace_back(vals);
    }

    return clusters;
}

std::vector<bool> read_annealing_output(std::istream &file) {
    std::vector<bool> solution;
    std::string line;

    std::getline(file, line);

    std::istringstream ss(line);
    bool active;
    char delimiter;
    std::vector<size_t> vals = std::vector<size_t>();

    while (ss >> active) {
        ss >> delimiter;
        solution.emplace_back(active);
    }

    return solution;
}

int assign_outliers(std::map<size_t, std::vector<point>> &clusters,
                    const std::vector<point> &outliers) {
    if (clusters.empty()) {
        std::cerr << "Cannot add outliers to empty clusters" << std::endl;
        return EXIT_FAILURE;
    }

    // Find the centroid of each cluster
    std::map<size_t, point> centroids;
    for (auto &[k, v] : clusters) {
        auto x = std::accumulate(
                     v.begin(), v.end(), 0.0,
                     [](auto accum, const auto &pt) { return accum + pt.x; }) /
                 v.size();
        auto y = std::accumulate(
                     v.begin(), v.end(), 0.0,
                     [](auto accum, const auto &pt) { return accum + pt.y; }) /
                 v.size();
        centroids[k] = point{x, y};
    }

    // Assign each outlier to the closest existing cluster
    for (const auto &outlier : outliers) {
        std::map<size_t, double> distances;
        for (const auto &[k, centroid] : centroids) {
            distances[k] = euclidean_distance(outlier, centroid);
        }
        auto smallest_distance = std::min_element(
            distances.begin(), distances.end(),
            [](const auto &l, const auto &r) { return l.second < r.second; });
        clusters[smallest_distance->first].emplace_back(outlier);
    }

    return EXIT_SUCCESS;
}

double euclidean_distance(const point &a, const point &b) {
    return std::hypot(a.x - b.x, a.y - b.y);
}

double silhouette(std::map<size_t, std::vector<point>> &clusters) {
    if (clusters.size() < 2) {
        std::cerr
            << "Silhouette score is ill-defined when using less than 2 clusters"
            << std::endl;
        return -1.0;
    }

    auto s_values = std::vector<double>();

    for (auto const &[this_label, this_points] : clusters) {
        auto cluster_size = this_points.size();
        for (auto const i : this_points) {
            double a, b;
            // Dissimilarities within the same cluster
            if (cluster_size == 1) {
                // a(i) can be defined as 0 when a cluster contains a single
                // element
                a = 0;
            } else {
                double num = 0;
                // We do not need to skip the distance between i and i itself
                // because it is 0
                for (auto const j : this_points) {
                    num += euclidean_distance(i, j);
                }
                a = num / (cluster_size - 1);
            }

            // Dissimilarities to other clusters
            auto distances = std::vector<double>();
            for (auto const &[other_label, other_points] : clusters) {
                if (other_label == this_label)
                    continue;

                double num = 0;
                for (auto j : other_points) {
                    num += euclidean_distance(i, j);
                }
                distances.emplace_back(num / other_points.size());
            }
            b = *std::min_element(distances.begin(), distances.end());

            double s = (b - a) / std::max(a, b);
            s_values.emplace_back(s);
        }
    }

    return std::accumulate(s_values.begin(), s_values.end(), 0.0) /
           s_values.size();
}

/*
Compute Silhoutte score using Euclidean distance as the distance metric:
https://en.wikipedia.org/wiki/Silhouette_(clustering)
All points are assumed to be assigned to a cluster.
*/
double silhouette(const std::vector<point> &pts,
                  const std::vector<int> &labels) {
    if (pts.empty()) {
        std::cerr << "Cannot compute Silhoutte score from empty data"
                  << std::endl;
        return -1.0;
    }

    auto clusters = std::map<size_t, std::vector<point>>();

    for (auto i = 0; i < pts.size(); i++) {
        auto label = labels[i];
        if (!clusters.contains(label)) {
            clusters[label] = std::vector<point>();
        }
        clusters[label].emplace_back(pts[i]);
    }
    return silhouette(clusters);
}

/*
Read input points file, intermediate clustering file and quantum job output to
reconstruct the corresponding clusters.
The resulting clusters and the unclassified points are inserted into the last
two arguments.
*/
int parse_quantum_job_output(
    const std::string &points_name, const std::string &clusters_name,
    const std::string &quantum_job_output_name,
    std::map<size_t, std::vector<point>> &solution_clusters,
    std::vector<point> &outliers) {

    std::ifstream file(points_name);
    if (!file) {
        std::cerr << "Error opening points file." << std::endl;
        return EXIT_FAILURE;
    }

    auto pts = read_points(file);

    std::ifstream clusters_file(clusters_name);
    if (!clusters_file) {
        std::cerr << "Error opening clusters file." << std::endl;
        return EXIT_FAILURE;
    }

    auto all_clusters = read_clusters(clusters_file);

    std::ifstream quantum_job_file(quantum_job_output_name);
    if (!quantum_job_file) {
        std::cerr << "Error opening quantum output file." << std::endl;
        return EXIT_FAILURE;
    }

    auto solution = read_annealing_output(quantum_job_file);

    std::unordered_set<size_t> outliers_idx;
    for (auto i = 0; i < pts.size(); i++) {
        outliers_idx.insert(i);
    }

    // Get indices of selected clusters
    for (size_t i = 0; i < solution.size(); i++) {
        if (solution[i]) {
            auto this_cluster_points = std::vector<point>();
            for (const auto pts_idx : all_clusters[i]) {
                this_cluster_points.emplace_back(pts[pts_idx]);
                outliers_idx.erase(pts_idx);
            }
            solution_clusters[i] = this_cluster_points;
        }
    }

    for (auto i : outliers_idx) {
        outliers.emplace_back(pts[i]);
    }

    return EXIT_SUCCESS;
}
