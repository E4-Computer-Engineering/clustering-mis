#include <cmath>
#include <cstdlib>
#include <functional>
#include <mpi.h>
#include <numeric>
#include <ostream>
#include <set>
#include <span>
#include <stdlib.h>

#include "points.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using ClusFuncType =
    std::function<int(int, const char *, const point *, int, int *)>;

// Short-circuited sets intersection
bool have_shared_elem(const std::set<int> &x, const std::set<int> &y) {
    auto i = x.begin();
    auto j = y.begin();
    while (i != x.end() && j != y.end()) {
        if (*i == *j)
            return true;
        else if (*i < *j)
            i++;
        else
            j++;
    }
    return false;
}

std::vector<std::vector<double>>
create_overlap_matrix(const std::vector<std::set<int>> &clusters) {
    // We want to give a different weight to each diagonal term. Otherwise, if
    // all clusters are considered equally good, the annealer will choose a
    // solution with many small clusters.
    std::vector<size_t> sizes(clusters.size());
    std::transform(clusters.begin(), clusters.end(), sizes.begin(),
                   [](const auto &cl) { return cl.size(); });
    // The biggest cluster will have a weight of 1. The others will be
    // normalized to be smaller, but in the range 0 < x < 1.
    auto max_size = *std::max_element(sizes.begin(), sizes.end());

    auto n = clusters.size();
    auto penalty = n;

    // Initialize empty matrix
    std::vector<std::vector<double>> res(n, std::vector<double>(n, 0.0));

    // Add penalty to overlapping clusters
    for (auto i = 0; i < n - 1; i++) {
        for (auto j = i + 1; j < n; j++) {
            if (have_shared_elem(clusters[i], clusters[j])) {
                res[i][j] = penalty;
            }
        }
    }

    // Set diagonal terms
    for (auto i = 0; i < n; i++) {
        res[i][i] = -(double)sizes[i] / max_size;
    }
    return res;
}

void write_matrix(std::ostream &out_stream,
                  const std::vector<std::vector<double>> &m) {
    for (auto i : m) {
        for (auto j = i.begin(); j != i.end(); j++) {
            if (j != i.begin()) {
                out_stream << " ";
            }
            out_stream << *j;
        }
        out_stream << std::endl;
    }
}

void print_matrix(const std::vector<std::vector<double>> &m) {
    write_matrix(std::cout, m);
}

int save_matrix(const std::string &file_name,
                const std::vector<std::vector<double>> &m) {
    std::ofstream file(file_name);
    if (!file) {
        std::cerr << "Error opening file." << std::endl;
        return EXIT_FAILURE;
    } else {
        write_matrix(file, m);
        return EXIT_SUCCESS;
    }
}

int save_clusters(const std::string &file_name,
                  std::vector<std::set<int>> cluster_elems) {
    std::ofstream file(file_name);
    if (!file) {
        std::cerr << "Error opening file." << std::endl;
        return EXIT_FAILURE;
    } else {
        for (auto cluster : cluster_elems) {
            for (auto j = cluster.begin(); j != cluster.end(); j++) {
                if (j != cluster.begin()) {
                    file << ",";
                }
                file << *j;
            }
            file << std::endl;
        }
        return EXIT_SUCCESS;
    }
}

void read_points(std::istream &file, std::vector<point> &points) {
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
}

// Return total number of clusters identified by each algorithm
int run_clustering_algorithms(int my_rank, int num_methods_proc,
                              const std::vector<point> &pts, int num_methods,
                              const std::vector<std::string> &methods,
                              std::vector<ClusFuncType> &functions,
                              std::vector<int> &assigned_clusters) {
    auto num_points = pts.size();

    // Each process can run more than one clustering algorithm.
    // This offset is used to track previous results and
    // to make every cluster number computed by this process unique.
    // e.g. Cluster 0 from algorithm 1 should not have the same label of
    // cluster 0 from algorithm 2.
    int clustering_offset = 0;

    for (int local_method_idx = 0; local_method_idx < num_methods_proc;
         local_method_idx++) {
        auto global_method_idx = num_methods_proc * my_rank + local_method_idx;

        if (global_method_idx >= num_methods) {
            break; // This method does not exist
        }

        std::string current_method_name = methods.at(global_method_idx);
        std::cout << "I am proc " << my_rank << " and I will deal with method "
                  << current_method_name << std::endl;

        auto res = functions[global_method_idx](
            my_rank, current_method_name.c_str(), pts.data(), num_points,
            assigned_clusters.data() + local_method_idx * num_points);

        auto begin = assigned_clusters.begin() + local_method_idx * num_points;
        auto end = begin + num_points;

        int max = *std::max_element(begin, end);
        int min = *std::min_element(begin, end);
        int current_num_cluster = (max - min) + 1;

        // If clustering indices started from e.g. 1, we should force them to
        // start from zero instead. We can add this adjustment term on top of
        // the other offset.
        int actual_offset = clustering_offset - min;
        for (auto it = begin; it != end; it++) {
            *it = *it + actual_offset;
        }

        clustering_offset += current_num_cluster;
    }
    return clustering_offset;
}

int main(int argc, char **argv) {
    int my_rank, num_processes;

    std::vector<std::string> methods = {"kmeans", "dbscan", "hclust"};
    std::vector<ClusFuncType> functions = {kmeans, dbscan, hclust};
    int num_methods = methods.size();

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Datatype MPI_POINT;
    create_mpi_point_type(&MPI_POINT);

    std::vector<point> pts;
    int num_points;
    // Read input file in rank 0
    if (my_rank == 0) {
        if (argc <= 1) {
            std::cerr << "No input file was specified." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        std::ifstream file(argv[1]);
        if (!file) {
            std::cerr << "Error opening file." << std::endl;
            ;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        read_points(file, pts);
        num_points = pts.size();
    }

    // Broadcast parsed input to other ranks
    MPI_Bcast(&num_points, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (my_rank != 0) {
        pts.resize(num_points);
    }

    MPI_Bcast(pts.data(), num_points, MPI_POINT, 0, MPI_COMM_WORLD);

    // The maximum number of methods to be run by each process
    auto num_methods_proc =
        std::ceil((float)num_methods / (float)num_processes);

    // Run the algorithms assigned to this process and flatten their results
    std::vector<int> assigned_clusters(num_methods_proc * num_points, 0);
    int ncl =
        run_clustering_algorithms(my_rank, num_methods_proc, pts, num_methods,
                                  methods, functions, assigned_clusters);

    std::vector<int> all_res;        // Aggregation of all clustering results
                                     // across all processes
    std::vector<int> cluster_counts; // The number of clusters from each rank
    std::vector<int> offsets; // The offset that should be applied to each
                              // cluster, depending on its rank
    if (my_rank == 0) {
        all_res.resize(num_points * num_processes * num_methods_proc);
        cluster_counts.resize(num_processes);
        offsets.resize(num_processes);
    }

    MPI_Gather(&ncl, 1, MPI_INT, cluster_counts.data(), 1, MPI_INT, 0,
               MPI_COMM_WORLD);

    int ncl_tot;
    // Each rank should adjust the indices assigned to its clusters,
    // preventing clusters from different algorithms to have identical IDs.
    if (my_rank == 0) {
        std::partial_sum(cluster_counts.begin(), cluster_counts.end(),
                         offsets.begin());
        ncl_tot = offsets.back();
        offsets.pop_back();
        offsets.insert(offsets.begin(), 0);
    }

    int indices_offset;
    MPI_Scatter(offsets.data(), 1, MPI_INT, &indices_offset, 1, MPI_INT, 0,
                MPI_COMM_WORLD);

    for (auto it = assigned_clusters.begin(); it != assigned_clusters.end();
         it++) {
        *it = *it + indices_offset;
    }

    MPI_Gather(assigned_clusters.data(), num_points * num_methods_proc, MPI_INT,
               all_res.data(), num_points * num_methods_proc, MPI_INT, 0,
               MPI_COMM_WORLD);

    if (my_rank == 0) {
        std::cout << "Total number of clusters is " << ncl_tot << std::endl;

        std::vector<std::set<int>> cluster_elems(ncl_tot);

        // Handle cases where the number of methods is not a multiple of the
        // number of processes (this removes the trailing zeros that can be seen
        // in all_res in such cases, that would otherwise clash with other zeros
        // corresponding to the 0-th cluster)
        all_res.resize(num_points * num_methods);

        // Create sets from each clustering algorithm
        for (size_t i = 0; i < num_methods; i++) {
            std::span method_data{all_res.begin() + i * num_points,
                                  all_res.begin() + (i + 1) * num_points};

            for (size_t index = 0; index < num_points; index++) {
                auto assigned_cluster = method_data[index];
                cluster_elems[assigned_cluster].insert(index);
            }
        }
        auto overlap_matrix = create_overlap_matrix(cluster_elems);

        if (argc < 3) {
            std::cout << "Overlap matrix:" << std::endl;
            print_matrix(overlap_matrix);
        } else {
            auto status = save_matrix(argv[2], overlap_matrix);
            if (status == EXIT_FAILURE) {
                std::cerr << "Unable to write the overlap matrix to file."
                          << std::endl;
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            if (argc > 3) {
                auto clus_status = save_clusters(argv[3], cluster_elems);
                if (clus_status == EXIT_FAILURE) {
                    std::cerr
                        << "Unable to write the obtained clusters to file."
                        << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                }
            }
        }
    }

    MPI_Finalize();
    return 0;
}
