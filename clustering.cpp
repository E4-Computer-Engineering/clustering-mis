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

std::vector<std::vector<int>>
create_overlap_matrix(const std::vector<std::set<int>> &clusters) {
    auto n = clusters.size();
    auto penalty = n;

    // Initialize empty matrix
    std::vector<std::vector<int>> res(n, std::vector<int>(n, 0));

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
        res[i][i] = -1;
    }
    return res;
}

void print_matrix(const std::vector<std::vector<int>> &m) {
    for (auto i : m) {
        for (auto j : i) {
            std::cout << j << " ";
        }
        std::cout << std::endl;
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

int main(int argc, char **argv) {
    int my_rank, num_processes;

    std::vector<std::string> methods = {"kmeans", "dbscan", "hclust"};
    int num_methods = methods.size();

    int (*functions[])(int, const char *, point *, int,
                       int *) = {kmeans, dbscan, hclust};

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Datatype MPI_POINT;
    create_mpi_point_type(&MPI_POINT);

    // FIXME this check suggests that a rank can compute multiple methods,
    // but the data is currently overwritten instead of accumulated
    if (num_methods % num_processes != 0) {
        std::cout << "Aborting, the number of clustering methods should be a "
                     "multiple of the number of processes."
                  << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    auto num_methods_proc = num_methods / num_processes;

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

    std::vector<int> clus(num_points, 0);
    for (int offset = 0; offset < num_methods_proc; offset++) {
        auto method_idx = num_methods_proc * my_rank + offset;

        std::string current_method_name = methods.at(method_idx);
        std::cout << "I am proc " << my_rank << " and I will deal with method "
                  << current_method_name << std::endl;
        auto res = functions[method_idx](my_rank, current_method_name.c_str(),
                                         pts.data(), num_points, clus.data());
    }

    int max = *std::max_element(clus.begin(), clus.end());
    int min = *std::min_element(clus.begin(), clus.end());
    int ncl = (max - min) + 1;

    std::vector<int> all_res; // Aggregation of all clustering results
    std::vector<int>
        cluster_counts; // The number of clusters from each clustering algorithm
    std::vector<int> offsets; // The offset that should be applied to each
                              // cluster, depending on its algorithm
    if (my_rank == 0) {
        all_res.resize(num_points * num_processes);
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
    // If clustering indices started from e.g. 1, we should force them to start
    // from zero instead. We can add this adjustment term on top of the global
    // offset.
    int actual_offset = indices_offset - min;

    for (auto it = clus.begin(); it != clus.end(); it++) {
        *it = *it + actual_offset;
    }

    MPI_Gather(clus.data(), num_points, MPI_INT, all_res.data(), num_points,
               MPI_INT, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        std::cout << "Total number of clusters is " << ncl_tot << std::endl;

        std::vector<std::set<int>> cluster_elems(ncl_tot);

        // Create sets from each clustering algorithm
        for (size_t i = 0; i < num_processes; i++) {
            std::span method_data{all_res.begin() + i * num_points,
                                  all_res.begin() + (i + 1) * num_points};

            for (size_t index = 0; index < num_points; index++) {
                auto assigned_cluster = method_data[index];
                cluster_elems[assigned_cluster].insert(index);
            }
        }
        auto overlap_matrix = create_overlap_matrix(cluster_elems);
        std::cout << "Overlap matrix:" << std::endl;
        print_matrix(overlap_matrix);
    }

    MPI_Finalize();
    return 0;
}
