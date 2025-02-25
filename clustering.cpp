#include <mpi.h>
#include <numeric>
#include <set>
#include <span>
#include <stdio.h>
#include <stdlib.h>

#include "points.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#define MAX_LINE_LENGTH 1024

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

int main(int argc, char **argv) {
    int myid, nproc;
    MPI_Status status;

    char line[MAX_LINE_LENGTH];
    int row = 0;

    int nptsincluster;
    std::vector<point> pts;

    double z1, z2;

    int nmethod_proc;
    std::vector<std::string> method = {"kmeans", "dbscan", "hclust"};
    int nmethod = method.size();

    int ncl_tot;

    int (*functions[])(int, const char *, point *, int,
                       int *) = {kmeans, dbscan, hclust};

    int res;

    MPI_Datatype MPI_POINT;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    create_mpi_point_type(&MPI_POINT);

    if (nmethod % nproc != 0) {
        printf("Error");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    nmethod_proc = nmethod / nproc;

    // Read input file in rank 0
    if (myid == 0) {
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

        std::string line;
        std::getline(file, line); // Skip first line

        while (std::getline(file, line)) {
            std::istringstream ss(line);
            point p;
            char comma;

            // Read x and y, assuming CSV format
            if (ss >> p.x >> comma >> p.y) {
                pts.emplace_back(p);
            }
        }
        nptsincluster = pts.size();
    }

    // Broadcast parsed input to other ranks
    MPI_Bcast(&nptsincluster, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (myid != 0) {
        pts.resize(nptsincluster);
    }

    //   MPI_Bcast(pts, nptsincluster * sizeof(point), MPI_BYTE, 0,
    //   MPI_COMM_WORLD);
    MPI_Bcast(pts.data(), nptsincluster, MPI_POINT, 0, MPI_COMM_WORLD);
    // Print received data in each process
    /*   for (int i = 0; i < nptsincluster; i++)
       {
         printf("I am proc %d and in pos %d, I received %f %f\n", myid, i,
       pts[i].x, pts[i].y);
       }
       printf("\n");
    */

    std::vector<int> clus(nptsincluster, 0);
    for (int im = 0; im < nmethod_proc; im++) {
        std::string mymethod = method[nmethod_proc * myid + im];
        printf("I am proc %d and I will deal with method %s\n", myid,
               mymethod.c_str());
        for (int i = 0; i < nmethod; i++) {
            if (mymethod == method[i]) {
                res = functions[i](myid, mymethod.c_str(), pts.data(),
                                   nptsincluster, clus.data());
            }
        }
    }

    int max = *std::max_element(clus.begin(), clus.end());
    int min = *std::min_element(clus.begin(), clus.end());
    int ncl = (max - min) + 1;

    std::vector<int> all_res; // Aggregation of all clustering results
    std::vector<int>
        cluster_counts; // The number of clusters from each clustering algorithm
    std::vector<int> offsets; // The offset that should be applied to each
                              // cluster, depending on its algorithm
    if (myid == 0) {
        all_res.resize(nptsincluster * nproc);
        cluster_counts.resize(nproc);
        offsets.resize(nproc);
    }

    MPI_Gather(&ncl, 1, MPI_INT, cluster_counts.data(), 1, MPI_INT, 0,
               MPI_COMM_WORLD);

    // Each rank should adjust the indices assigned to its clusters,
    // preventing clusters from different algorithms to have identical IDs.
    if (myid == 0) {
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

    MPI_Gather(clus.data(), nptsincluster, MPI_INT, all_res.data(),
               nptsincluster, MPI_INT, 0, MPI_COMM_WORLD);

    if (myid == 0) {
        printf("Total number of clusters is %d\n", ncl_tot);

        std::vector<std::set<int>> cluster_elems(ncl_tot);

        // Create sets from each clustering algorithm
        for (size_t i = 0; i < nproc; i++) {
            std::span method_data{all_res.begin() + i * nptsincluster,
                                  all_res.begin() + (i + 1) * nptsincluster};

            for (size_t index = 0; index < nptsincluster; index++) {
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
