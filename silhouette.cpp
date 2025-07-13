#include "common.h"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <vector>

int main(int argc, char **argv) {
    /*
    We expect these arguments:
    1) Points input file name
    2) Intermediate clustering file name
    3) Quantum job output file name
    4) Silhouette score output file (optional)
    */

    if (argc < 4) {
        std::cerr << "Provide all input files." << std::endl;
        return EXIT_FAILURE;
    }
    std::map<size_t, std::vector<point>> current_clusters;
    std::vector<point> outliers;

    auto status = parse_quantum_job_output(argv[1], argv[2], argv[3],
                                           current_clusters, outliers);

    if (status == EXIT_FAILURE) {
        std::cerr << "Something went wrong while reading files." << std::endl;
        return EXIT_FAILURE;
    }

    // TODO remove outliers management if unnecessary
    // We currently want all points to be assigned to a cluster
    auto outliers_status = assign_outliers(current_clusters, outliers);

    auto current_silhoutte_score = silhouette(current_clusters);

    if (argc == 4) {
        std::cout << current_silhoutte_score << std::endl;
    } else {
        std::ofstream out_file(argv[4]);
        out_file << current_silhoutte_score << std::endl;
    }

    return EXIT_SUCCESS;
}
