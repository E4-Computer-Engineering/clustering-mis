#include <istream>
#include <map>
#include <mpi.h>
#include <vector>

#ifndef POINT_H
#define POINT_H

// Define the structure
typedef struct {
    double x;
    double y;
} point;

void create_mpi_point_type(MPI_Datatype *mpi_point_type);

#ifdef __cplusplus
extern "C" {
#endif

int kmeans(int myrank, const char *str, const point *pts, int np, int *res, int seed);
#ifdef __cplusplus
}
#endif

int dbscan(int myrank, const char *str, const point *pts, int np, int *res_out, int seed);
int hclust(int myrank, const char *str, const point *pts, int np, int *res_out, int seed);

std::vector<point> read_points(std::istream &file);
std::vector<std::vector<size_t>> read_clusters(std::istream &file);
std::vector<bool> read_annealing_output(std::istream &file);
int assign_outliers(std::map<size_t, std::vector<point>> &clusters,
                    const std::vector<point> &outliers);
double euclidean_distance(const point &a, const point &b);

double silhouette(std::map<size_t, std::vector<point>> &clusters);

/*
Compute Silhoutte score using Euclidean distance as the distance metric:
https://en.wikipedia.org/wiki/Silhouette_(clustering)
All points are assumed to be assigned to a cluster.
*/
double silhouette(const std::vector<point> &pts,
                  const std::vector<int> &labels);

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
    std::vector<point> &outliers);

#endif // POINT_H
