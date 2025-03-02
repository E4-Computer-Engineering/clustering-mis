#include <mpi.h>

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

int kmeans(int myrank, const char *str, const point *pts, int np, int *res);
#ifdef __cplusplus
}
#endif

int dbscan(int myrank, const char *str, const point *pts, int np, int *res_out);
int hclust(int myrank, const char *str, const point *pts, int np, int *res_out);

#endif // POINT_H
