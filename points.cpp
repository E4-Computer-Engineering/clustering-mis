#include "points.h"
#include <stddef.h>  // For offsetof

void create_mpi_point_type(MPI_Datatype *mpi_point_type) {
    int block_lengths[2] = {1, 1}; // One double per field
    MPI_Aint offsets[2];  // Byte offsets of each field
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE}; // Field types

    // Compute offsets of struct members
    offsets[0] = offsetof(point, x);
    offsets[1] = offsetof(point, y);

    // Create the custom MPI struct datatype
    MPI_Type_create_struct(2, block_lengths, offsets, types, mpi_point_type);
    MPI_Type_commit(mpi_point_type);
}

