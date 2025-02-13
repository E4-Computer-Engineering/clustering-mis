#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "points.h"

#define MAX_ROWS 5
#define MAX_LINE_LENGTH 1024

int main(int argc, char **argv)
{
   int myid, nproc;
   MPI_Status status;

   char line[MAX_LINE_LENGTH];
   int row = 0;

   int nptsincluster;

   point *pts;
   int *clus = malloc(nptsincluster * sizeof(int));

   double z1, z2;

   MPI_Datatype MPI_POINT;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   create_mpi_point_type(&MPI_POINT);

   nptsincluster = MAX_ROWS -1;

   pts = calloc(nptsincluster, sizeof(point));

   if (myid == 0)
   {
        /* Read file */
        FILE *file_in = fopen("cluster_points_test.csv", "r");
        if (!file_in) {
           perror("Unable to open file!");
           return EXIT_FAILURE;
        }

        // Read header line and ignore it
        fgets(line, sizeof(line), file_in);

        // Read data line by line

        while (fgets(line, sizeof(line), file_in) && row < nptsincluster) {

            char *token = strtok(line, ",");
            if (token) z1 = atof(token);

            token = strtok(NULL, ",");
            if (token) z2 = atof(token);

            token = strtok(NULL, ",");
            if (token) clus[row] = atoi(token);

            printf("Leggo %f, %f\n", z1, z2);

	    pts[row].x = z1;
            pts[row].y = z2;
             
            /* Pointer to raw data */
            //config.objs[row] = &(pts[row]);

            row++;
        }

        fclose(file_in);

        /*End of read file*/
   	   
   }

//   MPI_Bcast(pts, nptsincluster * sizeof(point), MPI_BYTE, 0, MPI_COMM_WORLD);
   
   MPI_Bcast(pts, nptsincluster, MPI_POINT, 0, MPI_COMM_WORLD);
   // Print received data in each process
   for (int i = 0; i < nptsincluster; i++)
   {  
     printf("I am proc %d and in pos %d, I received %f %f\n", myid, i, pts[i].x, pts[i].y);
   }
   printf("\n");

   free(pts);
   MPI_Finalize();
   return 0;
}
