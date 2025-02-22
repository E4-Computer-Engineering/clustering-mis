#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "points.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#define MAX_LINE_LENGTH 1024

int main(int argc, char **argv)
{
   int myid, nproc;
   MPI_Status status;

   char line[MAX_LINE_LENGTH];
   int row = 0;

   int nptsincluster;
   std::vector<point> pts;
   std::vector<int> clus;
   double z1, z2;

   int nmethod_proc;
   std::vector<std::string> method = {"kmeans", "dbscan"};
   int nmethod = method.size();

   int (*functions[])(int, const char *, point *, int) = {kmeans, dbscan};

   int res;

   MPI_Datatype MPI_POINT;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   create_mpi_point_type(&MPI_POINT);

   if (nmethod % nproc != 0)
   {
      printf("Error");
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   nmethod_proc = nmethod / nproc;

   if (myid == 0)
   {
      std::ifstream file("cluster_points_article.csv");
      if (!file)
      {
         std::cerr << "Error opening file.\n";
         MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }

      std::string line;
      std::getline(file, line); // Skip first line

      while (std::getline(file, line))
      {
         std::istringstream ss(line);
         point p;
         char comma;

         if (ss >> p.x >> comma >> p.y)
         { // Read x and y, assuming CSV format
            pts.emplace_back(p);
         }
      }
      nptsincluster = pts.size();
   }

   MPI_Bcast(&nptsincluster, 1, MPI_INT, 0, MPI_COMM_WORLD);
   std::cout << nptsincluster << std::endl;
   //   MPI_Bcast(pts, nptsincluster * sizeof(point), MPI_BYTE, 0, MPI_COMM_WORLD);
   if (myid != 0)
   {
      pts.resize(nptsincluster);
   }

   MPI_Bcast(pts.data(), nptsincluster, MPI_POINT, 0, MPI_COMM_WORLD);
   // Print received data in each process
   /*   for (int i = 0; i < nptsincluster; i++)
      {
        printf("I am proc %d and in pos %d, I received %f %f\n", myid, i, pts[i].x, pts[i].y);
      }
      printf("\n");
   */

   for (int im = 0; im < nmethod_proc; im++)
   {
      std::string mymethod = method[nmethod_proc * myid + im];
      printf("I am proc %d and I will deal with method %s\n", myid, mymethod.c_str());
      for (int i = 0; i < 2; i++)
      {
         if (mymethod == method[i])
         {
            res = functions[i](myid, mymethod.c_str(), pts.data(), nptsincluster);
         }
      }
   }

   MPI_Finalize();
   return 0;
}
