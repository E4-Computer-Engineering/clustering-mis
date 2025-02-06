#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

#include "kmeans.h"

#define MAX_ROWS 789  // Adjust based on expected number of rows
#define MAX_LINE_LENGTH 1024

typedef struct point
{
	double x;
	double y;
} point;


static double pt_distance(const Pointer a, const Pointer b)
{
	point *pa = (point*)a;
	point *pb = (point*)b;

	double dx = (pa->x - pb->x);
	double dy = (pa->y - pb->y);

	return dx*dx + dy*dy;
}

static void pt_centroid(const Pointer * objs, const int * clusters, size_t num_objs, int cluster, Pointer centroid)
{
	int i;
	int num_cluster = 0;
	point sum;
	point **pts = (point**)objs;
	point *center = (point*)centroid;

	sum.x = sum.y = 0.0;

	if (num_objs <= 0) return;

	for (i = 0; i < num_objs; i++)
	{
		/* Only process objects of interest */
		if (clusters[i] != cluster) continue;

		sum.x += pts[i]->x;
		sum.y += pts[i]->y;
		num_cluster++;
	}
	if (num_cluster)
	{
		sum.x /= num_cluster;
		sum.y /= num_cluster;
		*center = sum;
	}
	return;
}

int main(int nargs, char **args)
{
	kmeans_config config;
	kmeans_result result;
	int i;
	point *pts;
	point *init;
	int print_results = 1;
	unsigned long start;

	int nptsincluster;
	int k = 7;

        int row = 0;
        char line[MAX_LINE_LENGTH];

	double z1, z2;

	srand(123);

	nptsincluster = MAX_ROWS -1;
        int *clus = malloc(nptsincluster * sizeof(int));

	/* Constants */
	config.k = k;
	config.num_objs = nptsincluster;
	config.max_iterations = 200;
	config.distance_method = pt_distance;
	config.centroid_method = pt_centroid;

	/* Inputs for K-means */
	config.objs = calloc(config.num_objs, sizeof(Pointer));
	config.centers = calloc(config.k, sizeof(Pointer));
	config.clusters = calloc(config.num_objs, sizeof(int));

	/* Storage for raw data */
	pts = calloc(config.num_objs, sizeof(point));
	init = calloc(config.k, sizeof(point));

        /* Read file */
        FILE *file_in = fopen("cluster_points_article.csv", "r");
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

	    pts[row].x = z1;
            pts[row].y = z2;

            /* Pointer to raw data */
            config.objs[row] = &(pts[row]);

            row++;
        }

        fclose(file_in);

	/*End of read file*/

	/* Populate the initial means vector with random start points */
	for (i = 0; i < config.k; i++)
	{
		int r = lround(config.num_objs * (1.0 * rand() / RAND_MAX));
		/* Populate raw data */
		init[i] = pts[r];
		/* Pointers to raw data */
		config.centers[i] = &(init[i]);

		if (print_results)
			printf("center[%d]\t%g\t%g\n", i, init[i].x, init[i].y);
	}

	/* run k-means! */
	start = time(NULL);
	result = kmeans(&config);

	printf("\n");
	printf("Iteration count: %d\n", config.total_iterations);
	printf("     Time taken: %ld seconds\n", (time(NULL) - start));
	printf(" Iterations/sec: %.3g\n", (1.0*config.total_iterations)/(time(NULL) - start));
	printf("\n");

	/* print results */
	
	if (print_results)
	{
		
	        FILE *file_out = fopen("output_data.txt", "w");
                if (!file_out) {
                   perror("Unable to open file!");
                   return EXIT_FAILURE;
                }

                // Write header
                fprintf(file_out, "Feature 1\tFeature 2\tCluster\n");
		
		for (i = 0; i < config.num_objs; i++)
		{
			point *pt = (point*)(config.objs[i]);

			if (config.objs[i])
				fprintf(file_out, "%g\t%g\t%d\n", pt->x, pt->y, config.clusters[i]);
			else
				fprintf(file_out, "N\tN\t%d\n", config.clusters[i]);
		}
	        fclose(file_out);
	}
         
	free(config.objs);
	free(config.clusters);
	free(config.centers);

	free(init);
	free(pts);
	free(clus);

        return 0;
}

