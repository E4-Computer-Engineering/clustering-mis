#include "kmeans_cl.h"
#include "common.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>

static double pt_distance(const Pointer a, const Pointer b) {
    point *pa = (point *)a;
    point *pb = (point *)b;

    double dx = (pa->x - pb->x);
    double dy = (pa->y - pb->y);

    return dx * dx + dy * dy;
}

static void pt_centroid(const Pointer *objs, const int *clusters,
                        size_t num_objs, int cluster, Pointer centroid) {
    int i;
    int num_cluster = 0;
    point sum;
    point **pts = (point **)objs;
    point *center = (point *)centroid;

    sum.x = sum.y = 0.0;

    if (num_objs <= 0)
        return;

    for (i = 0; i < num_objs; i++) {
        /* Only process objects of interest */
        if (clusters[i] != cluster)
            continue;

        sum.x += pts[i]->x;
        sum.y += pts[i]->y;
        num_cluster++;
    }
    if (num_cluster) {
        sum.x /= num_cluster;
        sum.y /= num_cluster;
        *center = sum;
    }
    return;
}

int kmeans(int myrank, const char *str, const point *pts, int n, int *res) {

    printf("I am rank %d and I am in the function %s\n", myrank, str);

    kmeans_config config;
    kmeans_result result;
    int i;
    point *init;
    bool print_results = false;
    unsigned long start;

    int k = 7;

    srand(123);

    /* Constants */
    config.k = k;
    config.num_objs = n;
    config.max_iterations = 200;
    config.distance_method = pt_distance;
    config.centroid_method = pt_centroid;

    /* Inputs for K-means */
    config.objs = calloc(config.num_objs, sizeof(Pointer));
    config.centers = calloc(config.k, sizeof(Pointer));
    config.clusters = calloc(config.num_objs, sizeof(int));

    /* Storage for raw data */

    init = calloc(config.k, sizeof(point));

    for (i = 0; i < n; i++) {
        /* Pointer to raw data */
        config.objs[i] = &(pts[i]);
    }

    /* Populate the initial means vector with random start points */
    for (i = 0; i < config.k; i++) {
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
    result = kmeans_cl(&config);

    /* print results */

    if (print_results) {
        FILE *file_out = fopen("kmeans_output_data.txt", "w");
        if (!file_out) {
            perror("Unable to open file!");
            return EXIT_FAILURE;
        }
        // Write header
        fprintf(file_out, "Feature 1\tFeature 2\tCluster\n");

        for (i = 0; i < config.num_objs; i++) {
            point *pt = (point *)(config.objs[i]);

            if (config.objs[i]) {
                // printf("here");
                fprintf(file_out, "%g\t%g\t%d\n", pt->x, pt->y,
                        config.clusters[i]);
            } else {
                // printf("there");
                fprintf(file_out, "N\tN\t%d\n", config.clusters[i]);
            }
        }
        fclose(file_out);
    }

    for (i = 0; i < n; i++) {
        res[i] = config.clusters[i];
    }
    free(config.objs);
    free(config.clusters);
    free(config.centers);

    free(init);

    return 0;
}
