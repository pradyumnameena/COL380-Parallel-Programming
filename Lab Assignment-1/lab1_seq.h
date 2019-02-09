#ifndef LAB1_SEQ_H
#define LAB1_SEQ_H

#include <vector>
using namespace std;

void compute_centroid(int* data_points,int* centroids,int* cluster_ids,int n,int k);

void centroid_update(int* data_points,int* centroids,int* cluster_ids,int* count_points,int n,int k);

void initialize(int* pointer,int n,int val);

void centroid_initialize(int* centroid,int* data_points,int num_cluster);

#endif

// many ops affectes time like taking sqrt took 24 or 14 seconds difference.
// instead of making fuction for distance compute at the same place reduced 7 seconds