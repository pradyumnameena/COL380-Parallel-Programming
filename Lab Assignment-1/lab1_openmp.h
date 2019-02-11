#ifndef LAB1_OPENMP_H
#define LAB1_OPENMP_H

#include <vector>
using namespace std;

void compute_centroid(int* data_points,float* centroids,int* cluster_ids,int n,int k);

void centroid_update(int* data_points,float* centroids,int* cluster_ids,int* count_points,int n,int k);

void initialize(int* pointer,int n,int val);

void centroid_initialize(float* centroid,int* data_points,int num_cluster);

#endif