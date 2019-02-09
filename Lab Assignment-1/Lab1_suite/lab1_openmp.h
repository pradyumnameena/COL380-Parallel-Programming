#ifndef LAB1_OPENMP_H
#define LAB1_OPENMP_H

#include <vector>
using namespace std;

vector<int> compute_centroid(vector<int> points,vector<int> centroids);

vector<int> centroid_update(int k,vector<int> points,vector<int> centroid_ids);

vector<int> get_points(int n,int* data_points);

vector<int> get_centroid_idx(int k);

vector<int> get_centroids(int k,int* data_points);

void print_vector(vector<int> vec);

#endif