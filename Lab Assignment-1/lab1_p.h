#ifndef LAB1_SEQ_H
#define LAB1_SEQ_H

#include <vector>
using namespace std;

vector<int> get_points(int n,int* data_points);

vector<int> get_centroid_idx(int k);

vector<int> get_centroids(int k,int* data_points);

void print_vector(vector<int> vec);

#endif