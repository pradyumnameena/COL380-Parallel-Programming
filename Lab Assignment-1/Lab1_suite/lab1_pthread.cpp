#include <stdlib.h>
#include "lab1_pthread.h"
#include "lab1_p.h"
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <vector>
using namespace std;

vector<int> get_points(int n,int* data_points){
	vector<int> vec;
	int idx = 0;
	while(idx<3*n){
		vec.push_back(*data_points);
		idx++;
		data_points++;
	}
	data_points-=(3*n);
	return vec;
}

vector<int> get_centroid_idx(int n){
	vector<int> vec;
	for(int i = 0;i<n;i++){
		vec.push_back(0);
	}
	return vec;
}

vector<int> get_centroids(int k,int* data_points){
	vector<int> vec;
	int idx = 0;
	while(idx<3*k){
		vec.push_back(*data_points);
		idx++;
		data_points++;
	}
	data_points-=(3*k);
	return vec;
}

void print_vector(vector<int> vec){
	int i = 0;
	int size = vec.size();
	for(i = 0;i<size;i++){
		printf("%d\n",vec.at(i));
	}
	return;
}

void kmeans_pthread(int num_threads, int N, int K, int* data_points, int** data_point_cluster, int** centroid_pointer, int* num_iterations){
	vector<int> points = get_points(N,data_points);
	vector<int> cluster_ids = get_centroid_idx(N);
	vector<int> centroids = get_centroids(K,data_points);
	vector<int> all_centroids;
	
	// for(int i = 0;i<centroids.size();i++){
	// 	all_centroids.push_back(centroids.at(i));
	// }
	
	// int count = 0;
	// while(count<1000){
	// 	cluster_ids = compute_centroid(points,centroids);
	// 	centroids = centroid_update(K,points,cluster_ids);
	// 	count++;
	// 	for(int i = 0;i<centroids.size();i++){
	// 		all_centroids.push_back(centroids.at(i));
	// 	}
	// }
	// *num_iterations = count;
	
	// *cluster_points = (int*)malloc(sizeof(int)*(N*4));
	// int* cluster_pointer = *cluster_points;
	// for(int i = 0;i<N;i++){
	// 	*cluster_pointer = points.at(3*i);
	// 	cluster_pointer++;
	// 	*cluster_pointer = points.at(3*i + 1);
	// 	cluster_pointer++;
	// 	*cluster_pointer = points.at(3*i + 2);
	// 	cluster_pointer++;
	// 	*cluster_pointer = cluster_ids.at(i);
	// 	cluster_pointer++;
	// }

	// *centroid_pointer = (int*)malloc(sizeof(int)*3*K*(1+count));
	// int* centroid_point = *centroid_pointer;
	// for(int i = 0;i<all_centroids.size();i++){
	// 	*centroid_point = all_centroids.at(i);
	// 	centroid_point++;
	// }
	return;
}