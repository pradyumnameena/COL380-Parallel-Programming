#include <stdlib.h>
#include "lab1_sequential.h"
#include "lab1_seq.h"
#include <stdio.h>
#include <math.h>
#include <vector>
using namespace std;

vector<int> compute_centroid(vector<int> points,vector<int> centroids){
	vector<int> vec;
	int num_points = points.size()/3;
	int num_clusters = centroids.size()/3;
	int point_idx = 0;
	int cluster_idx = 0;
	float min_dist = 4000;
	float dist = 0;
	int min_idx = 0;
	int dist_x,dist_y,dist_z;

	for(point_idx = 0;point_idx<num_points;point_idx++){
		min_dist = 4000;
		for(cluster_idx = 0;cluster_idx<num_clusters;cluster_idx++){
			dist_x = (points.at(3*point_idx) - centroids.at(3*cluster_idx))*(points.at(3*point_idx) - centroids.at(3*cluster_idx));
			dist_y = (points.at(3*point_idx + 1) - centroids.at(3*cluster_idx + 1))*(points.at(3*point_idx + 1) - centroids.at(3*cluster_idx + 1));
			dist_z = (points.at(3*point_idx + 2) - centroids.at(3*cluster_idx + 2))*(points.at(3*point_idx + 2) - centroids.at(3*cluster_idx + 2));
			dist = sqrt(dist_x+dist_y+dist_z);
			if(dist<min_dist){
				min_dist = dist;
				min_idx = cluster_idx;
			}
		}
		vec.push_back(min_idx);
	}
	return vec;
}

vector<int> centroid_update(int num_clusters,vector<int> points,vector<int> centroid_ids){
	vector<int> vec;
	vector<int> num_points;
	int number_points = points.size()/3;
	int point_idx = 0;
	int cluster_idx = 0;
	num_points.resize(num_clusters);
	vec.resize(3*num_clusters);
	
	for(point_idx = 0;point_idx<number_points;point_idx++){
		vec.at(3*centroid_ids.at(point_idx))+=points.at(3*point_idx);
		vec.at(3*centroid_ids.at(point_idx) + 1)+=points.at(3*point_idx + 1);
		vec.at(3*centroid_ids.at(point_idx) + 2)+=points.at(3*point_idx + 2);
		num_points.at(centroid_ids.at(point_idx))+=1;
	}

	for(int i = 0;i<num_clusters;i++){
		vec.at(3*i)/=num_points.at(i);
		vec.at(3*i + 1)/=num_points.at(i);
		vec.at(3*i + 2)/=num_points.at(i);
	}
	return vec;
}

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

void kmeans_sequential(int N,int K,int* data_points,int** cluster_points,int** centroid_pointer,int* num_iterations){
	vector<int> points = get_points(N,data_points);
	vector<int> cluster_ids = get_centroid_idx(N);
	vector<int> centroids = get_centroids(K,data_points);
	vector<int> all_centroids;
	for(int i = 0;i<centroids.size();i++){
		all_centroids.push_back(centroids.at(i));
	}
	
	int count = 0;
	while(count<10){
		cluster_ids = compute_centroid(points,centroids);
		centroids = centroid_update(K,points,cluster_ids);
		count++;
		for(int i = 0;i<centroids.size();i++){
			all_centroids.push_back(centroids.at(i));
		}
	}
	*num_iterations = count;
	
	*cluster_points = (int*)malloc(sizeof(int)*(N*4));
	int* cluster_pointer = *cluster_points;
	for(int i = 0;i<N;i++){
		*cluster_pointer = points.at(3*i);
		cluster_pointer++;
		*cluster_pointer = points.at(3*i + 1);
		cluster_pointer++;
		*cluster_pointer = points.at(3*i + 2);
		cluster_pointer++;
		*cluster_pointer = cluster_ids.at(i);
		cluster_pointer++;
	}

	*centroid_pointer = (int*)malloc(sizeof(int)*3*K*(1+count));
	int* centroid_point = *centroid_pointer;
	for(int i = 0;i<all_centroids.size();i++){
		*centroid_point = all_centroids.at(i);
		centroid_point++;
	}
	return;
}