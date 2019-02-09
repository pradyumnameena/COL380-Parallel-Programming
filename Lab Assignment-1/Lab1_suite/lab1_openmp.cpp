#include <stdlib.h>
#include "lab1_omp.h"
#include "lab1_openmp.h"
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <vector>
using namespace std;
int number_threads = 0;

vector<int> compute_centroid(vector<int> points,vector<int> centroids){
	int num_points = points.size()/3;
	int num_clusters = centroids.size()/3;
	int area = num_points/number_threads;
	vector<int> vec(num_points,0);
	int arr[num_points];

	#pragma omp parallel num_threads(number_threads)
	{	
		int tid = omp_get_thread_num();
		int start = tid*area;
		int end = start + area;
		int point_idx = 0;
		int cluster_idx = 0;
		float min_dist = 4000;
		float dist = 0;
		int min_idx = 0;
		int dist_x,dist_y,dist_z;

		for(point_idx = start;point_idx<end && point_idx<num_points;point_idx++){
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
			arr[point_idx] = min_idx;
		}

		#pragma omp barrier
	}
	for(int i = 0;i<num_points;i++){
		vec.at(i) = arr[i];
	}
	return vec;
}

vector<int> centroid_update(int num_clusters,vector<int> points,vector<int> centroid_ids){
	vector<int> vec(3*num_clusters,0);
	int number_points = points.size()/3;
	int area = number_points/number_threads;
	int num_points[num_clusters];
	int cluster[3*num_clusters];

	for(int i = 0;i<num_clusters;i++){
		num_points[i] = 0;
	}

	for(int i = 0;i<3*num_clusters;i++){
		cluster[i] = 0;
	}

	#pragma omp parallel num_threads(number_threads)
	{	
		int tid = omp_get_thread_num();
		int start = tid*area;
		int end = start + area;
		int point_idx = 0;
		int cluster_idx = 0;

		for(point_idx = start;point_idx<number_points && point_idx<end;point_idx++){
			cluster[3*centroid_ids.at(point_idx)] += points.at(3*point_idx);
			cluster[3*centroid_ids.at(point_idx) + 1] += points.at(3*point_idx + 1);
			cluster[3*centroid_ids.at(point_idx) + 2] += points.at(3*point_idx + 2);
			num_points[centroid_ids.at(point_idx)]+=1;
		}
	}

	for(int i = 0;i<num_clusters;i++){
		cluster[3*i]/=num_points[i];
		cluster[3*i + 1]/=num_points[i];
		cluster[3*i + 2]/=num_points[i];
	}

	for(int i = 0;i<3*num_clusters;i++){
		vec.at(i) = cluster[i];
	}

	return vec;
}

vector<int> get_points(int n,int* data_points){
	vector<int> vec(3*n,0);
	int idx = 0;
	while(idx<3*n){
		vec.at(idx) = *data_points;
		idx++;
		data_points++;
	}
	data_points-=(3*n);
	return vec;
}

vector<int> get_centroid_idx(int n){
	vector<int> vec(n,0);
	return vec;
}

vector<int> get_centroids(int k,int* data_points){
	vector<int> vec(3*k,0);
	int idx = 0;
	while(idx<3*k){
		vec.at(idx) = *data_points;
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

void kmeans_omp(int num_threads, int N, int K, int* data_points, int** cluster_points, int** centroid_pointer, int* num_iterations){
	vector<int> points = get_points(N,data_points);
	vector<int> cluster_ids = get_centroid_idx(N);
	vector<int> centroids = get_centroids(K,data_points);
	vector<int> all_centroids(1001*3*K);
	number_threads = num_threads;
	
	for(int i = 0;i<centroids.size();i++){
		all_centroids.push_back(centroids.at(i));
	}
	
	int count = 0;
	while(count<1000){
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