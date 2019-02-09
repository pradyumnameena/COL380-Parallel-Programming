#include <stdlib.h>
#include "lab1_openmp.h"
#include "lab1_omp.h"
#include <omp.h>
#include <stdio.h>
using namespace std;
int thread_count = 0;

void compute_centroid(int* data_points,int* centroids,int* cluster_ids,int n,int k){
	int area = n/thread_count;
	
	#pragma omp parallel num_threads(thread_count)
	{	
		int point_idx = 0;
		int cluster_idx = 0;
		int dist = 0;
		int min_dist = 0;
		int min_idx = 0;
		int dist_x,dist_y,dist_z;
		int tid = omp_get_thread_num();
		int start = tid*area;
		int end = start+area;
		
		for(point_idx = start;point_idx<end && point_idx<n;point_idx++){
			min_dist = 6000000;
			min_idx = 0;
			for(cluster_idx = 0;cluster_idx<k;cluster_idx++){
				dist_x = *(data_points+3*point_idx) - *(centroids+3*cluster_idx);
				dist_y = *(data_points+3*point_idx + 1) - *(centroids+3*cluster_idx + 1);
				dist_z = *(data_points+3*point_idx + 2) - *(centroids+3*cluster_idx + 2);
				dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
				if(dist<min_dist){
					min_dist = dist;
					min_idx = cluster_idx;
				}
			}
			*(cluster_ids+point_idx) = min_idx;
		}
	}
	return;
}

void centroid_update(int* data_points,int* centroids,int* cluster_ids,int* count_points,int n,int k){
	// initialize(centroids,3*k,0);
	int i = 0;
	while(i<3*k){
		*centroids = 0;
		i++;
		centroids++;
	}
	centroids-=3*k;
	// initialize(count_points,3*k,0);
	i = 0;
	while(i<k){
		*count_points = 0;
		i++;
		count_points++;
	}
	count_points-=k;
	int area = n/thread_count;

	#pragma omp parallel num_threads(thread_count)
	{
		int point_idx = 0;
		int cluster_idx = 0;
		int idx = 0;
		int tid = omp_get_thread_num();
		int start = tid*area;
		int end = start+area;

		for(point_idx = start;point_idx<n && point_idx<end;point_idx++){
			idx = *(cluster_ids+point_idx);
			*(centroids + 3*idx)+=*(data_points+3*point_idx);
			*(centroids + 3*idx + 1)+=*(data_points+3*point_idx + 1);
			*(centroids + 3*idx + 2)+=*(data_points+3*point_idx + 2);
			*(count_points+idx)+=1;
		}
	}
	int cluster_idx = 0;

	for(cluster_idx = 0;cluster_idx<k;cluster_idx++){
		if(*(count_points+cluster_idx)!=0){
			*(centroids + 3*cluster_idx)/=*(count_points+cluster_idx);
			*(centroids + 3*cluster_idx + 1)/=*(count_points+cluster_idx);
			*(centroids + 3*cluster_idx + 2)/=*(count_points+cluster_idx);
		}
	}
	return;
}

void initialize(int* pointer,int n,int val){
	int i = 0;
	while(i<n){
		*pointer = val;
		i++;
		pointer++;
	}
	pointer-=n;
	return;
}

void centroid_initialize(int* centroid,int* data_points,int num_cluster){
	int i = 0;
	while(i<3*num_cluster){
		*centroid = *data_points;
		centroid++;
		data_points++;
		i++;
	}
	centroid-=3*num_cluster;
	data_points-=3*num_cluster;
	return;
}

void kmeans_omp(int num_threads,int N,int K,int* data_points,int** cluster_points,int** centroid_pointer,int* num_iterations){
	double start_time = omp_get_wtime();
	thread_count = num_threads;
	int* centroid_ids;
	int* centroid;
	int* count_points;
	int centroid_idx = 0;
	int count = 0;
	int max_iterations = 200;
	vector<int> all_centroids((max_iterations+1)*K*3,0);

	centroid_ids = (int*)malloc(sizeof(int)*N);
	centroid = (int*)malloc(sizeof(int)*3*K);
	count_points = (int*)malloc(sizeof(int)*K);
	initialize(centroid_ids,N,0);
	centroid_initialize(centroid,data_points,K);

	// adding the initial centroid's coordinates
	while(centroid_idx<3*K){
		all_centroids.at(centroid_idx) = *centroid;
		centroid_idx++;
		centroid++;
	}
	centroid-=3*K;
	
	// main algorithm
	while(count<max_iterations){
		compute_centroid(data_points,centroid,centroid_ids,N,K);
		centroid_update(data_points,centroid,centroid_ids,count_points,N,K);
		
		for(int i = 0;i<3*K;i++){
			all_centroids.at(centroid_idx) = *centroid;
			centroid++;
			centroid_idx++;
		}
		centroid-=3*K;
		count++;
	}
	*num_iterations = count;
	
	// adding the points along with their centroid index in the pointer's location
	*cluster_points = (int*)malloc(sizeof(int)*(N*4));
	int* cluster_pointer = *cluster_points;
	int iter_counter = 0;
	while(iter_counter<N){
		*cluster_pointer = *data_points;
		cluster_pointer++;
		data_points++;
		*cluster_pointer = *data_points;
		cluster_pointer++;
		data_points++;
		*cluster_pointer = *data_points;
		cluster_pointer++;
		data_points++;
		*cluster_pointer = *centroid_ids;
		cluster_pointer++;
		centroid_ids++;
		iter_counter++;
	}

	// adding all centroids of various iterations into the pointer's location
	*centroid_pointer = (int*)malloc(sizeof(int)*3*K*(1+count));
	int* centroid_point = *centroid_pointer;
	for(int i = 0;i<all_centroids.size();i++){
		*centroid_point = all_centroids.at(i);
		centroid_point++;
	}
	double end_time = omp_get_wtime();
	double computation_time = (end_time - start_time);
	printf("Time Taken by openmp algorithm: %lf \n", computation_time);
	return;
}