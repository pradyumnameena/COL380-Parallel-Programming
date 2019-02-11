#include <stdlib.h>
#include "lab1_pthread.h"
#include "lab1_p.h"
#include <stdio.h>
#include <omp.h>
#include <pthread.h>
using namespace std;

int thread_count = 0;
int* data_points_pointer;
float* centroids_pointer;
int* cluster_ids_pointer;
int* num_points_pointer;
float* helper_pointer;
int num_points;
int num_clusters;
int area = 0;

void *compute_centroid_thread(void *thread_id){
	long t;
	t = (long)thread_id;
	int tid = (int)t;
	int point_idx = 0;
	int cluster_idx = 0;
	int start = tid*area;
	int end = start+area;
	int min_idx = 0;
	float dist = 0;
	float min_dist = 0;
	float dist_x,dist_y,dist_z;
	
	for(point_idx = start;point_idx<end && point_idx<num_points;point_idx++){
		min_dist = 60000000;
		min_idx = 0;
		for(cluster_idx = 0;cluster_idx<num_clusters;cluster_idx++){
			dist_x = *(data_points_pointer+3*point_idx) - *(centroids_pointer+3*cluster_idx);
			dist_y = *(data_points_pointer+3*point_idx + 1) - *(centroids_pointer+3*cluster_idx + 1);
			dist_z = *(data_points_pointer+3*point_idx + 2) - *(centroids_pointer+3*cluster_idx + 2);
			dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
			if(dist<min_dist){
				min_dist = dist;
				min_idx = cluster_idx;
			}
		}
		*(cluster_ids_pointer+point_idx) = min_idx;
	}
	pthread_exit(0);
}

void compute_centroid(){
	pthread_t thread_arr[thread_count];
	int rc = 0;
	for(int i = 0;i<thread_count;i++){
		rc = pthread_create(&thread_arr[i], NULL, compute_centroid_thread, (void *)i);
	}

	for(int i = 0;i<thread_count;i++){
		pthread_join(thread_arr[i], NULL);
	}
	return;
}

void *centroid_update_thread(void *thread_id){
	long t;
	t = (long)thread_id;
	int tid = (int)t;
	int point_idx = 0;
	int cluster_idx = 0;
	int start = tid*area;
	int end = start+area;
	int idx = 0;
	
	for(point_idx = start;point_idx<end && point_idx<num_points;point_idx++){
		idx = *(cluster_ids_pointer+point_idx);
		*(helper_pointer + num_clusters*3*tid + 3*idx) += *(data_points_pointer + 3*point_idx);
		*(helper_pointer + num_clusters*3*tid + 3*idx + 1) += *(data_points_pointer + 3*point_idx + 1);
		*(helper_pointer + num_clusters*3*tid + 3*idx + 2) += *(data_points_pointer + 3*point_idx + 2);
		*(num_points_pointer + idx)+=1;
	}
	pthread_exit(0);
}

void *centroid_update_thread_V2(void *thread_id){
	long t;
	t = (long)thread_id;
	int tid = (int)t;
	int counter = 0;
	for(counter = 0;counter<thread_count;counter++){
		*(centroids_pointer + 3*tid)+=*(helper_pointer + 3*num_clusters*counter + 3*tid);
		*(centroids_pointer + 3*tid + 1)+=*(helper_pointer + 3*num_clusters*counter + 3*tid + 1);
		*(centroids_pointer + 3*tid + 2)+=*(helper_pointer + 3*num_clusters*counter + 3*tid + 2);
	}
	
	if(*(num_points_pointer+tid)!=0){
		*(centroids_pointer + 3*tid)/=*(num_points_pointer+tid);
		*(centroids_pointer + 3*tid + 1)/=*(num_points_pointer+tid);
		*(centroids_pointer + 3*tid + 2)/=*(num_points_pointer+tid);
	}
	pthread_exit(0);
}

void centroid_update(){
	int i = 0;
	while(i<3*num_clusters){
		*centroids_pointer = 0;
		i++;
		centroids_pointer++;
	}
	centroids_pointer-=3*num_clusters;
	
	i = 0;
	while(i<num_clusters){
		*num_points_pointer = 0;
		i++;
		num_points_pointer++;
	}
	num_points_pointer-=num_clusters;

	i = 0;
	while(i<3*num_clusters*thread_count){
		*helper_pointer = 0;
		i++;
		helper_pointer++;
	}
	helper_pointer-=(3*num_clusters*thread_count);

	pthread_t thread_arr[thread_count];
	int rc = 0;
	for(int i = 0;i<thread_count;i++){
		rc = pthread_create(&thread_arr[i], NULL, centroid_update_thread, (void *)i);
	}

	for(int i = 0;i<thread_count;i++){
		pthread_join(thread_arr[i], NULL);
	}
	
	pthread_t thread_arr2[num_clusters];
	rc = 0;
	for(int i = 0;i<num_clusters;i++){
		rc = pthread_create(&thread_arr2[i], NULL, centroid_update_thread_V2, (void *)i);
	}

	for(int i = 0;i<num_clusters;i++){
		pthread_join(thread_arr2[i], NULL);
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

void centroid_initialize(float* centroid,int* data_points,int num_cluster){
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

void kmeans_pthread(int num_threads,int N,int K,int* data_points,int** cluster_points,float** centroid_pointer,int* num_iterations){
	double start_time = omp_get_wtime();
	
	// declaration of various datatypes to be used
	int* centroid_ids;
	float* centroid;
	int* count_points;
	int centroid_idx = 0;
	int count = 0;
	int max_iterations = 300;
	vector<float> all_centroids((max_iterations+1)*K*3,0);
	
	// malloc and initializing the centroids and cluster_ids array
	helper_pointer = (float*)malloc(sizeof(float)*3*K*(thread_count));
	centroid_ids = (int*)malloc(sizeof(int)*N);
	centroid = (float*)malloc(sizeof(float)*3*K);
	count_points = (int*)malloc(sizeof(int)*K);
	initialize(centroid_ids,N,0);
	centroid_initialize(centroid,data_points,K);
	
	// setting up the global variables
	data_points_pointer = data_points;
	centroids_pointer = centroid;
	cluster_ids_pointer = centroid_ids;
	num_points_pointer = count_points;
	area = N/num_threads;
	num_points = N;
	num_clusters = K;
	thread_count = num_threads;

	// adding the initial centroid's coordinates
	while(centroid_idx<3*K){
		all_centroids.at(centroid_idx) = *centroid;
		centroid_idx++;
		centroid++;
	}
	centroid-=3*K;
	
	// main algorithm
	while(count<max_iterations){
		compute_centroid();
		centroid_update();
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
	*centroid_pointer = (float*)malloc(sizeof(float)*3*K*(1+count));
	float* centroid_point = *centroid_pointer;
	for(int i = 0;i<all_centroids.size();i++){
		*centroid_point = all_centroids.at(i);
		centroid_point++;
	}
	double end_time = omp_get_wtime();
	double computation_time = (end_time - start_time);
	printf("Time Taken by pthread algorithm: %lf \n", computation_time);
	return;
}