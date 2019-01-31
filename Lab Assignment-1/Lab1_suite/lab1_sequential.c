#include <stdlib.h>
#include "lab1_sequential.h"
#include <stdio.h>
#include <math.h>

float distance(int* a,int* b){
	float dist = 0;
	int i = 3;
	while(i>0){
		dist+=(*a - *b)*(*a - *b);
		a++;
		b++;
		i-=1;
	}
	return sqrt(dist);
}

void assign_centroid(int n,int k,int* data_points,int* cluster_points,int* centroids,int* pointer){
	int point_idx = 0;
	int cluster_idx = 0;
	float dist = 0;
	float min_dist = 0;
	int idx = 0;
	for(point_idx = 0;point_idx<n;point_idx++){
		min_dist = 4000;
		for(cluster_idx = 0;cluster_idx<k;cluster_idx++){
			dist = distance(data_points+3*point_idx,centroids+3*cluster_idx);
			if(dist<min_dist){
				min_dist = dist;
	 			idx = cluster_idx;
			}
		}
		*(pointer + idx)+=1;
		*(cluster_points+4*point_idx+3) = idx;
	}
	return;
}

void compute_centroid(int n,int k,int* data_points,int* cluster_points,int* centroids,int* pointer){
	int point_idx = 0;
	int cluster_idx = 0;
	int idx = 0;

	for(point_idx = 0;point_idx<n;point_idx++){
		cluster_idx = *(cluster_points+4*point_idx+3);
		*(centroids + cluster_idx)+=*(data_points + 3*point_idx);
		*(centroids + cluster_idx + 1)+=*(data_points + 3*point_idx + 1);
		*(centroids + cluster_idx + 2)+=*(data_points + 3*point_idx + 2);
	}

	for(idx = 0;idx<k;idx++){
		*(centroids+3*idx) = (int)((*(centroids+3*idx))/(*(pointer+idx)));
		*(centroids+3*idx+1) = (int)((*(centroids+3*idx+1))/(*(pointer+idx)));
		*(centroids+3*idx+2) = (int)((*(centroids+3*idx+2))/(*(pointer+idx)));
	}
	return;
}

void initialize(int n,int k,int* data_points,int* cluster_points,int* centroids){
	int cluster_idx = 0;
	int point_idx = 0;
	int count = 3;
	
	while(cluster_idx<3*k){
		*centroids = *data_points;
		centroids++;
		data_points++;
		cluster_idx++;
	}
	centroids-=(3*k);
	data_points-=(3*k);

	while(point_idx<n){
		// initializing cluster_points array
		count = 3;
		while(count!=0){
			*cluster_points = *data_points;
			cluster_points++;
			data_points++;
			count-=1;
		}
		*cluster_points = 0;
		cluster_points++;
		point_idx++;
	}
	data_points-=(3*n);
	cluster_points-=(4*n);
	return;
}

void settozero(int* pointer, int k){
	int i = k;
	while(i!=0){
		*pointer = 0;
		pointer++;
		i-=1;
	}
}

void kmeans_sequential(int N,int K,int* data_points,int** cluster_points,int** centroids,int* num_iterations){
	*cluster_points = (int*)malloc(sizeof(int)*(N*4));
	*centroids = (int*)malloc(sizeof(int)*(K*3));
	int* cluster_pointer = *cluster_points;
	int* centroid_pointer = *centroids;
	int* pointer = (int*)malloc(sizeof(int)*K);
	int count = 2;

	initialize(N,K,data_points,cluster_pointer,centroid_pointer);
	while(count<10){
		settozero(pointer,K);
		assign_centroid(N,K,data_points,cluster_pointer,centroid_pointer,pointer);
		centroid_pointer = realloc(centroid_pointer,sizeof(int)*(K*3)*count);
		compute_centroid(N,K,data_points,cluster_pointer,centroid_pointer,pointer);
		count++;
	}
	*num_iterations = count;
	return;
}