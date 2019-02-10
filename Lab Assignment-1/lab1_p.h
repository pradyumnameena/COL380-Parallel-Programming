#ifndef LAB1_P_H
#define LAB1_P_H

#include <vector>
using namespace std;

void *compute_centroid_thread(void* thread_id);

void compute_centroid();

void *centroid_update_thread(void *thread_id);

void *centroid_update_thread_V2(void *thread_id);

void *centroid_update_thread_V3(void *thread_id);

void centroid_update();

void initialize(int* pointer,int n,int val);

void centroid_initialize(float* centroid,int* data_points,int num_cluster);

#endif