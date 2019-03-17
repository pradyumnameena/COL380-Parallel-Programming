#include <iostream>
#include <cstdlib>
#include <omp.h>
#include "lab2_omp.h"
using namespace std;

float* matrix_mult(int m1,int n1,int m2,int n2,float* mat1,float* mat2);

float* matrix_transpose(int m,int n,float* mat);

float* matrix_mult(int m1,int n1,int m2,int n2,float* mat1,float* mat2){
	float* res_matrix = (float*)malloc(sizeof(float)*m1*n2);
	for(int i = 0;i<m1*n2;i++){
		*(res_matrix+i) = 0;
	}

	for(int i = 0;i<m1;i++){
		for(int j = 0;j<n2;j++){
			for(int k = 0;k<n1;k++){
				*(res_matrix+i*n2+j)+= (*(mat1+i*n1+k)) * (*(mat2+k*n2+j));
			}
		}
	}
	return res_matrix;
}

float* matrix_transpose(int m,int n,float* mat){
	float* res_matrix = (float*)malloc(sizeof(float)*m*n);
	for(int i = 0;i<m;i++){
		for(int j = 0;j<n;j++){
			*(res_matrix+j*m+i) = *(mat+i*n+j);
		}
	}
	return res_matrix;
}

void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{	
	float* D_transpose = matrix_transpose(M,N,D);
	float* matrix_prod = matrix_mult(M,N,N,M,D,D_transpose);
	free(D_transpose);
	
	// computing the eigen values and sorting them in descending order
		while(true){
			for(i = 0;i<M;i++){
				for(j = 0;j<M;j++){

				}
			}
		}
}

void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    
}
