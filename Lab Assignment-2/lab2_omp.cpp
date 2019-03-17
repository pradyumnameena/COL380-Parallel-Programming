#include <iostream>
#include <cstdlib>
#include <omp.h>
#include "lab2_omp.h"
using namespace std;

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
}

void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    // sigma has square root of eigen values hence all will be positive
    int i = 0;
    int j = 0;
    int idx = 0;
    int index = 0;
    float sum = 0;
    float temp = 0.0;
    float* reduced_variance = (float*)malloc(sizeof(float)*N);
    int* eigen_value_sorted_index = (int*)malloc(sizeof(int)*N);

    // sum of eigen values
    for(idx = 0;idx<N;idx++){
    	sum+=(*(SIGMA+idx));
    }
    
    // reduction variances
    #pragma omp parallel for
    {
    	for(idx = 0;idx<N;idx++){
	    	*(reduced_variance+idx) = (*(SIGMA+idx))/sum;
	    	*(eigen_value_sorted_index+idx) = idx;
	    }
    }

    // sorting the eigen values while keeping track of the corresponding column
    for(i = 0;i<N-1;i++){
    	for(j = i+1;j<N;j++){
    		if((*(reduced_variance+j)) > (*(reduced_variance+i))){
    			temp = *(reduced_variance+j);
    			*(reduced_variance+j) = *(reduced_variance+i);
    			*(reduced_variance+i) = temp;
    			*(eigen_value_sorted_index+i) = j;
    			*(eigen_value_sorted_index+j) = i;
    		}
    	}
    }

    // eigen_value_sorted_index stores ki is index par kaunse index ka eigenvector ki cooresponding value hai
    idx = 0;
    float summed_variance = 0;
    float retention_reqd = retention*0.01;
    for(idx = 0;idx<N;idx++){
    	summed_variance+=(*(reduced_variance+idx));
    	if(summed_variance>=retention_reqd){
    		*K = idx+1;
    		break;
    	}
    }
    // check ki agar 100 daalen as retention toh gadbad kyon aa rha hai

    // idx+1 is the number of features to be used
    int col = 0;
    int to_be_used = 0;
    float* reduced_eigen_vector_matrix = (float*)malloc(sizeof(float)*(idx+1)*M);
    
    // function to retrieve the corresponding eigen vectors
    for(index = 0;index<=idx;index++){
    	to_be_used = *(eigen_value_sorted_index+index);
    	for(col=0;col<M;col++){
    		*(reduced_eigen_vector_matrix+col*(1+idx)+index) = *(U+col*N+to_be_used);
    	}
    }
    
    *(D_HAT) = matrix_mult(M,N,M,1+idx,D,reduced_eigen_vector_matrix);
}
