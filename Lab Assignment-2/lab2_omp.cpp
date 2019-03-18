#include <iostream>
#include <cstdlib>
#include <math.h>
#include <omp.h>
#include "lab2_omp.h"
using namespace std;

void matrix_mult(int m1,int n1,int m2,int n2,float* mat1,float* mat2,float* res_matrix){
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
	return;
}

void matrix_transpose(int m,int n,float* mat,float* res_matrix){
	for(int i = 0;i<m;i++){
		for(int j = 0;j<n;j++){
			*(res_matrix+j*m+i) = *(mat+i*n+j);
		}
	}
	return;
}

void copy_matrix_values(float* mat1,float* mat2,int num_vals){
	#pragma omp parallel for
	{
		for(int i = 0;i<num_vals;i++){
			*(mat2+i) = *(mat1+i);
		}
	}
	return;
}

void identity_matrix(float* mat,int M){
	int i = 0;
	#pragma omp parallel
	{
		#pragma omp for
		{
			for(i = 0;i<M*M;i++){
				*(mat+i) = 0;
			}
		}

		#pragma omp for
		{
			for(i = 0;i<M;i++){
				*(mat+i*M+i) = 1;
			}
		}
	}
	return;
}

void givens_rotation(int m,float* matrix,float* Q,float* R){
	int i = 0;
	int j = 0;
	float c = 0;
	float s = 0;
	float* identity = (float*)malloc(sizeof(float)*m*m);
	float* Q_new = (float*)malloc(sizeof(float)*m*m);
	float* R_new = (float*)malloc(sizeof(float)*m*m);
	float* identity_transpose = (float*)malloc(sizeof(float)*m*m);
	identity_matrix(Q,m);
	identity_matrix(identity,m);
	copy_matrix_values(matrix,R,m*m);

	for(i = 0;i<m-1;i++){
		for(j = m-1;j>i;j--){
			s = *(R+j*m+i);
			c = *(R+(j-1)*m+i);
			s = s/(sqrt(c*c + s*s));
			c = c/(sqrt(c*c + s*s));
			*(identity+j*m+i) = c;
			*(identity+j*m+i+1) = s;
			*(identity+(j+1)*m+i) = -1*s;
			*(identity+(j+1)*m+i+1) = c;
			matrix_mult(m,m,m,m,identity,R,R_new);
			copy_matrix_values(R_new,R,m*m);
			matrix_transpose(m,m,identity,identity_transpose);
			matrix_mult(m,m,m,m,Q,identity_transpose,Q_new);
			copy_matrix_values(Q_new,Q,m*m);
		}
	}

	free(Q_new);
	free(R_new);
	free(identity);
	free(identity_transpose);
	return;
}

void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{	
	float* D_transpose = (float*)malloc(sizeof(float)*N*M);
	float* matrix_prod = (float*)malloc(sizeof(float)*M*M);
	matrix_transpose(M,N,D,D_transpose);
	matrix_mult(M,N,N,M,D,D_transpose,matrix_prod);
	free(D_transpose);
	
	float* QR_decomposition_Q = (float*)malloc(sizeof(float)*M*M);
	float* QR_decomposition_R = (float*)malloc(sizeof(float)*M*M);
	float* D_original = (float*)malloc(sizeof(float)*M*M);
	float* E_original = (float*)malloc(sizeof(float)*M*M);
	float* E_new = (float*)malloc(sizeof(float)*M*M);
	
	int row = 0;
	int col = 0;
	int num_iters = 0;
	bool condition = true;
	int iterations_limit = 1;
	
	copy_matrix_values(matrix_prod,D_original,M*M);
	identity_matrix(E_original,M);
	identity_matrix(E_new,M);

	while(condition){
		givens_rotation(M,D_original,QR_decomposition_Q,QR_decomposition_R);
		matrix_mult(M,M,M,M,QR_decomposition_R,QR_decomposition_Q,D_original);
		matrix_mult(M,M,M,M,E_original,QR_decomposition_Q,E_new);
		copy_matrix_values(E_new,E_original,M*M);
		condition = (num_iters<iterations_limit-1);
		num_iters++;
	}
	return;
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
    
    // trivial/corner case
    if(retention==100){
    	idx = N-1;
    }

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
 
    matrix_mult(M,N,M,1+idx,D,reduced_eigen_vector_matrix,*(D_HAT));
}