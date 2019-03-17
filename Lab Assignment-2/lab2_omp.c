// #include <malloc.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

void matrix_mult(int m1,int n1,int m2,int n2,int* mat1,int* mat2,int* result);

void matrix_mult(int m1,int n1,int m2,int n2,int* mat1,int* mat2,int* result){
	if(n1!=m2){
		printf("Different dimensions. Product not possible\n");
		return;
	}

	// initializing
	for(int i = 0;i<m1*n2*n1;i++){
		*(result + i) = 0;
	}

	// computing
	for(int i = 0;i<m1;i++){
		for(int j = 0;j<n2;j++){
			for(int k = 0;k<n1;k++){
				*(result + i*n2 + k) += (*(mat1 + i*n1 + k)) * (*(mat2 + j*n1 + k));
			}
		}
	}
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{
	float* matrix_prod = (float*)malloc(sizeof(float)*(M*M));
	int i = 0;
	int j = 0;
	int k = 0;
	int offset1 = 0;
	int offset2 = 0;
	int offset3 = 0;
	int thread_count = 4;
	int num_vals = M*M;
	float val = 0.0;
	int* V_T_pointer = *V_T;
	int* sigma_pointer = *SIGMA;

	// For QR algorithm
	int* gram_schmidt_Q = (float*)malloc(sizeof(float)*(M*M));
	int* gram_schmidt_R = (float*)malloc(sizeof(float)*(M*M));
	int* gram_schmidt_E = (float*)malloc(sizeof(float)*(M*M));
	int* gram_schmidt_D = (float*)malloc(sizeof(float)*(M*M));


	#pragma omp parallel shared(matrix_prod) num_threads(thread_count)
	{	
		// initializing memory for dot product values
		#pragma omp for
		{
			for(i = 0;i<num_vals;i++){
				*(matrix_prod + i) = 0;
			}
		}

		// computing the dot product
		#pragma omp for private(i,j,k,offset1,offset2,offset3,val)
		{
			for(i = 0;i<M;i++){
				offset1 = i*M;
				offset2 = i*N;
				for(j = 0;j<M;j++){
					offset3 = j*N;
					for(k = 0;k<N;k++){
						val += (*(D+offset2+k)) * (*(D+offset2+k));
					}
					*(matrix_prod + offset1 + j) = val;
				}
			}
		}		
	}

	// computing the eigen values and sorting them in descending order
		while(true){
			for(i = 0;i<M;i++){
				for(j = 0;j<M;j++){

				}
			}
		}
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    
}
