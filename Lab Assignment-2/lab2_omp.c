// #include <malloc.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

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

		// computing the eigen values and sorting them in descending order

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
