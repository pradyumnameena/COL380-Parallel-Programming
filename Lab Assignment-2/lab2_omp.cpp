#include <iostream>
#include <cstdlib>
#include <math.h>
#include <omp.h>
#include "lab2_omp.h"
using namespace std;

void matrix_mult(int m1,int n1,int m2,int n2,float* mat1,float* mat2,float* res_matrix){
	#pragma omp parallel for
	{
		for(int i = 0;i<m1*n2;i++){
			*(res_matrix+i) = 0;
		}
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

void matrix_mult_with_transpose(int m1,int n1,int m2,int n2,float* mat1,float* mat2,float* result_matrix){

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
				*(mat+i*M+i) = 1.0;
			}
		}
	}
	return;
}

void get_column_vector(float* matrix,float* col_vector,int length,int m){
	int idx = 0;
	for(idx = 0;idx<length;idx++){
		*(col_vector+idx) = *(matrix+(m-length+idx)*m + m-length);
	}
	return;
}

float get_norm(float* mat,int length){
	float rv = 0.0;
	int idx = 0;
	for(idx = 0;idx<length;idx++){
		rv+=((*(mat+idx))*(*(mat+idx)));
	}
	return sqrt(rv);
}

void matrix_scalar_mult(float* mat,float val,int num_elements){
	int idx = 0;
	#pragma omp parallel for
	{
		for(idx = 0;idx<num_elements;idx++){
			(*(mat+idx))*=val;
		}
	}
	return;
}

void subtract(float* mat1,float* mat2,int num_elements,int m){
	int i = 0;
	int j = 0;
	int length = int(sqrt(num_elements));
	int offset = length-m;

	for(i = m-length;i<m;i++){
		for(j = m-length;j<m;j++){
			*(mat1+i*m+j) -= (*(mat2+(offset+i)*length + (offset+j)));
		}
	}
	return;
}

void householders(int m,float* matrix,float* Q,float* R){
	int col = 0;
	float norm = 0;
	float* transpose;
	float* col_vector;
	float* matrix_prod;
	float* identity_mat = (float*)malloc(sizeof(float)*m*m);
	float* matrix_new = (float*)malloc(sizeof(float)*m*m);
	float* matrix_new_copy = (float*)malloc(sizeof(float)*m*m);

	float* Q_new = (float*)malloc(sizeof(float)*m*m);
	float* R_new = (float*)malloc(sizeof(float)*m*m);
	copy_matrix_values(matrix,R,m*m);
	copy_matrix_values(matrix,matrix_new,m*m);
	copy_matrix_values(matrix,matrix_new_copy,m*m);

	for(col = 0;col<m;col++){
		col_vector = (float*)malloc(sizeof(float)*(m-col));
		get_column_vector(matrix_new,col_vector,m-col,m);
		
		norm = get_norm(col_vector,m-col);
		if((*(col_vector))>0){
			*(col_vector) = (*(col_vector)) + norm;
		}else{
			*(col_vector) = (*(col_vector)) - norm;
		}

		norm = get_norm(col_vector,m-col);
		matrix_scalar_mult(col_vector,(1.0)/norm,(m-col));
		transpose = (float*)malloc(sizeof(float)*(m-col));
		matrix_prod = (float*)malloc(sizeof(float)*(m-col)*(m-col));
		matrix_transpose(m-col,1,col_vector,transpose);
		matrix_mult(m-col,1,1,m-col,col_vector,transpose,matrix_prod);

		identity_matrix(identity_mat,m);
		matrix_scalar_mult(matrix_prod,2.0,(m-col)*(m-col));
		subtract(identity_mat,matrix_prod,(m-col)*(m-col),m);
		
		free(transpose);
		free(col_vector);
		free(matrix_prod);

		matrix_mult(m,m,m,m,identity_mat,matrix_new,matrix_new_copy);
		copy_matrix_values(matrix_new_copy,matrix_new,m*m);
		matrix_mult(m,m,m,m,Q,identity_mat,Q_new);
		copy_matrix_values(Q_new,Q,m*m);
		matrix_mult(m,m,m,m,identity_mat,R,R_new);
		copy_matrix_values(R_new,R,m*m);
	}

	free(Q_new);
	free(R_new);
	free(identity_mat);
	free(matrix_new_copy);
	free(matrix_new);
	return;
}

float compute_error(float* mat1,float* mat2,int num){
	float rv = 0;
	int idx = 0;
	for(idx = 0;idx<num;idx++){
		*(mat2+idx) -= (*(mat1+idx));
		*(mat2+idx) *= (*(mat2+idx));
		rv += (*(mat2+idx));
	}
	return sqrt(rv);
}

void get_columns(float* mat,float* mat2,int length,int num_cols){
	int i = 0;
	int j = 0;

	for(i = 0;i<length;i++){
		for(j = 0;j<num_cols;j++){
			*(mat2+length*i+j) = *(mat+i*length+j);
		}
	}
	return;
}

void matrix_print(int rows,int cols,float* pointer){
	int i = 0;
	int j = 0;
	for(i = 0;i<rows;i++){
		for(j = 0;j<cols;j++){
			printf("%f, ",*(pointer+i*cols+j));
		}
		printf("\n");
	}
}

void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{	
	float* D_transpose = (float*)malloc(sizeof(float)*N*M);
	float* matrix_prod = (float*)malloc(sizeof(float)*N*N);
	matrix_transpose(M,N,D,D_transpose);
	matrix_mult(N,M,M,N,D_transpose,D,matrix_prod);

	float* QR_decomposition_Q = (float*)malloc(sizeof(float)*N*N);
	float* QR_decomposition_R = (float*)malloc(sizeof(float)*N*N);
	float* D_original = (float*)malloc(sizeof(float)*N*N);
	float* D_old = (float*)malloc(sizeof(float)*N*N);
	float* E_original = (float*)malloc(sizeof(float)*N*N);
	float* E_new = (float*)malloc(sizeof(float)*N*N);

	int num_iters = 0;
	float epsilon = 0.001;
	float error = 10;
	int max_iterations = 30000;

	copy_matrix_values(matrix_prod,D_original,N*N);
	free(matrix_prod);
	identity_matrix(E_original,N);
	identity_matrix(E_new,N);
	
	while(error>epsilon && num_iters<max_iterations){
		identity_matrix(QR_decomposition_Q,N);
		identity_matrix(QR_decomposition_R,N);
		householders(N,D_original,QR_decomposition_Q,QR_decomposition_R);
		copy_matrix_values(D_original,D_old,N*N);
		matrix_mult(N,N,N,N,QR_decomposition_R,QR_decomposition_Q,D_original);
		matrix_mult(N,N,N,N,E_original,QR_decomposition_Q,E_new);
		copy_matrix_values(E_new,E_original,N*N);
		error = compute_error(D_original,D_old,N*N);
		num_iters++;
	}
	printf("Number of iterations = %d\n",num_iters);
	printf("Error = %f\n",error);

	printf("Eigen vector matrix\n");
	matrix_print(N,N,E_new);
	printf("Eigen values matrix\n");
	matrix_print(N,N,D_original);
	printf("*************\n");

	free(QR_decomposition_Q);
	free(QR_decomposition_R);
	free(E_original);
	free(D_old);
	
	int idx = 0;
	float* eigen_values = (float*)malloc(sizeof(float)*N);
    int* eigen_value_sorted_index = (int*)malloc(sizeof(int)*N);
    
	for(idx = 0;idx<N;idx++){
		*(D_original+idx*N+idx) = sqrt(abs((*(D_original+idx*N+idx))));
		*(eigen_values+idx) = *(D_original+idx*N+idx);
		*(eigen_value_sorted_index+idx) = idx;
	}

	int i = 0;
	int j = 0;
    float temp = 0;
    for(i = 0;i<N-1;i++){
    	for(j = i+1;j<N;j++){
    		if((*(eigen_values+j)) > (*(eigen_values+i))){
    			temp = *(eigen_values+j);
    			*(eigen_values+j) = *(eigen_values+i);
    			*(eigen_values+i) = temp;
    			*(eigen_value_sorted_index+i) = j;
    			*(eigen_value_sorted_index+j) = i;
    		}
    	}
    }
    
    idx = 0;
    int diag_idx = 0;
    float* eigen_values_pointer = *(SIGMA);
    float* Sigma_inverse = (float*)malloc(sizeof(float)*M*N);

    for(idx = 0;idx<N;idx++){
		*(eigen_values_pointer+idx) = *(eigen_values+idx);
	}

	free(eigen_values);

    idx = 0;
	for(idx = 0;idx<N*M;idx++){
		*(Sigma_inverse+idx) = 0;
	}
	
	for(diag_idx = 0;diag_idx<N;diag_idx++){
		*(Sigma_inverse+diag_idx*N+diag_idx) = 1/(*(eigen_values_pointer+diag_idx));
	}
	
	idx = 0;
    int row = 0;
    int index = 0;
    float* eigen_vector_pointer = *(U);
    for(idx = 0;idx<N;idx++){
    	index = *(eigen_value_sorted_index+idx);
    	for(row = 0;row<N;row++){
    		*(eigen_vector_pointer+row*N+idx) = *(E_new+row*N+index);
    	}
    }

	free(eigen_value_sorted_index);
	float* U_transpose = (float*)malloc(sizeof(float)*N*N);
    matrix_transpose(N,N,eigen_vector_pointer,U_transpose);
    matrix_mult(M,N,N,N,Sigma_inverse,U_transpose,Sigma_inverse);
    free(U_transpose);
    matrix_mult(M,N,N,M,Sigma_inverse,D_transpose,*(V_T));
    free(D_transpose);
    free(Sigma_inverse);
	return;
}

void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{	
	int idx = 0;
    float sum = 0;
    float retention_reqd = retention*0.01;
    float* reduced_variance = (float*)malloc(sizeof(float)*N);

    // sum of eigen values
    for(idx = 0;idx<N;idx++){
    	sum+=(*(SIGMA+idx));
    }

    // reduction variances
    #pragma omp parallel for
    {
    	for(idx = 0;idx<N;idx++){
	    	*(reduced_variance+idx) = (*(SIGMA+idx))/sum;
	    }
    }

    sum = 0.0;
    for(idx = 0;idx<N;idx++){
    	sum+=(*(reduced_variance+idx));
    	if(sum>=retention_reqd){
    		break;
    	}
    }

    if(retention==100){
    	idx = N-1;
    }
    *K = 1+idx;

    free(reduced_variance);
    float* W = (float*)malloc(sizeof(float)*(1+idx)*N);
    get_columns(U,W,N,1+idx);
    float* D_HAT_pointer = (float*)malloc(sizeof(float)*(1+idx)*M);
    *D_HAT = D_HAT_pointer;
    matrix_mult(M,N,N,1+idx,D,W,D_HAT_pointer);
    free(W);
    return;
}



