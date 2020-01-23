#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>

using namespace std;
float limit = 10000;

// For random initialization
vector<double*> initialize(int n){
	int i,j;
	vector<double*> matrix;

	for(i = 0;i<n;i++){
		double* address = (double*) malloc(n*sizeof(double));
		matrix.push_back(address);
		for(j = 0;j<n;j++){
			*(matrix[i] + j) = drand48() * limit;
		}
	}
	return matrix;
}

/* Reads matrix from a file
	Format -
	 1. first line should contain n - the dimension of matrix
	 2. next n lines should contain n space separated matrix values */
vector<double*> read_data(){
	int i,j,n;
	vector<double*> matrix;
	string file_path = "./data.txt";

	ifstream file;
	file.open(file_path);
	if(file.is_open()){
		file >> n;
		for(i = 0;i<n;i++){
			double* address = (double*) malloc(n*sizeof(double));
			matrix.push_back(address);
			for(j = 0;j<n;j++){
				file >> *(matrix[i]+j);
			}
		}
	}
	
	file.close();
	return matrix;
}

// Copies a matrix
vector<double*> matrix_copy(vector<double*> m1){
	int i,j,n;
	n = m1.size();
	vector<double*> matrix;
	for(i = 0;i<n;i++){
		double* address = (double*) malloc(n*sizeof(double));
		matrix.push_back(address);
		for(j = 0;j<n;j++){
			*(matrix[i] + j) = *(m1[i] + j);
		}
	}
	return matrix;
}

// Matrix multiplication for square matrices
vector<double*> matrix_mult(vector<double*> m1,vector<double*> m2){
	int i,j,k,n;
	double temp;
	vector<double*> matrix;
	
	n = m1.size();
	for(i = 0;i<n;i++){
		double* address = (double*) malloc(n*sizeof(double));
		matrix.push_back(address);
		for(j = 0;j<n;j++){
			temp = 0.0;
			for(k = 0;k<n;k++){
				temp += (*(m1[i] + k) * *(m2[k] + j)); 
			}
			*(matrix[i] + j) = temp;
		}
	}
	return matrix;
}

// Displaying any square matrix
void display(vector<double*> vec){
	int i,j,n;
	n = vec.size();
	for(i = 0;i<n;i++){
		for(j = 0;j<n;j++){
			cout << *(vec[i] + j) << " ";
		}
		cout << endl;
	}
	cout << "**************THE END**************" << endl;
}

// Calculates L2,1 norm for m1-m2 matrix
double checker(vector<double*> m1, vector<double*> m2){
	double v1,v2,temp;
	double ans = 0.0;
	int i,j,n;
	n = m1.size();

	for(i = 0;i<n;i++){
		temp = 0.0;
		for(j = 0;j<n;j++){
			v1 = *(m1[j] + i);
			v2 = *(m2[j] + i);
			temp+=((v1-v2)*(v1-v2));
		}
		ans+=sqrt(temp);
	}

	return ans;
}

// LU decomposition function
void LUdecomp(vector<double*> mat, vector<int>& perm, vector<double*> lower, vector<double*> upper){
	double max;
	double* add;
	int i,j,n,temp,idx,col;
	
	n = mat.size();

	// Loop on columns
	for(col = 0;col<n;col++){
		// Find index of row with maximum entry in the column
		max = 0.0;
		idx = col;
		for(i = col;i<n;i++){
			if(fabs(*(mat[i] + col))>max){
				max = abs(*(mat[i] + col));
				idx = i;
			}
		}

		// Error -> Singular Matrix
		if(max==0.0){
			cout << "Singular Matrix" << endl;
			return;
		}
		
		// Swapping in permutation vector
		temp = perm[col];
		perm[col] = perm[idx];
		perm[idx] = temp;

		// Swapping the rows
		add = mat[col];
		mat[col] = mat[idx];
		mat[idx] = add;

		// L[col,1:col-1] <-> L[idx,1:col-1]
		for(i = 0;i<col;i++){
			max = *(lower[col] + i);
			*(lower[col] + i) = *(lower[idx] + i);
			*(lower[idx] + i) = max;
		}

		// U[col,col] = A[col,col]
		*(upper[col] + col) = *(mat[col] + col);

		/*
			L[i,col] = A[i,col]/U[col,col]
			U[col,i] = A[col,i]
		*/
		for(i = col+1;i<n;i++){
			*(lower[i] + col) = (*(mat[i] + col))/(*(upper[col] + col));
			*(upper[col] + i) = *(mat[col] + i);
		}

		// A[i,j] -= (L[i,k]*U[k,j])
		for(i = col+1;i<n;i++){
			for(j = col+1;j<n;j++){
				*(mat[i] + j) = *(mat[i] + j) - (*(lower[i] + col))*(*(upper[col] + j));
			}
		}
	}
}

int main(int argc, char const *argv[]){
	// Permutation vector, lower and upper matrix, permutation matrix, Matrix, mMatrix copy
	vector<int> p;
	vector<double*> lower;
	vector<double*> upper;
	vector<double*> perm;
	struct timeval start, end;
	double time_taken;
	
	vector<double*> matrix = read_data();
	// vector<double*> matrix = initialize(100);
	vector<double*> matrix_cp = matrix_copy(matrix);
	
	int i,j,k;
	int n = matrix.size();

	// P[i] = i
	for(i = 0;i<n;i++){
		p.push_back(i);
	}

	// Allocating the memory and inserting into vectors
	for(i = 0;i<n;i++){
		double* add1 = (double*)malloc(n*sizeof(double));
		double* add2 = (double*)malloc(n*sizeof(double));
		double* add3 = (double*)malloc(n*sizeof(double));
		lower.push_back(add1);
		upper.push_back(add2);
		perm.push_back(add3);
	}

	// Initialising all the matrices
	for(i = 0;i<n;i++){
		for(j = 0;j<n;j++){
			*(upper[i] + j) = 0.0;
			*(lower[i] + j) = 0.0;
			*(perm[i] + j) = 0.0;
		}
		*(lower[i] + i) = 1.0;
	}

	gettimeofday(&start, NULL);
	LUdecomp(matrix_cp,p,lower,upper);
	gettimeofday(&end, NULL);

	time_taken = (end.tv_sec - start.tv_sec) * 1e6;
	time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;

	// Generating the permutation matrix
	for(i = 0;i<n;i++){
		*(perm[i] + p[i]) = 1.0;
	}

	// Calculating the L2,1 norm
	vector<double*> PA = matrix_mult(perm,matrix);
	vector<double*> LU = matrix_mult(lower,upper);
	double ans = checker(PA,LU);
	cout << "L2,1 norm = " << ans << endl;
	cout << "Time taken by program is : " << time_taken << " sec" << endl; 
}