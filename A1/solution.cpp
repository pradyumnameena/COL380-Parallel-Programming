#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <omp.h>

using namespace std;

void display(vector<float*> vec){
	int n = vec.size();
	for(int i = 0;i<n;i++){
		for(int j = 0;j<n;j++){
			cout << *(vec[i] + j) << " ";
		}
		cout << "" << endl;
	}
	cout << "***********" << endl;
}

vector<float*> matrix_mult(vector<float*> m1,vector<float*> m2){
	int n = m1.size();
	vector<float*> matrix;
	for(int i = 0;i<n;i++){
		float* address = (float*) malloc(n*sizeof(float));
		matrix.push_back(address);
		for(int j = 0;j<n;j++){
			*(matrix[i] + j) = 0;
			for(int k = 0;k<n;k++){
				*(matrix[i] + j) = *(matrix[i] + j) + (*(m1[i] + k) * *(m2[k] + j));
			}
		}
	}
	return matrix;
}

vector<float*> matrix_diff(vector<float*> m1,vector<float*> m2){
	int n = m1.size();
	vector<float*> matrix;
	for(int i = 0;i<n;i++){
		float* address = (float*) malloc(n*sizeof(float));
		matrix.push_back(address);
		for(int j = 0;j<n;j++){
			*(matrix[i] + j) = (*(m1[i] + j)) - (*(m2[i] + j));
		}
	}
	return matrix;
}

vector<float*> copy_matrix(vector<float*> m1){
	int n = m1.size();
	vector<float*> matrix;
	for(int i = 0;i<n;i++){
		float* address = (float*) malloc(n*sizeof(float));
		matrix.push_back(address);
		for(int j = 0;j<n;j++){
			*(matrix[i] + j) = *(m1[i] + j);
		}
	}
	return matrix;
}

float checker(vector<float*> m){
	float ans = 0.0;
	int n = m.size();
	vector<float> matrix;

	for(int i = 0;i<n;i++){
		matrix.push_back(0.0);
		for(int j = 0;j<n;j++){
			matrix[i] = matrix[i] + (*(m[j]+i)) * (*(m[j]+i));
		}
	}

	for(int i = 0;i<n;i++){
		ans+=sqrt(matrix[i]);
	}

	return ans;
}

vector<float*> read_data(){
	int n;
	vector<float*> matrix;
	string file_path = "./data.txt";

	ifstream file;
	file.open(file_path);
	if(file.is_open()){
		file >> n;
		for(int i = 0;i<n;i++){
			float* address = (float*) malloc(n*sizeof(float));
			matrix.push_back(address);
			for(int j = 0;j<n;j++){
				file >> *(matrix[i]+j);
			}
		}
	}
	
	file.close();
	return matrix;
}

vector<float*> initialize(int n){
	vector<float*> matrix;
	for(int i = 0;i<n;i++){
		float* address = (float*) malloc(n*sizeof(float));
		matrix.push_back(address);
		for(int j = 0;j<n;j++){
			*(matrix[i] + j) = 1.0*i*j;
		}
	}
	return matrix;
}

void LUdecomp(vector<float*> mat, vector<int>& perm, vector<float*> lower, vector<float*> upper){
	int temp;
	float max;
	float* add;
	int idx = -1;
	int n = mat.size();

	for(int k = 0;k<n;k++){
		max = 0.0;
		idx = k;
		for(int i = k;i<n;i++){
			if(abs(*(mat[i] + k))>max){
				max = abs(*(mat[i] + k));
				idx = i;
			}
		}

		if(max==0.0){
			return;
		}
		
		temp = perm[k];
		perm[k] = perm[idx];
		perm[idx] = temp;

		add = mat[k];
		mat[k] = mat[idx];
		mat[idx] = add;

		for(int i = 0;i<k;i++){
			max = *(lower[k] + i);
			*(lower[k] + i) = *(lower[idx] + i);
			*(lower[idx] + i) = max;
		}

		*(upper[k] + k) = *(mat[k] + k);

		for(int i = k+1;i<n;i++){
			*(lower[i] + k) = (*(mat[i] + k))/(*(upper[k] + k));
			*(upper[k] + i) = *(mat[k] + i);
		}

		for(int i = k+1;i<n;i++){
			for(int j = k+1;j<n;j++){
				*(mat[i] + j) = *(mat[i] + j) - (*(lower[i] + k))*(*(upper[k] + j));
			}
		}
	}

	for(int i = 0;i<n;i++){
		cout << perm[i] << endl;
	}
}

// void LUdecomp_OMP(vector<float*> mat, vector<int> perm, vector<float*> lower, vector<float*> upper,int num_threads){

// 	// return
// }

int main(int argc, char const *argv[]){
	vector<int> perm;
	vector<float*> lower;
	vector<float*> upper;
	vector<float*> identity;
	vector<float*> matrix = read_data();
	vector<float*> temp = copy_matrix(matrix);
	
	int n = matrix.size();
	for(int i = 0;i<n;i++){
		perm.push_back(i);
	}

	for(int i = 0;i<n;i++){
		float* add1 = (float*) malloc(n*sizeof(float));
		float* add2 = (float*) malloc(n*sizeof(float));
		float* add3 = (float*) malloc(n*sizeof(float));
		lower.push_back(add1);
		upper.push_back(add2);
		identity.push_back(add3);
	}

	for(int i = 0;i<n;i++){
		for(int j = 0;j<n;j++){
			*(upper[i] + j) = 0.0;
			*(lower[i] + j) = 0.0;
			*(identity[i] + j) = 0.0;
			if(i==j){
				*(lower[i] + j) = 1.0;
			}
		}
	}

	// int num_threads = 1;
	LUdecomp(temp,perm,lower,upper);
	// LUdecomp_OMP(matrix,perm,lower,upper,num_threads);

	for(int i = 0;i<n;i++){
		cout << perm[i] << endl;
		*(identity[i] + perm[i]) = 1.0;
	}

	
	vector<float*> prod1 = matrix_mult(identity,matrix);
	vector<float*> prod2 = matrix_mult(lower,upper);
	vector<float*> diff = matrix_diff(prod1,prod2);
	float ans = checker(diff);
	cout << ans << endl;
}