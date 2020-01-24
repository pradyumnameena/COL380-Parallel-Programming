#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>
#include <pthread.h>

using namespace std;
double limit = 1000;
int thread_count = 1;

struct threadArgs{
int tid;
int end;
int start;
int col;
int idx;
vector<double*> *array;
};

struct threadArgs2{
int tid;
int end;
int start;
int col;
int idx;
vector<double*> *l;
vector<double*> *u;
vector<double*> *mat;
};

void *func(void* threadarg){
	struct threadArgs *args = (struct threadArgs*) threadarg;
	int tid = args->tid;
	int max, col, idx, i;
	col = args->col,idx = args->idx;
	for(i = start; i < end; i++){
		// max = *(lower[col] + (i);
		// *(lower[col] + i) = *(lower[idx] + i);
		// *(lower[idx] + i) = max;
		max = *(args->array[col] + i);
		*(args->array[col] + i) = *(args->array[idx] + i);
		*(args->array[idx] + i) = max;
	}
}

void *func2(void* threadarg){
	struct threadArgs* args = (struct threadArgs2*) threadarg;
	int tid = args->tid;
	int col, idx, i;
	col = args->col,idx = args->idx;
	for(i = start; i < end; i++){
		// *(lower[i] + col) = (*(mat[i] + col))/(*(upper[col] + col));
		// *(upper[col] + i) = *(mat[col] + i);
		*(args->l[i] + col) = *(args->mat[i] + col) / (*(args->u[col] + col));
		*(args->u[col] + i) = *(args->mat[col] + i);
	}
}

void *func3(void* threadarg){
	struct threadArgs* args = (struct threadArgs2*) threadarg;
	int tid = args->tid;
	int col, idx, i, j;
	col = args->col,idx = args->idx;
	for(i = start; i < end; i++){
		for(j = start; j < end; j++){
			// *(mat[i] + j) = *(mat[i] + j) - (*(lower[i] + col))*(*(upper[col] + j));
			*(args->mat[i] + j) = *(args->mat[i] + j) - (*(args->l[i] + col)) * (*(args->upper[col] + j));
		}
	}
}

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

double checker(vector<int>& perm, vector<double*> mat, vector<double*> lower, vector<double*> upper){
	double v1,v2,temp,ans;
	int i,j,k,n;

	n = mat.size();
	ans = 0.0;
	for(j = 0;j<n;j++){
		temp = 0.0;
		for(i = 0;i<n;i++){
			v1 = *(mat[perm[i]] + j);
			v2 = 0.0;
			for(k = 0;k<n;k++){
				v2+=((*(lower[i] + k)) * (*(upper[k] + j)));
			}
			temp+=((v2-v1)*(v2-v1));
		}
		ans+=sqrt(temp);
	}

	return ans;
}

void LUdecomp(vector<double*> mat, vector<int>& perm, vector<double*> lower, vector<double*> upper){
	double max;
	double* add;
	int i,j,n,temp,idx,col;
	
	n = mat.size();

	for(col = 0;col<n;col++){
		cout << col << endl;
		max = 0.0;
		idx = col;
		for(i = col;i<n;i++){
			if(fabs(*(mat[i] + col))>max){
				max = fabs(*(mat[i] + col));
				idx = i;
			}
		}

		if(max==0.0){
			cout << "Singular Matrix" << endl;
			return;
		}
		
		temp = perm[col];
		perm[col] = perm[idx];
		perm[idx] = temp;

		add = mat[col];
		mat[col] = mat[idx];
		mat[idx] = add;

		*(upper[col] + col) = *(mat[col] + col);

		for(i = 0;i<col;i++){
			max = *(lower[col] + i);
			*(lower[col] + i) = *(lower[idx] + i);
			*(lower[idx] + i) = max;
		}

		for(i = col+1;i<n;i++){
			*(lower[i] + col) = (*(mat[i] + col))/(*(upper[col] + col));
			*(upper[col] + i) = *(mat[col] + i);
		}

		for(i = col+1;i<n;i++){
			for(j = col+1;j<n;j++){
				*(mat[i] + j) = *(mat[i] + j) - (*(lower[i] + col))*(*(upper[col] + j));
			}
		}
	}
}

int main(int argc, char const *argv[]){
	vector<int> perm;
	vector<double*> lower;
	vector<double*> upper;
	struct timeval start, end;
	double time_taken,time_taken2;
	
	int n = atoi(argv[1]);
	thread_count = atoi(argv[2]);
	int check = atoi(argv[3]);

	// vector<double*> matrix = read_data();
	vector<double*> matrix = initialize(n);
	vector<double*> matrix_cp = matrix_copy(matrix);
	
	int i,j,k;

	for(i = 0;i<n;i++){
		perm.push_back(i);
	}

	for(i = 0;i<n;i++){
		double* add1 = (double*)malloc(n*sizeof(double));
		double* add2 = (double*)malloc(n*sizeof(double));
		lower.push_back(add1);
		upper.push_back(add2);
	}

	for(i = 0;i<n;i++){
		for(j = 0;j<n;j++){
			*(upper[i] + j) = 0.0;
			*(lower[i] + j) = 0.0;
		}
		*(lower[i] + i) = 1.0;
	}

	gettimeofday(&start, NULL);
	LUdecomp(matrix_cp,perm,lower,upper);
	gettimeofday(&end, NULL);

	time_taken = (end.tv_sec - start.tv_sec) * 1e6;
	time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
	cout << "Time taken by program is : " << time_taken << " sec" << endl; 

	if(check==1){
		gettimeofday(&start, NULL);
		double ans = checker(perm,matrix,lower,upper);
		gettimeofday(&end, NULL);
		
		time_taken2 = (end.tv_sec - start.tv_sec) * 1e6;
		time_taken2 = (time_taken2 + (end.tv_usec - start.tv_usec)) * 1e-6;
		
		cout << "L2,1 norm = " << ans << endl;
		cout << "Time taken for checking : " << time_taken2 << " sec" << endl;
	}
}