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

struct threadArgs1{
	int col;
	int idx;
	int end;
	int start;
	vector<double*> *array;
};

struct threadArgs2{
	int col;
	int idx;
	int end;
	int start;
	vector<double*> *lower;
	vector<double*> *upper;
	vector<double*> *mat;
};

void *func1(void* threadarg){
	struct threadArgs1 *args = (struct threadArgs1*)threadarg;
	int col = args->col;
	int idx = args->idx;
	int end = args->end;
	int start = args->start;
	double temp;
	
	for(int i = start;i<end;i++){
		temp = *((*args->array)[col] + i);
		*((*args->array)[col] + i) = *((*args->array)[idx] + i);
		*((*args->array)[idx] + i) = temp;
	}
	
	pthread_exit(0);
}

void *func2(void* threadarg){
	struct threadArgs2 *args = (struct threadArgs2*)threadarg;
	int col = args->col;
	int idx = args->idx;
	int end = args->end;
	int start = args->start;
	
	for(int i = start;i<end;i++){
		*((*args->lower)[i] + col) = *((*args->mat)[i] + col) / (*((*args->upper)[col] + col));
		*((*args->upper)[col] + i) = *((*args->mat)[col] + i);
	}
	
	pthread_exit(0);
}

void *func3(void* threadarg){
	struct threadArgs2 *args = (struct threadArgs2*)threadarg;
	int col = args->col;
	int idx = args->idx;
	int end = args->end;
	int start = args->start;
	
	for(int i = start;i<end;i++){
		for(int j = col+1;j<idx;j++){
			*((*args->mat)[i] + j) = *((*args->mat)[i] + j) - (*((*args->lower)[i] + col))*(*((*args->upper)[col] + j));
		}
	}
	
	pthread_exit(0);
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

vector<double*> read_data(string file_path,int n){
	int i,j;
	vector<double*> matrix;

	ifstream file;
	file.open(file_path);
	if(file.is_open()){
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

void write_data(vector<double*> matrix,string name){
	int n = matrix.size();
	ofstream file;
	file.open("./" + name);
	
	if(file.is_open()){
		for(int i = 0;i<n;i++){
			for(int j = 0;j<n;j++){
				file << *(matrix[i] + j);
				if(j==n-1){
					if(i!=n-1){
						file << "\n";
					}
				}else{
					file << " ";
				}
			}
		}
	}

	file.close();
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
			for(k = 0;k<n;k++){
				v1-=((*(lower[i] + k)) * (*(upper[k] + j)));
			}
			temp+=(v1*v1);
		}
		ans+=sqrt(temp);
	}

	return ans;
}

void LUdecomp(vector<double*> mat, vector<int>& perm, vector<double*> lower, vector<double*> upper){
	double max;
	double* add;
	int i,j,n,temp,idx,col,res;
	int bounds,part;
	
	n = mat.size();
	
	pthread_t threads[thread_count];
	struct threadArgs1 threadArr1[thread_count];
	struct threadArgs2 threadArr2[thread_count];

	for(int i = 0;i<thread_count;i++){
		threadArr1[i].array = &lower;
		threadArr2[i].lower = &lower;
		threadArr2[i].upper = &upper;
		threadArr2[i].mat = &mat;
	}
	
	for(col = 0;col<n;col++){
		res = 0;
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

		part = col/thread_count;
		bounds = 0;

		for(int thread_id = 0;thread_id<thread_count;thread_id++){
			threadArr1[thread_id].col = col;
			threadArr1[thread_id].idx = idx;
			threadArr1[thread_id].start = bounds;
			threadArr1[thread_id].end = bounds + part;
			if(thread_id==thread_count-1){
				threadArr1[thread_id].end = col;
			}
			bounds+=part;
			res = pthread_create(&threads[thread_id],NULL,func1,&threadArr1[thread_id]);
		}

		for(int thread_id = 0;thread_id<thread_count;thread_id++){
			pthread_join(threads[thread_id],NULL);
		}

		part = (n-1-col)/thread_count;
		bounds = col+1;

		for(int thread_id = 0;thread_id<thread_count;thread_id++){
			threadArr2[thread_id].col = col;
			threadArr2[thread_id].idx = idx;
			threadArr2[thread_id].start = bounds;
			threadArr2[thread_id].end = bounds + part;
			if(thread_id==thread_count-1){
				threadArr2[thread_id].end = n;
			}
			bounds+=part;
			res = pthread_create(&threads[thread_id],NULL,func2,&threadArr2[thread_id]);
		}

		for(int thread_id = 0;thread_id<thread_count;thread_id++){
			pthread_join(threads[thread_id],NULL);
		}
		
		part = (n-1-col)/thread_count;
		bounds = col+1;

		for(int thread_id = 0;thread_id<thread_count;thread_id++){
			threadArr2[thread_id].col = col;
			threadArr2[thread_id].idx = n;
			threadArr2[thread_id].start = bounds;
			threadArr2[thread_id].end = bounds + part;
			if(thread_id==thread_count-1){
				threadArr2[thread_id].end = n;
			}
			bounds+=part;
			res = pthread_create(&threads[thread_id],NULL,func3,&threadArr2[thread_id]);
		}

		for(int thread_id = 0;thread_id<thread_count;thread_id++){
			pthread_join(threads[thread_id],NULL);
		}
	}
}

int main(int argc, char const *argv[]){
	vector<int> perm;
	vector<double*> lower;
	vector<double*> upper;
	vector<double*> perm_final;
	struct timeval start, end;
	double time_taken,time_taken2;
	
	int n = atoi(argv[1]);
	thread_count = atoi(argv[2]);
	int check = atoi(argv[3]);
	string file_path = argv[4];
	int write = atoi(argv[5]);

	vector<double*> matrix = read_data(file_path,n);
	// vector<double*> matrix = initialize(n);
	vector<double*> matrix_cp = matrix_copy(matrix);
	
	int i,j,k;

	for(i = 0;i<n;i++){
		perm.push_back(i);
	}

	for(i = 0;i<n;i++){
		double* add1 = (double*)malloc(n*sizeof(double));
		double* add2 = (double*)malloc(n*sizeof(double));
		double* add3 = (double*)malloc(n*sizeof(double));
		lower.push_back(add1);
		upper.push_back(add2);
		perm_final.push_back(add3);
	}

	for(i = 0;i<n;i++){
		for(j = 0;j<n;j++){
			*(upper[i] + j) = 0.0;
			*(lower[i] + j) = 0.0;
			*(perm_final[i] + j) = 0.0;
		}
		*(lower[i] + i) = 1.0;
	}

	gettimeofday(&start, NULL);
	LUdecomp(matrix_cp,perm,lower,upper);
	gettimeofday(&end, NULL);

	time_taken = (end.tv_sec - start.tv_sec) * 1e6;
	time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
	cout << "LU Decomposition in " << time_taken << " sec" << endl;

	cout << check << endl;

	if(check==1){
		gettimeofday(&start, NULL);
		double ans = checker(perm,matrix,lower,upper);
		gettimeofday(&end, NULL);
		
		time_taken2 = (end.tv_sec - start.tv_sec) * 1e6;
		time_taken2 = (time_taken2 + (end.tv_usec - start.tv_usec)) * 1e-6;
		
		cout << "L2,1 norm = " << ans << endl;
		cout << "Time taken for checking : " << time_taken2 << " sec" << endl;
	}

	// Writing into files
	
	if(write==1){
		for(int i = 0;i<n;i++){
			*(perm_final[i] + perm[i]) = 1.0;
		}

		string P_name = "P_" + to_string(n) + "_" + to_string(thread_count) + ".txt";
		string L_name = "L_" + to_string(n) + "_" + to_string(thread_count) + ".txt";
		string U_name = "U_" + to_string(n) + "_" + to_string(thread_count) + ".txt";
		
		write_data(perm_final,P_name);
		write_data(lower,L_name);
		write_data(upper,U_name);
	}
}