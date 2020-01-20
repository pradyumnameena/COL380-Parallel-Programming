#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

void display(vector<float*> vec){
	int n = vec.size();
	for(int i = 0;i<n;i++){
		for(int j = 0;j<n;j++){
			cout << *(vec[i] + j) << " ";
		}
		cout << "" << endl;
	}
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

float checker(vector<float*> m){
	float ans = 0.0;
	int n = m.size();
	vector<float> matrix;

	for(int i = 0;i<n;i++){
		matrix.push_back(0.0);
		for(int j = 0;j<n;j++){
			matrix[i] = matrix[i] + (*(m[i]+j)) * (*(m[i]+j));
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

int main(int argc, char const *argv[]){
	vector<float*> vec1 = read_data();
	vector<float*> vec2 = read_data();
	vector<float*> vec3 = matrix_mult(vec1, vec2);
	display(vec3);
}