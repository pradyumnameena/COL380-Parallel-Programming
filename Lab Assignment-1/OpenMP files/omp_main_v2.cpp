#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <omp.h>
#include <string.h>
using namespace std;

class point{
	public: 
		vector<int> coordinates;
		int cluster_index;

	point(vector<int> coord,int idx){
		coordinates = coord;
		cluster_index = idx;
	}

	void change_index(int nidx){
		cluster_index = nidx;
	}

	double compute_dist(vector<double> centroid){
		double dist = 0.0;
		int idx = 0;
		float diff = 0.0;
		for(idx = 0;idx<=2;idx++){
			diff = double(coordinates.at(idx)) - centroid.at(idx);
			dist+=(diff*diff);
		}
		return dist;
	}
};

class centroid{
	public: 
		vector<double> coordinates;
		int cluster_index;

	centroid(vector<double> coord,int idx){
		coordinates = coord;
		cluster_index = idx;
	}

	void change_index(int nidx){
		cluster_index = nidx;
	}
};

void initialize(vector<point*> points,vector<centroid*> centroids){
	ifstream read_file;
	string file_path = "/Volumes/Macintosh_HD2/IIT_Stuff/Academics/6th sem/COL380/Assignments/Lab-1/data.txt";
	read_file.open(file_path);
	char output[100];
	int n,d,k;
    char *tok;
    vector<int> p1;
    
    if(read_file.is_open()){
        read_file.getline(output,100);
		n = stoi(output);
		read_file.getline(output,100);
		d = stoi(output);
		read_file.getline(output,100);
		k = stoi(output);
        
		while(!read_file.eof()){
			read_file.getline(output,100);
            tok = strtok(output," ");
            while(tok!=NULL){
                p1.push_back(stoi(tok));
                tok = strtok(NULL," ");
            }
            point* p = new point(p1,0);
            points.push_back(p);
            p1.clear();
		}
        read_file.close();
	}

	// initializing the clusters
	return;
}

int main(){
	vector<point*> points;
	vector<centroid*> centroids;
	initialize(points,centroids); 
	int num_points = points.size();
	int num_centroids = centroids.size();
	
	// vector<int> a;
	// vector<double> b;
	// a.push_back(1);
	// a.push_back(2);
	// a.push_back(3);
	// b.push_back(1.000001);
	// b.push_back(2.000001);
	// b.push_back(3.000001);
	// point a1(a,0);
	// centroid b1(b,0);

	// parallel part
	int num_threads = 4; // initially same as number of cores on the system

	#pragma omp parallel
	{	
		// cout << omp_get_thread_num() << endl;
		// cout << a1.compute_dist(b1.coordinates) << endl;
		int tid = omp_get_thread_num();
		int start = 
	}
	return 0;
}