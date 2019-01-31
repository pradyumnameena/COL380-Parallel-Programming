#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <list>
#include <map>
#include <math.h>
#include <string.h>
using namespace std;

float distance(vector<int> point,vector<float> centroid){
	float dist = 0.0;
	int idx = 0;
	int size = point.size();
	float diff = 0.0;

	for(idx = 0;idx<size;idx++){
		diff = (float)point.at(idx) - centroid.at(idx);
		dist+=(diff*diff);
	}
	return sqrt(dist);
}

vector<float> add(vector<float> p1,vector<int> p2){
	int idx = 0;
	int size = p1.size();
	vector<float> rv;
	for(idx = 0;idx<size;idx++){
		rv.push_back(p1.at(idx)+(float)p2.at(idx));
	}
	return rv;
}

vector<float> divide(vector<float> p,int m){
	int idx = 0;
	int size = p.size();
	vector<float> rv;
	for(idx = 0;idx<size;idx++){
		rv.push_back(p.at(idx)/m);
	}
	return rv;
}

vector<float> convert(vector<int> vec){
	vector<float> rv;
	int size = vec.size();
	int idx = 0;
	for(idx = 0;idx<size;idx++){
		rv.push_back((float)vec.at(idx));
	}
	return rv;
}

void print_points(vector<vector<int> > points, map<vector<int>,int> mapping, bool with_map){
	int num_points = points.size();
	int idx = 0;
	if(with_map){
		for(idx = 0;idx<num_points;idx++){
			cout << "(" << points.at(idx).at(0) << "," << points.at(idx).at(1) << "," << points.at(idx).at(2) << ") -> " << mapping.at(points.at(idx)) << endl;
		}
	}else{
		for(idx = 0;idx<num_points;idx++){
			cout << "(" << points.at(idx).at(0) << "," << points.at(idx).at(1) << "," << points.at(idx).at(2) << ")" << endl;
		}
	}
}

void print_centroids(vector<vector<float> > centroids){
	int num_centroids = centroids.size();
	int idx = 0;
	for(idx = 0;idx<num_centroids;idx++){
		cout << "(" << centroids.at(idx).at(0) << "," << centroids.at(idx).at(1) << "," << centroids.at(idx).at(2) << ")" << endl;
	}
}

float objective_function(vector<vector<int> > points, map<vector<int>,int> mapping){
	float cost = 0;
	return cost;
}

void nearest_centroid(vector<vector<int> > points,vector<vector<float> > centroids,map<vector<int>,int> mapping){
	int num_points = points.size();
	int num_clusters = centroids.size();
	int point_idx = 0;
	int centroid_idx = 0;
	float min_dist = 0.0;
	int min_idx = -1;
	float dist = 0.0;
	for(point_idx = 0;point_idx<num_points;point_idx++){
		min_dist = 3000;
		for(centroid_idx = 0;centroid_idx<num_clusters;centroid_idx++){
			dist = distance(points.at(point_idx),centroids.at(centroid_idx));
			if(dist<min_dist){
				min_dist = dist;
				min_idx = centroid_idx;
			}
		}
		mapping.erase(points.at(point_idx));
		mapping.insert(pair<vector<int>,int>(points.at(point_idx),min_idx));
	}
	return;
}

void centroid_update(vector<vector<int> > points,vector<vector<float> > centroids,map<vector<int>,int> mapping){
	int num_features = 3;
	int num_points = points.size();
	int k = centroids.size();
	int idx = 0;
	int index = 0;
	centroids.clear();
	vector<float> vec;
	vector<int> num;
	vec.assign(num_features,0.0);
	num.assign(k,0);
	centroids.assign(k,vec);
	
	for(idx = 0;idx<num_points;idx++){
		index = mapping.at(points.at(idx));
		centroids[index] = add(centroids.at(index),points.at(idx));
		num[index]+=1;
	}

	for(idx = 0;idx<k;idx++){
		centroids[idx] = divide(centroids.at(idx),num.at(idx));
	}
	return;
}

void algo(vector<vector<int> > points, vector<vector<float> > centroids, map<vector<int>,int> mapping){
	int num_iterations = 100;
	while(num_iterations!=0){
		nearest_centroid(points,centroids,mapping);
		centroid_update(points,centroids,mapping);
		num_iterations-=1;
	}
}

void initialize(vector<vector<int> > points, vector<vector<float> > centroids, map<vector<int>,int> mapping){
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
            points.push_back(p1);
            mapping.insert(pair<vector<int>,int>(p1,0));
            p1.clear();
		}
        read_file.close();
	}

	// initializing the clusters
	while(k!=0){
		centroids.push_back(convert(points.at(k)));
		k-=1;
	}
	return;
}

int main(){
	vector<vector<int> > points;
	vector<vector<float> > centroids;
	map<vector<int>,int> mapping;
	initialize(points,centroids,mapping);
	algo(points,centroids,mapping);
	return 0;
}