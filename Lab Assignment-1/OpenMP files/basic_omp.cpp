#include <iostream>
#include <cstdlib>
#include <random>
#include <omp.h>
using namespace std;

int main(){
	int arr[64];
	//bool found = false;
	//int span = 8;
	int num_threads = 8;
    int tid;

	for(int i = 0;i<=63;i++){
        arr[i] = rand()%100;
	} 
    
    omp_set_num_threads(num_threads);
	#pragma omp parallel shared(arr) private(tid)
	{
        tid = omp_get_thread_num();
        int start = 8*tid;
        int end = 8*(tid)+7;
        int sum = 0;
        for(int i = start;i<=end;i++){
            sum+=arr[i];
        }
        cout << "Sum by thread " << tid << " is " << sum << endl;
	}
	return 0;
}