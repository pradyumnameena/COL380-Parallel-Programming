#!/bin/bash
if [ "$1" == "s" ]; then
	clang++ sequential.cpp
elif [ "$1" == "p" ]; then
	clang++ -Xpreprocessor -pthread pthread.cpp -lpthread
elif [ "$1" == "o" ]; then
	clang++ -Xpreprocessor -fopenmp openmp.cpp -lomp
else
	echo "Supply appropriate and correct arguments"
	exit
fi 
