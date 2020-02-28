#!/bin/bash
if [ "$1" == "s" ]; then
	clang sequential.c
elif [ "$1" == "b" ]; then
	mpicc mpi_p2p_b.c
elif [ "$1" == "n" ]; then
	mpicc mpi_p2p_n.c
elif [ "$1" == "c" ]; then
	mpicc mpi_c.c
else
	echo "Supply appropriate and correct arguments"
	exit
fi 
