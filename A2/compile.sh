#!/bin/bash
if [ "$1" == "s" ]; then
	clang sequential.c
else
	echo "Supply appropriate and correct arguments"
	exit
fi 
