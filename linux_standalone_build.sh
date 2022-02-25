#!/bin/bash

echo Enter the number of threads you want to use to build the project
echo To use a default \(4\) value just press Enter.

read n_threads

if [ -z "$n_threads" ] || [ -z "${n_threads}" ];
then
	n_threads=4
fi

if [[ $n_threads =~ ^-?[0-9]+$ ]] && [ $n_threads -gt 0 ];
then
	echo Start building using $n_threads threads
	mkdir build
	cd build
	cmake .. -DJUST_BUILD=1 -D CMAKE_C_COMPILER=gcc-7 -D CMAKE_CXX_COMPILER=g++-7
	make -j$n_threads CC=gcc-7 CPP=g++-7 CXX=g++-7 LD=g++-7
else
	echo Incorrent number of threads
fi
