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
	cmake .. -DJUST_BUILD=1
	make -j$n_threads
else
	echo Incorrent number of threads
fi
