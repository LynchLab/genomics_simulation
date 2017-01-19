#!/bin/bash

for 
	mpirun -n $X mapgd proview -H ../sequences/temp-header.txt -n ../sequences/name-file.txt -bs > /dev/null
	mpirun -n $X mapgd allele -i ../analysis_files/pro.txt.gz -c 1 -g 5 > /dev/null
	mpirun -n $X mapgd relatedness -i ../analysis_files/mapgd_calls.gcf.gz > /dev/null
	mpirun -n $X mapgd linkage -p ../analysis_files/pro.txt.gz -m ../analysis_files/mapgd_calls.txt.gz  > /dev/null

