#!/bin/bash

name=$1

mapgd allele -i ../analysis_files/pro.txt.gz -c 1 -g 5 -n | mapgd filter -q 0.001 -p 3 -g 5 | gzip - > ../analysis_files/mapgd_calls_newton.txt.gz
zcat ../analysis_files/mapgd_calls_newton.txt.gz | tail -n +6 | sed '$d' >  ../analysis_files/mapgd_calls_newton-trim.csv
