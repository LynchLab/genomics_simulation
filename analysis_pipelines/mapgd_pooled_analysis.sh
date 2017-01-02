#!/bin/bash

name=$1

samtools mpileup -B ../sequences/pool.sort.bam -f ../sequences/$name -q 1 -Q 0 > ../analysis_files/mpileup.txt
cat ../analysis_files/mpileup.txt | mapgd proview -H ../sequences/temp-header.txt -bs > ../analysis_files/pro.txt
mapgd pool -i ../analysis_files/pro.txt -a 2 > ../analysis_files/mapgd_calls.txt
python get_frequencies_from_pol.py mapgd_calls.txt > ../analysis_files/mapgd_calls-trim.csv
