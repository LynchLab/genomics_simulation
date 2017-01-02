#!/bin/bash

name=$1

samtools mpileup -B ../sequences/seq*.sort.rmdup.bam -f ../sequences/$name -q 1 -Q 0 | gzip - > ../analysis_files/mpileup.txt
mapgd proview -H ../sequences/temp-header.txt -n ../sequences/name-file.txt -bs | gzip - > ../analysis_files/pro.txt.gz
mapgd allele -i ../analysis_files/pro.txt.gz -c 1 -g 5 | mapgd filter -q 0.001 -p 3 -g 5 | gzip - > ../analysis_files/mapgd_calls.txt.gz
zcat ../analysis_files/mapgd_calls.txt.gz | tail -n +6 | sed '$d' >  ../analysis_files/mapgd_calls-trim.csv
mapgd genotype -p ../analysis_files/pro.txt.gz -m ../analysis_files/mapgd_calls.txt.gz | mapgd relatedness > ../analysis_files/mapgd_relatedness.out
