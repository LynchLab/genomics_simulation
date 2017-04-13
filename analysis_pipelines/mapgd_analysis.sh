#!/bin/bash

name=$1

samtools mpileup -B ../sequences/seq*.sort.rmdup.bam -f ../sequences/$name -q 5 -Q 7 | gzip - > ../analysis_files/mpileup.txt.gz
mapgd proview -H ../sequences/temp-header.txt -n ../sequences/name-file.txt -bs | gzip - > ../analysis_files/pro.txt.gz
echo "calling alleles"
mapgd allele -i ../analysis_files/pro.txt.gz -c 1 -g 20 -e 0.0001 | mapgd filter -q 0.001 -p 1 -g 30 -N 20 | gzip - > ../analysis_files/mapgd_calls.txt.gz
zcat ../analysis_files/mapgd_calls.txt.gz | tail -n +6 | sed '$d' >  ../analysis_files/mapgd_calls-trim.csv
echo "estimating genotypes"
mapgd genotype -p ../analysis_files/pro.txt.gz -m ../analysis_files/mapgd_calls.txt.gz | gzip - > ../analysis_files/genotype.gcf.gz
#echo "estimating ld"
#mapgd linkage -p ../analysis_files/pro.txt.gz -m ../analysis_files/mapgd_calls.txt.gz -D 500 > ../analysis_files/mapgd_linkage.out
#echo "estimating relatedness"
#cat ../analysis_files/genotype.gcf.gz | gunzip - | mapgd relatedness > ../analysis_files/mapgd_relatedness.out
#mapgd quant -r ../analysis_files/mapge_relatedness.out -p ../pedigree.txt 
