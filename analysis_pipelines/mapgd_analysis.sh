#!/bin/bash

name=$1
LD_DIST=$2

samtools mpileup -B ../sequences/seq*.sort.rmdup.bam -f ../sequences/$name -q 5 -Q 7 | gzip - > ../analysis_files/mpileup.txt.gz
mapgd proview -H ../sequences/temp-header.txt -n ../sequences/name-file.txt -bs | gzip - > ../analysis_files/pro.txt.gz
echo "calling alleles"
mapgd allele -i ../analysis_files/pro.txt.gz -c 1 -g 2 -e 0.0001 | mapgd filter -q 0.001 -p 10 -g 2 -N 1 | gzip - > ../analysis_files/mapgd_calls.txt.gz
#mapgd allele -i ../analysis_files/pro.txt.gz -c 1 -g 20 -e 0.0001 | mapgd filter -q 0.001 -p 1 -g 2 -N 1 | gzip - > ../analysis_files/mapgd_calls.txt.gz
zcat ../analysis_files/mapgd_calls.txt.gz | tail -n +6 | sed '$d' >  ../analysis_files/mapgd_calls-trim.csv
mapgd allele -i ../analysis_files/pro.txt.gz -c 1 -g 2 -e 0.0001 | mapgd filter -q 0.01 -p 10 -g 10 -N 1 | gzip - > ../analysis_files/mapgd_calls.txt.gz
echo "estimating genotypes"
mapgd genotype -p ../analysis_files/pro.txt.gz -m ../analysis_files/mapgd_calls.txt.gz | gzip - > ../analysis_files/genotype.gcf.gz
echo "estimating ld"
mapgd linkage -i ../analysis_files/genotype.gcf.gz -D $LD_DIST | gzip - > ../analysis_files/mapgd_linkage.out.gz
zcat ../analysis_files/mapgd_linkage.out.gz | tail -n +6 | sed '$d' >  ../analysis_files/mapgd_linkage-trim.csv
echo "estimating relatedness"
cat ../analysis_files/genotype.gcf.gz | gunzip - | mapgd relatedness > ../analysis_files/mapgd_relatedness.out
echo "estimating quantitive componenets"

#mapgd quant -r ../analysis_files/mapge_relatedness.out -p ../pedigree.txt 
