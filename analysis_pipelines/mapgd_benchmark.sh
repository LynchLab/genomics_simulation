#!/bin/bash

name=$1

for x in 1 2 4
do
	export OMP_NUM_THREADS=$x
	export OMP_THREAD_LIMIT=$x
	echo -n "mapgd, allele, $x threads "
	(time (samtools mpileup -B ../sequences/seq*.sort.rmdup.bam -f ../sequences/$name -q 1 -Q 0 2> /dev/null | mapgd proview -H ../sequences/temp-header.txt -n ../sequences/name-file.txt -bs | mapgd allele -c 1 -g 5 > /dev/null )) 2>&1 | head -2 | tail -n +2
#	(time mapgd allele -i ../analysis_files/pro.txt.gz -c 1 -g 5 > /dev/null) 2>&1 | head -2 | tail -n +2
#	echo -n "mapgd, relatedness, $x threads " 
#	(time mapgd relatedness -i ../analysis_files/genotype.gcf.gz > /dev/null) 2>&1 | head -2 | tail -n +2
done
#mapgd linkage -p ../analysis_files/pro.txt.gz -m ../analysis_files/mapgd_calls.txt.gz  > ../analysis_files/mapgd_linkage.out
#cat genotype.gcf.gz | gunzip - | mapgd relatedness > ../analysis_files/mapgd_relatedness.out
#mapgd quant -r ../analysis_files/mapge_relatedness.out -p ../pedigree.txt 
