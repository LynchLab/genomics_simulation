#!/bin/bash

GATK_PATH=~/src/gatk
PICARD_PATH=~/src/picard
root=`echo $1 | cut -d '.' -f 1`
ref=$1

for file in ../sequences/seq_*.sort.rmdup.bam;
do
	FILE="$FILE -I $file"
done


for x in 1 1,2 1,2,3,4
do
	export OMP_NUM_THREADS=$x
	export OMP_THREAD_LIMIT=$x
	echo -n "GATK, UnifiedGenotyper, $x threads " 
	(time taskset -c $x java -jar $GATK_PATH/GenomeAnalysisTK.jar -R ../sequences/$root.fa -T UnifiedGenotyper $FILE -o /dev/null -stand_call_conf 5 2> /dev/null ) 2>&1 | head -2 | tail -n +2
#	taskset -c $x java -jar $GATK_PATH/GenomeAnalysisTK.jar -R ../sequences/$root.fa -T UnifiedGenotyper $FILE -o /dev/null -stand_call_conf 5 
done
#mapgd linkage -p ../analysis_files/pro.txt.gz -m ../analysis_files/mapgd_calls.txt.gz  > ../analysis_files/mapgd_linkage.out
#cat genotype.gcf.gz | gunzip - | mapgd relatedness > ../analysis_files/mapgd_relatedness.out
#mapgd quant -r ../analysis_files/mapge_relatedness.out -p ../pedigree.txt 
