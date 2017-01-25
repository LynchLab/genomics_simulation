#!/bin/bash
name=$1
cd ../sequences/
#angsd -out ../analysis_files/angsd_calls.vcf -bam bam_list.txt -anc $name -GL 1 -doSnpStat -doMaf 1 -doMajorMinor 1 -minMaf 0.001 -nThreads 10
angsd -HWE_pval_F 1 -out ../analysis_files/angsd_calls.vcf -bam bam_list.txt -minMapQ 7 -uniqueOnly 1 -doMajorMinor 5 -minMaf 0.001 -doMaf 1 -anc $name -GL 1 -nThreads 10
cd ../analysis_pipelines/
cat ../analysis_files/angsd_calls.vcf.hweF.gz | gunzip - | cut -d '	' -f 2,5,7 > ../analysis_files/angsd_frequencies.csv
ibs
