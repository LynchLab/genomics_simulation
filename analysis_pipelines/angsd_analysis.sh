#!/bin/bash
name=$1
cd ../sequences/
angsd -out ../analysis_files/angsd_calls.vcf -bam bam_list.txt -anc $name -GL 1 -doMaf 1 -doMajorMinor 5 -minMaf 0.001 -nThreads 10
angsd -GL 1 -out ../analysis_files/angsd_calls -doGlf 1 -bam bam_list.txt -anc $name
cd ../analysis_pipelines/
ibs
