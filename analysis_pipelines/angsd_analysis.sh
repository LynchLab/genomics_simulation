#!/bin/bash
cd ../sequences/
angsd -out ../analysis_files/angsd_calls.vcf -bam bam_list.txt -anc reference.fa -GL 1 -doMaf 1 -doMajorMinor 5 -minMaf 0.001 -nThreads 10
angsd -GL 1 -out ../analysis_files/angsd_calls -doGlf 1 -bam bam_list.txt -anc referebce.fa
cd ../analysis_pipeline/
ibs
