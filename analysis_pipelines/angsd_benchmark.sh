#!/bin/bash
name=$1
cd ../sequences/
for x in 1 2 4
do
	export OMP_NUM_THREADS=$x
	export OMP_THREAD_LIMIT=$x
	echo -n "angds, doMaf, $x threads "
	(time (angsd -out ../analysis_files/angsd_calls.vcf -bam bam_list.txt -anc $name -GL 1 -doMaf 1 -doMajorMinor 5 -minMaf 0.001 -nThreads $x 2> /dev/null)) 2>&1 | head -2 | tail -n +2
done
cd ../analysis_pipelines/
