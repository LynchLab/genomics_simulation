#!/bin/bash

name=$1

for x in 1 2 4
do
	export OMP_NUM_THREADS=$x
	export OMP_THREAD_LIMIT=$x
	echo -n "bcftools, call, threads $x "
	(time (bcftools mpileup -f ../sequences/$name ../sequences/seq_*.sort.rmdup.bam -q 5 -Q 5 2> /dev/null | bcftools call -m  --threads $x 2> /dev/null > /dev/null ) )  2>&1 | head -2 | tail -n +2
done
