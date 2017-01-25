#!/bin/bash

name=$1

bcftools mpileup -f ../sequences/$name ../sequences/seq_*.sort.rmdup.bam -q 5 -Q 7 | bcftools call -m > ../analysis_files/bcftools_calls.vcf
python get_frequencies_from_vcf.py ../analysis_files/bcftools_calls.vcf > ../analysis_files/bcf_frequencies.csv
