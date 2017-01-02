#!/bin/bash

cd analysis_pipelines

#./mapgd_analysis_newton.sh $REF
./bcftools_analysis.sh $REF
./angsd_analysis.sh $REF
./gatk_analysis.sh $REF

#./gcta_analysis.sh $REF

python get_frequencies.py ../sequences/states.txt ../sequences/polymorphisms.map > ../analysis_files/true_frequencies.csv
python get_frequencies_from_vcf.py ../analysis_files/gatk_calls.vcf > ../analysis_files/gatk_frequencies.csv
cat ../analysis_files/angsd_calls.vcf.mafs.gz | gunzip - | cut -d '	' -f 2,6 > ../analysis_files/angsd_frequencies.csv

Rscript make_figure_1.rscript
