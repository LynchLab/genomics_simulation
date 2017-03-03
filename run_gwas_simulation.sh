#!/bin/bash

POPULATION=1000
POP2=$((POPULATION/2))
SAMPLE=1000
TIME=$((2*POPULATION))
TIMEX=$((9*PULATION+POPULATION/2))
REF="reference.fa"
K=2000
SNPS=$(($K*10))  

echo "simulating population"
cd sequences
../population_simulation/pedigree_sim $POPULATION $TIME $SNPS 50 0.01 r t 1 -1 2> var | states_to_vcf.py states.txt
cat name-file.txt | cut -d '	' -f 1-$((SAMPLE+2)) > name-file2.txt
mv name-file2.txt name-file.txt
rm -rf pedigree.txt.gz
gzip pedigree.txt
cd ..

cd analysis_pipelines

#./mapgd_analysis_newton.sh $REF
./mapgd_gwas.sh $REF
./angsd_gwas.sh $REF
./gcta_gwas.sh 
#./gcta_analysis.sh $REF

Rscript make_figure_gwas_1.rscript	#Bias RMSE of allele frequenceis
