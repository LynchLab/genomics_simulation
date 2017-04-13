#!/bin/bash

POPULATION=400
POP2=$((POPULATION/2))
SAMPLE=400
TIME=$((2*POPULATION))
TIMEX=$((9*PULATION+POPULATION/2))
REF="reference.fa"
K=500
SNPS=$(($K*10))  

echo "simulating population"
cd sequences
../population_simulation/pedigree_sim $POPULATION $TIME $SNPS 50 0.003 r b 2> var | gzip - > states.txt.gz
#cat name-file.txt | cut -d '	' -f 1-$((SAMPLE+2)) > name-file2.txt
#mv name-file2.txt name-file.txt
#rm -rf pedigree.txt.gz
#gzip pedigree.txt
#cd ..

cd ../analysis_pipelines

#zcat ../sequences/states.txt.gz | python-2.7.9 -u states_to_vcf.py $POPULATION 0 > ../analysis_files/states.vcf
zcat ../sequences/states.txt.gz | ../fast_correl/fast_fast_correl_full $POPULATION 0 | gzip - > ../analysis_files/mapgd_relatedness.out.gz
zcat ../analysis_files/mapgd_relatedness.out.gz | head -3 | tail -1 | cut -f 3 -d '	' > inbred
zcat ../analysis_files/mapgd_relatedness.out.gz | head -$((POPULATION+1)) | tail -$((POPULATION-1)) | cut -f 5 -d '	' >> inbred
zcat ../sequences/states.txt.gz | python-2.7.9 -u states_to_pheno_w_inbred.py $POPULATION 0 $SNPS 500 inbred 1 -1 0.0001 > ../analysis_files/plink.pheno

./gcta_gwas.sh 
./plink_gwas.sh 
./mapgd_gwas.sh
./PRSice.sh

Rscript make_figure_gwas_1.rscript	#Bias RMSE of allele frequenceis
