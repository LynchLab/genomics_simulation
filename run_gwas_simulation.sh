#!/bin/bash

source settings.sh

echo $SAMPLE

echo "simulating phenotypes"

cd ./analysis_pipelines

mapgd writevcf2 -s ../sequences/states.bin -n ../analysis_files/name-file.txt -o ../analysis_files/states.vcf
zcat ../sequences/states.bin | ../fast_correl/call_relatedness $SAMPLE | gzip - > ../analysis_files/mapgd_relatedness.out.gz

cat ../analysis_files/t_final.txt | head -n -1| tail -n +3| awk -F '	' -v OFS='	' '{print $1, $1, $2}' > ../analysis_files/plinki_$A_$D_$E.pheno

./gcta_gwas.sh ../analysis_files/plinki_$A_$D_$E.pheno 
./gvcblup_gwas.sh
#./plink_gwas.sh 
#./mapgd_gwas.sh
#./PRSice.sh

#Rscript make_figure_gwas_1.rscript	#Bias RMSE of allele frequenceis
