#!/bin/bash

source settings.sh

echo $SAMPLE

echo "simulating phenotypes"

cd ./analysis_pipelines

zcat ../sequences/states.txt.gz | python -u text_to_bin.py | python-2.7.9 -u states_to_vcf.py $SAMPLE ../sequences/name-file.txt > ../analysis_files/states.vcf
zcat ../sequences/states.txt.gz | python -u text_to_bin.py | gzip - > ../sequences/states_bin.txt.gz
zcat ../sequences/states_bin.txt.gz | ../fast_correl/call_relatedness $SAMPLE | gzip - > ../analysis_files/mapgd_relatedness.out.gz
zcat ../analysis_files/mapgd_relatedness.out.gz | head -3 | tail -1 | cut -f 3 -d '	' > ../analysis_files/inbred
zcat ../analysis_files/mapgd_relatedness.out.gz | head -$((SAMPLE+1)) | tail -$((SAMPLE-1)) | cut -f 5 -d '	' >> ../analysis_files/inbred

zcat ../sequences/states_bin.txt.gz | python-2.7.9 -u states_to_pheno_w_inbred.py ../sequences/name-file.txt $SNPS 500 ../analysis_files/inbred $A $D $E > ../analysis_files/plinki_$A_$D_$E.pheno

./gcta_gwas.sh 
./gvcblup_gwas.sh
#./plink_gwas.sh 
#./mapgd_gwas.sh
#./PRSice.sh

#Rscript make_figure_gwas_1.rscript	#Bias RMSE of allele frequenceis
