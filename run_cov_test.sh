#!/bin/bash
source settings.sh

cd analysis_files/
#pedigree_sim -N 1200 -s 40000 -g 9001 -bnGmdt -k 1000 -v 1000 -e 0.008 > ../sequences/states.bin
#pedigree_sim -N 2700 -s 10000 -g 9001 -lbnGmdt -k 1000 -v 5000 -e 0.008 -Q 0 -S 0.003 > ../sequences/states.bin
pedigree_sim -y a -T 8998 -N 2700 -s 10000 -g 9001 -lbnGmdt -k 1000 -v 5000 -e 0.008 -Q 0 -S 0.003 > ../sequences/states.bin
#pedigree_sim -N 9408 -s 40000 -g 9001 -bnGmdt -k 1000 -v 1000 -e 0.008 > ../sequences/states.bin
cd ../analysis_pipelines/
cat ../sequences/states.bin | call_relatedness ../analysis_files/name-file.txt > ../analysis_files/mapgd_relatedness.out
#cat ../sequences/states.bin | get_kappa_f ../analysis_files/name-file.txt > ../analysis_files/const.txt
mapgd writevcf2 -s ../sequences/states.bin -n ../analysis_files/name-file.txt -o ../analysis_files/states.vcf

plink --vcf ../analysis_files/states.vcf --make-bed --maf 0.00001 --out ../analysis_files/plink
plink --vcf ../analysis_files/states.vcf --genome --maf 0.00001 --out ../analysis_files/plink
gcta64 --bfile ../analysis_files/plink --make-grm-d-gz --out ../analysis_files/gcta
gcta64 --bfile ../analysis_files/plink --make-grm-gz --out ../analysis_files/gcta

#gcta64 --bfile ../analysis_files/plink --make-grm-part 3 1 --out ../analysis_files/gcta
#gcta64 --bfile ../analysis_files/plink --make-grm-part 3 2 --out ../analysis_files/gcta
#gcta64 --bfile ../analysis_files/plink --make-grm-part 3 3 --out ../analysis_files/gcta

cd make_figures/
Rscript make_figure_cov\(beta\,beta\).rscript
cd ..
