#!/bin/bash
source settings.sh

cd analysis_files/
#pedigree_sim -N 1200 -s 40000 -g 9001 -bnGmdt -k 1000 -v 1000 -e 0.008 > ../sequences/states.bin
#pedigree_sim -N 2700 -s 10000 -g 9001 -lbnGmdt -k 1000 -v 5000 -e 0.008 -Q 0 -S 0.003 > ../sequences/states.bin
#8998
#pedigree_sim -y r -T 8998 -N 2700 -s 10000 -g 9001 -lbnGmdt -k 1000 -v 20000 -e 0.008 -Q 0 -S 0.005 > ../sequences/states.bin

#		1	2	3	4	5	6	7	8	9	10
#1 mean_a==0	x		x		x		x		x
#2 mean_d==0		x	x			x	x			x
#4 dom					x	x	x	x	
#8 rec									x	x	x
#16 d=0  


pedigree_sim -y g -N 300 -s 100 -g 1001 -blnGmdt -k 1000 -v 3200 -e 0.008 -Q 4 -S 0.01 > ../sequences/states.bin
#pedigree_sim -T 8998 -N 10092 -s 10000 -g 9001 -bnGmdt -k 1000 -v 5000 -e 0.008 -Q 0 -S 0.003 > ../sequences/states.bin
#pedigree_sim -N 10092 -s 40000 -g 9001 -bnGmdt -k 1000 -v 1000 -e 0.008 > ../sequences/states.bin
cd ../analysis_pipelines/

#cat ../sequences/states.bin | call_relatedness ../analysis_files/name-file.txt > ../analysis_files/mapgd_relatedness.out
cat ../sequences/states.bin | ../gsl_correl/call_relatedness -m 0 > ../analysis_files/mapgd_relatedness_0.out
cat ../sequences/states.bin | ../gsl_correl/call_relatedness -m 1 > ../analysis_files/mapgd_relatedness_1.out
cp ../analysis_files/mapgd_relatedness_1.out ../analysis_files/mapgd_relatedness.out
cat ../sequences/states.bin | ../gsl_correl/call_relatedness -m 2 > ../analysis_files/mapgd_relatedness_2.out
cat ../sequences/states.bin | ../gsl_correl/call_relatedness -m 3 > ../analysis_files/mapgd_relatedness_3.out
cat ../sequences/states.bin | ../gsl_correl/call_relatedness -m 4 > ../analysis_files/mapgd_relatedness_4.out
cat ../sequences/states.bin | ../gsl_correl/call_relatedness -m 5 > ../analysis_files/mapgd_relatedness_5.out

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
