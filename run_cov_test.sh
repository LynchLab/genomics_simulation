#!/bin/bash
source settings.sh

rm -rf ./analysis_files/estimates.txt


#		1	2	3	4	5	6	7	8	9	10
#1 mean_a==0	x		x		x		x		x
#2 mean_d==0		x	x			x	x			x
#4 dom					x	x	x	x	
#8 rec									x	x	x
#16 d=0 r


for i in {1..1000}
do

cd analysis_files/

#pedigree_sim -y r -N 588 -s 500 -g 5001 -blnGmdt -k 1000 -v 16000 -e 0.03 --mu_d 0 --mu_a 0 --sigma_a 1 --sigma_d 1 --rho 0 > ../sequences/states.bin
pedigree_sim -y r -N 588 -s 500 -g 5001 -blnGmdt -k 1000 -v 16000 -e 0.005 --mu_d 0.0 --mu_a 0.0 --sigma_a 1.0 --sigma_d 0.5 --rho 0.0 > ../sequences/states.bin

cd ../analysis_pipelines/

cat ../sequences/states.bin | ../gsl_correl/call_relatedness -m 2 > ../analysis_files/mapgd_relatedness.out
#cat ../sequences/states.bin | ../gsl_correl/call_relatedness -m 1 > ../analysis_files/mapgd_relatedness_1.out
#cat ../sequences/states.bin | ../gsl_correl/call_relatedness -m 2 > ../analysis_files/mapgd_relatedness_2.out
#cat ../sequences/states.bin | ../gsl_correl/call_relatedness -m 3 > ../analysis_files/mapgd_relatedness_3.out
#cat ../sequences/states.bin | ../gsl_correl/call_relatedness -m 4 > ../analysis_files/mapgd_relatedness_4.out
#cat ../sequences/states.bin | ../gsl_correl/call_relatedness -m 5 > ../analysis_files/mapgd_relatedness_5.out

#cp ../analysis_files/mapgd_relatedness_1.out ../analysis_files/mapgd_relatedness.out
#cp ../analysis_files/mapgd_relatedness_2.out ../analysis_files/mapgd_relatedness.out

mapgd writevcf2 -s ../sequences/states.bin -n ../analysis_files/name-file.txt -o ../analysis_files/states.vcf

plink --vcf ../analysis_files/states.vcf --make-bed --maf 0.00001 --out ../analysis_files/plink
plink --vcf ../analysis_files/states.vcf --genome --maf 0.00001 --out ../analysis_files/plink
gcta64 --bfile ../analysis_files/plink --make-grm-d-gz --out ../analysis_files/gcta
gcta64 --bfile ../analysis_files/plink --make-grm-gz --out ../analysis_files/gcta

#gcta64 --bfile ../analysis_files/plink --make-grm-part 3 1 --out ../analysis_files/gcta
#gcta64 --bfile ../analysis_files/plink --make-grm-part 3 2 --out ../analysis_files/gcta
#gcta64 --bfile ../analysis_files/plink --make-grm-part 3 3 --out ../analysis_files/gcta

cd make_figures/
Rscript make_csv.rscript
cd ..
python ./python_utilities/add_salt.py ../analysis_files/ 1 > ../analysis_files/plink.pheno
exit 0

python ./python_utilities/add_salt.py ../analysis_files/ 0.1 > ../analysis_files/plink.pheno
./gcta_gwas.sh ../analysis_files/plink.pheno
python python_utilities/ln_maximize2.py ../analysis_files/ > ../analysis_files/mapgd.hsq
python python_utilities/join.py ../analysis_files/true_var.hsq ../analysis_files/test.hsq ../analysis_files/mapgd.hsq  >> ../analysis_files/estimates.txt

python ./python_utilities/add_salt.py ../analysis_files/ 0.25 > ../analysis_files/plink.pheno
./gcta_gwas.sh ../analysis_files/plink.pheno
python python_utilities/ln_maximize2.py ../analysis_files/ > ../analysis_files/mapgd.hsq
python python_utilities/join.py ../analysis_files/true_var.hsq ../analysis_files/test.hsq ../analysis_files/mapgd.hsq  >> ../analysis_files/estimates.txt

python ./python_utilities/add_salt.py ../analysis_files/ 0.50 > ../analysis_files/plink.pheno
./gcta_gwas.sh ../analysis_files/plink.pheno
python python_utilities/ln_maximize2.py ../analysis_files/ > ../analysis_files/mapgd.hsq
python python_utilities/join.py ../analysis_files/true_var.hsq ../analysis_files/test.hsq ../analysis_files/mapgd.hsq  >> ../analysis_files/estimates.txt

python ./python_utilities/add_salt.py ../analysis_files/ 0.75 > ../analysis_files/plink.pheno
./gcta_gwas.sh ../analysis_files/plink.pheno
python python_utilities/ln_maximize2.py ../analysis_files/ > ../analysis_files/mapgd.hsq
python python_utilities/join.py ../analysis_files/true_var.hsq ../analysis_files/test.hsq ../analysis_files/mapgd.hsq  >> ../analysis_files/estimates.txt

python ./python_utilities/add_salt.py ../analysis_files/ 0.9 > ../analysis_files/plink.pheno
./gcta_gwas.sh ../analysis_files/plink.pheno
python python_utilities/ln_maximize2.py ../analysis_files/ > ../analysis_files/mapgd.hsq
python python_utilities/join.py ../analysis_files/true_var.hsq ../analysis_files/test.hsq ../analysis_files/mapgd.hsq  >> ../analysis_files/estimates.txt

cd ../

done
#cd make_figures/
#Rscript make_figure_cov\(beta\,beta\).rscript
#Rscript make_h2_image.rscript
#cd ..
