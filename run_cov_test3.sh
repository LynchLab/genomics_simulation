#!/bin/bash
source settings.sh

rm -rf ./analysis_files/estimates.txt

#@NS:15427	GS:15627	BS:1418352456

n=1200
#n=12
N=1200

S=10000
V=320000

#N=15472
#S=15600
#V=499200

G=1

state_file=../sequences/states.txt

FIRST=true

rm ./analysis_files/estimates.txt

for i in {1..10}
do
for mua in 0.0 -0.25 0.25
do
for mud in 0.0 -0.25 0.25
do
for rh in 0.0 -0.25 0.25  
do
for sd in 1.0 0.0625 0.25 
do
for sa in 1.0 0.0625 0.25  
do
for e in 0.5 0.75 0.25
do
for F in 0 8 32 16 
do
echo $ma $md $sa $sd $rh
cd analysis_files/

if [ $F -eq 0 ]  
then
	pedigree_sim -i ms_state.bin -y r -N $n -1 $N -s $S -g $G -plunMmGt -k $G -v $V -e 0.005 --mu_d $mud --mu_a $mua --sigma_a $sa --sigma_d $sd --rho $rh > $state_file
else 
	pedigree_sim -i ms_state.bin -y g -f $F -N $n -1 $N -s $S -g $G -plunMmGt -k $G -v $V -e 0.005 --mu_d $mud --mu_a $mua --sigma_a $sa --sigma_d $sd --rho $rh > $state_file
fi

cd ../analysis_pipelines/

cat $state_file | ../gsl_correl/call_relatedness_new > ../analysis_files/mapgd_relatedness.out 2> mapgd_relatedness.err

python ./python_utilities/add_salt.py ../analysis_files/ $e > ../analysis_files/plink.pheno
mapgd readphen -i ../analysis_files/plink.pheno > ../analysis_files/mapgd.phe

echo "Running stupid crap..."

#python ./python_utilities/make_true_A.py ../analysis_files/cov.txt > ../analysis_files/true_A.csv
#python ./python_utilities/make_true_D.py ../analysis_files/cov.txt > ../analysis_files/true_D.csv
#python ./python_utilities/make_true_G.py ../analysis_files/cov.txt > ../analysis_files/true_G.csv

#python ./python_utilities/make_true_A2.py ../analysis_files/cov2.txt > ../analysis_files/true_A2.csv
#python ./python_utilities/make_true_D2.py ../analysis_files/cov2.txt > ../analysis_files/true_D2.csv
#python ./python_utilities/make_true_G2.py ../analysis_files/cov2.txt > ../analysis_files/true_G2.csv

mapgd writevcf2 -s $state_file -n ../analysis_files/name-file.txt -o ../analysis_files/states.vcf

plink --vcf ../analysis_files/states.vcf --make-bed --maf 0.00001 --out ../analysis_files/plink
plink --vcf ../analysis_files/states.vcf --genome --maf 0.00001 --out ../analysis_files/plink
gcta64 --bfile ../analysis_files/plink --make-grm-d-gz  --out ../analysis_files/gcta
gcta64 --bfile ../analysis_files/plink --make-grm-gz  --out ../analysis_files/gcta

#mapgd readbed -n $N ../analysis_files/plink.bed -b > ../analysis_files/states.bin
#cat ../analysis_files/states.bin | ../gsl_correl/call_relatedness_new > ../analysis_files/mapgd_relatedness_2.out

mle -i ../analysis_files/mapgd_relatedness.out -p ../analysis_files/mapgd.phe  > ../analysis_files/output_i_${i}_f_${F}_ma_${mua}_md_${mud}_sa_${sa}_sd_${sd}_r_${rh}_e_${e}.txt 2> ../analysis_files/mle.err

file=output_i_${i}_f_${F}_ma_${mua}_md_${mud}_sa_${sa}_sd_${sd}_r_${rh}_e_${e}.txt

A=$(tail -n 1 ../analysis_files/$file | awk '{$1=$1};1' | cut -f 4 -d ' ')
D=$(tail -n 1 ../analysis_files/$file | awk '{$1=$1};1' | cut -f 6 -d ' ')
a=$(tail -n 1 ../analysis_files/$file | awk '{$1=$1};1' | cut -f 3 -d ' ')
d=$(tail -n 1 ../analysis_files/$file | awk '{$1=$1};1' | cut -f 5 -d ' ')
rho=$(tail -n 1 ../analysis_files/$file | awk '{$1=$1};1' | cut -f 7 -d ' ')

#A=$mua
#D=$mud
#a=$sa
#d=$sd
#rho=$rh

echo $A $D $a $d $rho

cat $state_file | ../gsl_correl/call_relatedness_add -A $A -D  $D -a $a -d $d -e $rho > ../analysis_files/mapgd_relatedness_add.out 2> ../analysis_files/mapgd.hsq
#cat $state_file | ../gsl_correl/call_relatedness_add -e ${rh} > ../analysis_files/mapgd_relatedness_add.out 2> ../analysis_files/mapgd.hsq

./gcta_gwas.sh ../analysis_files/plink.pheno

if [ $FIRST ]
then
	python python_utilities/join.py ../analysis_files/true_var.hsq ../analysis_files/gcta_greml_10_egenvec.hsq ../analysis_files/gcta_greml.hsq ../analysis_files/mapgd.hsq  ../analysis_files/$file p >> ../analysis_files/estimates.txt
	FIRST=false
else
	python python_utilities/join.py ../analysis_files/true_var.hsq ../analysis_files/gcta_greml_10_egenvec.hsq ../analysis_files/gcta_greml.hsq ../analysis_files/mapgd.hsq  ../analysis_files/$file >> ../analysis_files/estimates.txt
fi
cd make_figures
Rscript make_figure_cov\(beta\,beta\).rscript 

dir=i:${i}_f:${F}_ma:${mua}_md:${mud}_sa:${sa}_sd:${sd}_rh:${rh}_e:${e}

mkdir $dir

mv *.png $dir/

cd ../

cd ../

exit

done
done
done
done
done
done
done
done
