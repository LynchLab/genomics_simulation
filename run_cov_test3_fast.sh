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

G=0

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

pedigree_sim -i states.txt -y r -s $S -g $G -olMm -k $G -v $V -e 0.005 --mu_d $mud --mu_a $mua --sigma_a $sa --sigma_d $sd --rho $rh #> $state_file

cd ../analysis_pipelines/

#cat $state_file | ../gsl_correl/call_relatedness_new > ../analysis_files/mapgd_relatedness.out 2> mapgd_relatedness.err

python ./python_utilities/add_salt.py ../analysis_files/ $e > ../analysis_files/plink.pheno
mapgd readphen -i ../analysis_files/plink.pheno > ../analysis_files/mapgd.phe

mle -i ../analysis_files/mapgd_relatedness.out -p ../analysis_files/mapgd.phe -h > ../analysis_files/output_i_${i}_f_${F}_ma_${mua}_md_${mud}_sa_${sa}_sd_${sd}_r_${rh}_e_${e}.txt 2> ../analysis_files/mle.err

file=output_i_${i}_f_${F}_ma_${mua}_md_${mud}_sa_${sa}_sd_${sd}_r_${rh}_e_${e}.txt

A=$(tail -n 1 ../analysis_files/$file | awk '{$1=$1};1' | cut -f 4 -d ' ')
D=$(tail -n 1 ../analysis_files/$file | awk '{$1=$1};1' | cut -f 6 -d ' ')
a=$(tail -n 1 ../analysis_files/$file | awk '{$1=$1};1' | cut -f 3 -d ' ')
d=$(tail -n 1 ../analysis_files/$file | awk '{$1=$1};1' | cut -f 5 -d ' ')
rho=$(tail -n 1 ../analysis_files/$file | awk '{$1=$1};1' | cut -f 7 -d ' ')

echo $A $D $a $d $rho

pwd

if [ $FIRST ]
then
	python python_utilities/join2.py ../analysis_files/true_var.hsq ../analysis_files/$file p >> ../analysis_files/estimates.txt
	FIRST=false
else
	python python_utilities/join2.py ../analysis_files/true_var.hsq ../analysis_files/$file >> ../analysis_files/estimates.txt
fi

cd ../

done
done
done
done
done
done
done
done
