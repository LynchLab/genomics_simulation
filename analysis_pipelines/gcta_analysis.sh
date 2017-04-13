gunzip ../sequences/pedigree.txt.gz
TARGET=`tail -1 ../sequences/pedigree.txt | cut -d '	' -f 1`
SIZE=`grep -c ^$TARGET ../sequences/pedigree.txt`
#python ../make_fam.py ../sequences/pedigree.txt | tail -n $SIZE |  sed 's/-9/-9.01/' > ../analysis_files/plink.fam
#cat ../analysis_files/plink.fam | cut -d '	' -f 1,2,6 > ../analysis_files/plink.pheno
gzip -f ../sequences/pedigree.txt

rm -rf states.vcf
gunzip ../analysis_files/states.vcf.gz
plink --vcf ../analysis_files/states.vcf --make-bed --allow-extra-chr --out ../analysis_files/plink
gzip -f ../analysis_files/states.vcf

gcta64 --bfile ../analysis_files/plink --make-grm-gz --out ../analysis_files/gcta
gcta64 --bfile ../analysis_files/plink --make-grm --out ../analysis_files/gcta
gcta64 --grm ../analysis_files/gcta --pheno ../analysis_files/plink.pheno --reml --out ../analysis_files/test --thread-num 10 --grm-cutoff 0.125 --reml-pred-rand

# --reml-pred-rand
#
#Predict the random effects by the BLUP (best linear unbiased prediction) method. This option is actually to predict the total genetic effect (called “breeding value” in animal genetics) of each individual attributed by the aggregative effect of the SNPs used to estimate the GRM. The total genetic effects of all the individuals will be saved in a plain ext file *.indi.blp.
#
#Output file format
#
#test.indi.blp (no header line; columns are family ID, individual ID, an intermediate variable, the total genetic effect,  another intermediate variable and the residual effect.
#
#If there are multiple GRMs fitted in the model, each GRM will insert additional two columns, , i.e. an intermediate variable and a total genetic effect, in front of the last two columns)
#
#01       0101    -0.012    -0.014    -0.010    -0.035
#
#02       0203    0.021     0.031    -0.027    -0.031
#
#03       0305    0.097     0.102    -0.026    -0.041
#
#…… 
#
#--blup-snp   test.indi.blp
#
#Calculate the BLUP solutions for the SNP effects (you have to specify the option --bfile to read the genotype data). This option takes the output of the option --reml-pred-rand as input  (*.indi.blp file) and transforms the BLUP solutions for individuals to the BLUP solutions for the SNPs, which can subsequently be used to predict the total genetic effect of individuals in an independent sample by PLINK --score option.
#
#Output file format
#
#test.snp.blp (columns are SNP ID, reference allele and BLUP of SNP effect; if there are multiple GRMs fitted in the model, each GRM will add an additional column to the file; the last column is for the residual effect)
#
#rs103645   A     0.00312    0.00451
#
#rs175292   G    -0.00021    0.00139
