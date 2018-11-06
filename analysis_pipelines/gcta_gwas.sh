M=0.15
N=10
include=include.txt
in=../analysis_files/gcta
out=../analysis_files/gcta_greml

gcta64 --grm-gz $in --pheno $1 --reml --out $out --thread-num 10 --grm-cutoff $M 
gcta64 --grm-gz $in --grm-cutoff $M --pca $N --out ../analysis_files/$N
gcta64 --grm-gz $in --pheno $1 --reml --out ${out}_${N}_egenvec --thread-num 10 --reml-pred-rand --qcovar ../analysis_files/$N.eigenvec
gcta64 --grm-gz $in --pheno $1 --HEreg --out ${out}_HE --thread-num 10 --grm-cutoff $M 
