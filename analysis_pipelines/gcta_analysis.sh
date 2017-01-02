plink --vcf ../analysis_files/gatk_calls.vcf --make-bed
python ../make_fam.py ../sequences/pedigree.txt | tail -n 500 > plink.fam
gcta64 --bfile plink --make-grm --out plink
cat plink.fam | cut -d '	' -f 1,2,5 > plink.pheno
gcta64 --grm plink --pheno plink.pheno --reml --out test --thread-num 10

