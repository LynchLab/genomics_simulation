rm -rf states.vcf
gunzip ../analysis_files/states.vcf.gz
plink --vcf ../analysis_files/states.vcf --make-bed --allow-extra-chr --out ../analysis_files/plink
gzip -f ../analysis_files/states.vcf

gcta64 --bfile ../analysis_files/plink --make-grm-gz --out ../analysis_files/gcta
gcta64 --bfile ../analysis_files/plink --make-grm --out ../analysis_files/gcta
