plink --vcf ../analysis_files/gatk_calls.vcf.gz --allow-extra-chr --r2 --ld-window-r2 0  --ld-window-kb 0.5 --genome --out ../analysis_files/plink
sed "s/^[ \t]*//" -i ../analysis_files/plink.ld
sed "s/ \+/ /g" -i ../analysis_files/plink.ld
