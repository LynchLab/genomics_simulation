LD_DIST=`echo print \$\(\( $1/1000. \)\) | ksh`
plink --vcf ../analysis_files/gatk_calls.vcf.gz --allow-extra-chr --r2 --ld-window-r2 0  --ld-window-kb $LD_DIST --genome --out ../analysis_files/plink
sed "s/^[ \t]*//" -i ../analysis_files/plink.ld
sed "s/ \+/ /g" -i ../analysis_files/plink.ld
