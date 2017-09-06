source settings.sh

cd sequences/
pedigree_sim -N 9408 -s 10000 -g 20001 -lbnGmdt -k 1000 -v 1000 -y r > ../sequences/states.bin
cd ../analysis_pipelines/
cat ../sequences/states.bin | call_relatedness ../sequences/name-file.txt > ../analysis_files/mapgd_relatedness.out
cat ../sequences/states.bin | python-2.7.9 -u states_to_vcf.py ../sequences/sequences/name-file.txt > ../analysis_files/states.vcf

plink --vcf ../analysis_files/states.vcf --make-bed --allow-extra-chr --out ../analysis_files/plink
gcta64 --bfile ../analysis_files/plink --make-grm --out ../analysis_files/gcta


