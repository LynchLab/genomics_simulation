#!/bin/bash

NUMBER=50
TIME=$((NUMBER))
REF="reference.fa"
COV=3
SIZE=50000 #bwa craps out at snp densities gt ~1:150
SNPS=500

echo "simulating reference"
python reference_simulation/mutation_simulation2.py -l 24 -S $SIZE > ./sequences/reference_mutations.txt
python reference_simulation/mutation_simulation.py -m ./sequences/reference_mutations.txt -s ./reference_simulation/seed.fa > ./sequences/reference.fa

echo "simulating population"
cd sequences
../population_simulation/pedigree_sim $NUMBER $TIME $SNPS 1 0 > states.txt 2> var
cd ..

echo "making individual genomes"
cd variant_simulation
bash state_to_fasta.sh ../sequences/reference.fa ../sequences/states.txt
cd ..
echo "sequencing..."

for x in $(seq -f "%03g" 0 1 $((NUMBER-1)) )
do 
	NAME=seq_$x

	echo $NAME

	./sequencing_simulation/art_illumina -ss HS25 -sam -i ./sequences/$NAME.0.fa -p -l 150 -f $COV -m 200 -s 10 -o temp.0 > /dev/null
	./sequencing_simulation/art_illumina -ss HS25 -sam -i ./sequences/$NAME.1.fa -p -l 150 -f $COV -m 200 -s 10 -o temp.1 > /dev/null

	cat temp.01.fq > ./sequences/temp.1.fq
	cat temp.11.fq >> ./sequences/temp.1.fq
	cat temp.02.fq > ./sequences/temp.2.fq
	cat temp.12.fq >> ./sequences/temp.2.fq

	bash ./alignment/run_alignment.sh sequences/$REF ./sequences/temp ./sequences/$NAME $(($NUMBER*$TIME+10#$x)) > /dev/null 2> /dev/null
#	bash ./alignment/run_alignment.sh sequences/$REF ./sequences/temp ./sequences/$NAME $(($NUMBER*$TIME+10#$x)) 

	rm temp.01.fq
	rm temp.11.fq
	rm temp.02.fq
	rm temp.12.fq

	rm temp.01.aln
	rm temp.11.aln
	rm temp.02.aln
	rm temp.12.aln

	rm temp.0.sam
	rm temp.1.sam

#	rm seqeuences/temp.1.fq
#	rm seqeuences/temp.2.fq

	cd sequences
	samtools index $NAME.sort.rmdup.bam
	cd ..
done

ls sequences/ | grep seq.*.sort.rmdup.bam$ > sequences/bam_list.txt

cd analysis_pipelines

#./mapgd_analysis_newton.sh $REF
./mapgd_analysis.sh $REF
exit
./bcftools_analysis.sh $REF
./angsd_analysis.sh $REF
./gatk_analysis.sh $REF

#./gcta_analysis.sh $REF

python get_frequencies.py ../sequences/states.txt ../sequences/polymorphisms.map > ../analysis_files/true_frequencies.csv
python get_frequencies_from_vcf.py ../analysis_files/gatk_calls.vcf > ../analysis_files/gatk_frequencies.csv
cat ../analysis_files/angsd_calls.vcf.mafs.gz | gunzip - | cut -d '	' -f 2,6 > ../analysis_files/angsd_frequencies.csv

Rscript make_figure_1.rscript
