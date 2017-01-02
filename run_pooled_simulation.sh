#!/bin/bash

NUMBER=100
TIME=$((NUMBER))
REF="reference.fa"
COV=0.5

echo "simulating reference"
python reference_simulation/mutation_simulation2.py -l 24 -S 5000000 > sequences/reference_mutations.txt
python reference_simulation/mutation_simulation.py -m sequences/reference_mutations.txt -s reference_simulation/seed.fa > sequences/reference.fa

echo "simulating population"
cd sequences
../population_simulation/pedigree_sim $NUMBER $TIME 6250 1 1 > states.txt 2> var
cd ..

echo "making individual genomes"
cd variant_simulation
bash state_to_fasta.sh ../sequences/reference.fa ../sequences/states.txt
cd ..
echo "sequencing..."
	
rm -rf ./sequences/temp.1.fq
rm -rf ./sequences/temp.2.fq

for x in $(seq -f "%03g" 0 1 $((NUMBER-1)) )
do 
	NAME=seq_$x

	echo $NAME

	./sequencing_simulation/art_illumina -ss HS25 -sam -i ./sequences/$NAME.0.fa -p -l 150 -f $COV -m 200 -s 10 -o temp.0 > /dev/null
	./sequencing_simulation/art_illumina -ss HS25 -sam -i ./sequences/$NAME.1.fa -p -l 150 -f $COV -m 200 -s 10 -o temp.1 > /dev/null

	rm *.aln
	rm *.sam

	cat temp.01.fq >> ./sequences/temp.1.fq
	cat temp.11.fq >> ./sequences/temp.1.fq
	cat temp.02.fq >> ./sequences/temp.2.fq
	cat temp.12.fq >> ./sequences/temp.2.fq
done 

bash ./alignment/run_alignment.sh sequences/$REF ./sequences/temp ./sequences/pool pool > /dev/null 2> /dev/null

cd sequences
	samtools index pool.sort.bam
	samtools index pool.sort.rmdup.bam
cd ..

ls sequences/ | grep seq.*.sort.rmdup.bam$ > sequences/bam_list.txt

cd analysis_pipelines

./mapgd_pooled_analysis.sh $REF
./breseq_pooled_analysis.sh $REF

python get_frequencies.py ../sequences/states.txt ../sequences/polymorphisms.map > ../analysis_files/true_frequencies.csv

Rscript make_figure_2.rscript
