#!/bin/bash

source settings.sh

echo "simulating population"
cd sequences
#../population_simulation/non_pedigree_sim $POPULATION $SNPS t 2> var | cut -d '	' -f 1-$((2*SAMPLE+1)) > states.txt
#../population_simulation/pedigree_sim $POPULATION $TIME $SNPS 0.002 0.002 s b 1 0.5 2> var | cut -d '	' -f 1-$((2*SAMPLE+1)) > states.txt
../population_simulation/pedigree_sim $POPULATION $TIME $SNPS 50 0.01 r t 2> var | cut -d '	' -f 1$SAMPLE_CHRM > states.txt
head -2 name-file.txt > temp-file.txt
head -3 name-file.txt | tail -n 1 | cut -d '	' -f 1$SAMPLE_NAME >> temp-file.txt
tail -n 1 name-file.txt >> temp-file.txt
mv temp-file.txt name-file.txt
rm -rf pedigree.txt.gz
	gzip pedigree.txt
cd ..

case $REFTYPE in
  [S]   )
		echo "simulating reference"
		python reference_simulation/mutation_simulation2.py -l 24 -S $SIZE > ./sequences/reference_mutations.txt
		python reference_simulation/mutation_simulation.py -m ./sequences/reference_mutations.txt -s ./reference_simulation/seed.fa > ./sequences/$REF
		;;
  [Y]   )
		echo "using yeast chromosome IV"
		cp ./real_genomes/S288C_Chromosome\ IV.fsa ./sequences/$REF
		;;
  [D]   )
                echo "using dmel3R "
                cp ./real_genomes/dmel-3R.fa ./sequences/$REF
                ;;
esac

echo "making individual genomes"
cd variant_simulation
bash state_to_fasta.sh ../sequences/reference.fa ../sequences/states.txt
cd ..
echo "sequencing..."

rm -rf ./sequences/states.txt.gz
gzip ./sequences/states.txt

for x in $(seq -f "%03g" 0 1 $((SAMPLE-1)) )
do 
	NAME=seq_$x

	echo $NAME

	gunzip ./sequences/$NAME.0.fa.gz
	gunzip ./sequences/$NAME.1.fa.gz

	./sequencing_simulation/art_illumina -qs -15 -qs2 -15 -ss HS25 -sam -i ./sequences/$NAME.0.fa -p -l 150 -f $COV -m 200 -s 10 -o temp.0 > /dev/null
	./sequencing_simulation/art_illumina -qs -15 -qs2 -15 -ss HS25 -sam -i ./sequences/$NAME.1.fa -p -l 150 -f $COV -m 200 -s 10 -o temp.1 > /dev/null

#	./sequencing_simulation/art_illumina -ss HS25 -sam -i ./sequences/$NAME.0.fa -p -l 150 -f $COV -m 200 -s 10 -o temp.0 > /dev/null
#	./sequencing_simulation/art_illumina -ss HS25 -sam -i ./sequences/$NAME.1.fa -p -l 150 -f $COV -m 200 -s 10 -o temp.1 > /dev/null

	rm ./sequences/$NAME.0.fa
	rm ./sequences/$NAME.1.fa

	cat temp.01.fq > ./sequences/temp.1.fq
	cat temp.11.fq >> ./sequences/temp.1.fq
	cat temp.02.fq > ./sequences/temp.2.fq
	cat temp.12.fq >> ./sequences/temp.2.fq

	bash ./alignment/run_mem_alignment.sh sequences/$REF ./sequences/temp ./sequences/$NAME $(($SAMPLE*$TIME+10#$x)) > /dev/null 2> /dev/null
#	bash ./alignment/run_alignment.sh sequences/$REF ./sequences/temp ./sequences/$NAME $(($SAMPLE*$TIME+10#$x)) > /dev/null 2> /dev/null
#	bash ./alignment/run_alignment.sh sequences/$REF ./sequences/temp ./sequences/$NAME $(($SAMPLE*$TIME+10#$x)) 

#	mv ./sequences/temp.1.fq ./sequences/$NAME.0.fq
#	mv ./sequences/temp.2.fq ./sequences/$NAME.1.fq

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

	cd sequences
	samtools index $NAME.sort.rmdup.bam
	cd ..
done

ls sequences/ | grep seq.*.sort.rmdup.bam$ > sequences/bam_list.txt

rm seqeuences/temp.1.fq
rm seqeuences/temp.2.fq

cd analysis_pipelines

#./mapgd_analysis_newton.sh $REF
./mapgd_analysis.sh $REF $LD_DIST
./bcftools_analysis.sh $REF
./angsd_analysis.sh $REF
./gatk_analysis.sh $REF
./plink_analysis.sh $LD_DIST
./gcta_analysis.sh 

gunzip ../sequences/states.txt
python get_frequencies.py ../sequences/states.txt ../sequences/polymorphisms.map > ../analysis_files/true_frequencies.csv
python get_ld.py ../sequences/states.txt ../sequences/polymorphisms.map $LD_DIST > ../analysis_files/true_ld.csv
gzip ../sequences/states.txt

Rscript Ackerman2017/make_figure_1a.rscript	#Bias RMSE of allele frequenceis
Rscript Ackerman2017/make_figure_1c.rscript	#Bias RMSE of inbreeding
Rscript Ackerman2017/make_figure_1d.rscript	#Bias RMSE of LD.
Rscript Ackerman2017/make_figure_1e.rscript	#ROC.

#./mapgd_benchmark.sh $REF > ../analysis_files/benchmark.csv
#./bcftools_benchmark.sh $REF >> ../analysis_files/benchmark.csv
#./angsd_benchmark.sh $REF >> ../analysis_files/benchmark.csv
#./gatk_benchmark.sh $REF >> ../analysis_files/benchmark.csv
