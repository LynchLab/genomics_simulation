REF=$1
VAR=$2
SAMPLE=$3

bgzip ../analysis_files/$VAR
tabix -f ../analysis_files/$VAR.gz
vg construct -r ../sequences/$REF -v ../analysis_files/$VAR.gz > ../analysis_files/graph.vg
vg index -x ../analysis_files/graph.xg -g ../analysis_files/graph.gcsa -k 11 ../analysis_files/graph.vg

for x in $(seq -f "%03g" 0 1 $((SAMPLE-1)) )
do
	echo "mapping $x"
	vg map -x ../analysis_files/graph.xg -g ../analysis_files/graph.gcsa -f ../sequences/seq_$x.0.fq -f ../sequences/seq_$x.1.fq | vg sift -q 7 - | ../sequences/seq_$x.gam
	vg surject -p "gb|BK006935.2|" -b -x ../analysis_files/graph.xg ../sequences/seq_$x.gam | samtools sort - | samtools rmdup - - > ../sequences/seq_$x-vg.bam
	vg index -N ../sequences/seq_$x.gam
done

#vg genotype -v $x.vg $x.gam.index/ > ../analysis_files/vg-calls.vcf
