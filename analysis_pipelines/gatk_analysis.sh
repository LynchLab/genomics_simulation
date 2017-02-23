GATK_PATH=~/src/gatk
PICARD_PATH=~/src/picard

HET=0.01

ref=$1
root=`echo $1 | cut -d '.' -f 1`
name_pre=`echo $name | cut -f 1 -d '.'`
rm -f ../sequences/$root.dict
java -jar $PICARD_PATH/picard.jar CreateSequenceDictionary R=../sequences/$root.fa O=../sequences/$root.dict

for file in ../sequences/seq_*.sort.rmdup.bam;
do
	FILE="$FILE -I $file"
done

java -jar $GATK_PATH/GenomeAnalysisTK.jar -R ../sequences/$root.fa -T UnifiedGenotyper $FILE -o ../analysis_files/gatk_calls.vcf -stand_call_conf 5

python get_frequencies_from_vcf.py ../analysis_files/gatk_calls.vcf > ../analysis_files/gatk_frequencies.csv
gzip ../analysis_files/gatk_calls.vcf
