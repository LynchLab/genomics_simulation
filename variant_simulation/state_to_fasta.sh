ref=$1
states=$2
SIZE="`head -1 $states | grep '	' -o | wc -l`"
POLY=../sequences/polymorphisms.map

python mutation_simulation2.py -n "`wc -l $states | cut -d ' ' -f 1`" -l "`head -1 $ref | cut -f 2 -d ':' `" -o True > temp
python mutation_sort.py temp > $POLY

rm temp

python state_to_fasta.py -N $((SIZE/2)) -s $states -m $POLY

for N in `seq 0 2 $((SIZE-2))` 
do
	name=`printf %03d $((N/2))`
	echo $name
	python mutation_simulation_snp.py -m seq_$name.0.poly -s $ref -N seq_$name > ../sequences/seq_$name.0.fa
	python mutation_simulation_snp.py -m seq_$name.1.poly -s $ref -N seq_$name > ../sequences/seq_$name.1.fa
done

rm *.poly
