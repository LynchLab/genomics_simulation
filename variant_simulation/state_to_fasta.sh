ref=$1
states=$2
SIZE="`head -1 $states | grep '	' -o | wc -l`"
POLY=../sequences/polymorphisms.map
REFSIZE=$((`tail -n +2 $1 | wc -c`-`tail -n +2 $1 | wc -l`))

python mutation_simulation2.py -n "`wc -l $states | cut -d ' ' -f 1`" -l $REFSIZE -o True > temp
python mutation_sort.py temp > $POLY

rm temp
rm -f poly.db
python state_to_fasta.py -N $((SIZE/2)) -s $states -m $POLY 

for N in `seq 0 2 $((SIZE-2))` 
do
	name=`printf %03d $((N/2))`
	echo $name

	echo "SELECT var FROM snps WHERE sample=$N;" | sqlite3 poly.db  | ./mutate -r $ref -n seq_$name  | gzip - > ../sequences/seq_$name.0.fa.gz
	echo "SELECT var FROM snps WHERE sample=$((N+1));" | sqlite3 poly.db  | ./mutate -r $ref -n seq_$name | gzip - > ../sequences/seq_$name.1.fa.gz
done

rm -f poly.db
