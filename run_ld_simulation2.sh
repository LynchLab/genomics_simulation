#!/bin/bash
source settings.sh

rm -rf ./analysis_files/estimates.txt


state_file=ld_state.txt
genmap=
#@NS:15427	GS:15627	BS:1418352456

n=19200
#n=7500
N=7500

S=1000
V=32000

G=5000
T=$((G-3))
C=1000
#C=2

cd  analysis_files2

pedigree_sim -T $T -y g -N $n -1 $N -s $S -g $G -c $C -plnGt -k $G -v $V -e 0.005  > $state_file

cd ../analysis_pipelines/
echo $((G*n+n-1))
echo python count_sub_pedigree.py ../analysis_files2/pedigree.txt $((G*n+n-1)) 2 $n 
python count_sub_pedigree.py ../analysis_files2/pedigree.txt $((G*n+n-1)) $((G-1)) $n > junk3
../gsl_correl/call_relatedness_lambda -i ../analysis_files2/$state_file -c $C -r 5 -? > junk #2> junk2
#../gsl_correl/call_relatedness_count -i ../analysis_files2/$state_file > junk_count #2> junk2
