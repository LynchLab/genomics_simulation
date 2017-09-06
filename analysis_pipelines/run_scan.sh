for a in 0 1
do
for b in 0 1
do
for c in -0.5 0 0.5
do
for x in {0..3}
do 
for y in 1 2 5
do 
for z in {1..100}
do 
echo $x $y
k=`echo "$y*10^$x" | bc`
zcat ../sequences/states_bin.txt.gz | python-2.7.9 -u states_to_pheno_gamma.py ../sequences/name-file.txt 2000 $k ../analysis_files/inbred $a $b $c > /dev/null 2>> temp${a}_${b}_${c}
done
done
done
done 
done 
done
