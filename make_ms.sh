file1=ms
file2=ms.state

segsites=50000
rho=5000  #32000
len=100000 #32000

N=2 #$((8*100))
#$ms 2 1 -t 50000 -r 5000 100000
 
../msdir/ms $N 1 -t $segsites -r $rho $len > $file1
