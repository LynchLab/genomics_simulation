POPULATION=200
POP2=$((POPULATION/2))

SAMPLE=50
SAMP=`seq 0 1 $((POPULATION-1)) | shuf | head -$SAMPLE`

for s in $SAMP
do
SAMPLE_CHRM+=","
SAMPLE_CHRM+=$((s*2+2))
SAMPLE_CHRM+=","
SAMPLE_CHRM+=$((s*2+3))
SAMPLE_NAME+=","
SAMPLE_NAME+=$((s+2))
done

TIME=$((1*POPULATION))
TIMEX=$((9*PULATION+POPULATION/2))
REF="reference.fa"
COV=3
K=100
SIZE=$(($K*640))        #make sure to use bwa mem!!
SNPS=$(($K*10))
REFTYPE='Y'

LD_DIST=1000

A=1
D=0
E=0.5
