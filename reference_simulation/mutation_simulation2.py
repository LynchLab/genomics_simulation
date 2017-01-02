import random
import sys
import argparse
import mutation

base=['a','c','g','t']

parser = argparse.ArgumentParser(description='make a mutation file.')
parser.add_argument('-n', metavar='--number', type=int, default=5000,
                        help='number of mutations simulated')
parser.add_argument('-l', metavar='--length', type=int, default=10000,
                        help='starting sequence length')
parser.add_argument('-o', metavar='--onlysnp', type=bool, default=False,
                        help='only simulate snp')
parser.add_argument('-S', metavar='--size', type=int, default=0,
                        help='final size')

args = parser.parse_args()

Time=args.n
seq_l=args.l
size=args.S

if (not args.o):
	kwts=[0.00005,	0.002,	0.002,	0.004,	0.002,	0.0001,	0.0,	1.]	
else:
	kwts=[0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.0,	1.]	
	
#twts=[0.2,	0.1,	0.1,	0.2,	0.1,	0.2,	0.0,	1.]	
wfns=["denovo", "copy", "tdcopy", "trans", "invert", "delete", "static",	"SNP"]
wts=[kwts[0]]

for k in kwts[1:]:
	wts.append(k+wts[-1])	
t=0
while t<Time or seq_l<size:
	if (seq_l>0):
		INT_A=random.randint(0, seq_l-1 )
		INT_C=random.randint(0, seq_l-1 )
	else:
		INT_A=0
		INT_C=0
	INT_B=INT_A+random.randint(1, 1000)
	BOOL=False
	N=random.randint(2,10)
	FNS_FLOAT=random.random()
	x=0
	while(wts[x]<FNS_FLOAT):
		x+=1
	if x==7:
		INT_B=random.randint(1, 3)
	seq_l=mutation.fns_size[x](seq_l, INT_A, INT_B, INT_C, BOOL, N)
	print wfns[x], INT_A+1, INT_B, INT_C, BOOL, N, seq_l
	t+=1
