import random
import sys
import argparse
import mutation

parser = argparse.ArgumentParser(description='simulates a reference by sequential mutations.')
parser.add_argument('-n', metavar='--number', type=int, default=10000,
                        help='number of mutations simulated')
parser.add_argument('-s', metavar='--sequence', type=str,
                   help='sequence file')
parser.add_argument('-m', metavar='--mutation', type=str,
                   help='mutation file')
parser.add_argument('-N', metavar='--name', type=str, default="simulated_sequence",
                   help='sequence name')

args = parser.parse_args()

seq=args.s
mut=args.m
seq_name=args.N

File=open(seq)
seq=''.join(File.read().split('\n')[1:])
File.close()

File=open(mut)

for line in File:
	[INT_A, INT_B, INT_C, BOOL, N]=[0,0,0,False,0]
	fields=line.split()
	name=fields[0]
	INT_A=fields[1]
	if len(fields)>2:
		INT_B=fields[2]
	if len(fields)>3:
		INT_C=fields[3]
	if len(fields)>4:
		BOOL=fields[4]
	if len(fields)>5:
		 N=fields[5]
	if len(fields)>6:
		 size=fields[6]
	INT_A, INT_B, INT_C, N=map(int, [INT_A, INT_B, INT_C, N])
	seq=mutation.fnsk[name](seq, INT_A-1, INT_B, INT_C, BOOL, N)		

print ">"+seq_name+":"+str(len(seq))
for x in range(0, len(seq), 90):
	print seq[x:(x+90)]
