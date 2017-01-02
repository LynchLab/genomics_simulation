import math
import sys
from scipy import stats
import numpy
import itertools

if len(sys.argv)==2:
	print sys.argv	
	File=open(sys.argv[1])
else:
	File=sys.stdin


print "run freq rho1 rho2 rho3 P OX EX OY EY"

D=[]
Sbuffer=[]
line=File.readline()
DIST=2


N=len(map(int, line.split('\t')[1:] ))/2
XMAX=N

for d1 in range(0, N):
	D.append([])
	for X in range (0, XMAX):
		D[d1].append([])
		for j in range (0, 2):
			D[d1][X].append(0)
for line in File:
	S=map(int, line.split('\t')[1:])
	s=[]
	for x in range(0, N*2, 2):
		s.append(int(S[x]!=S[x+1]))
	for x in range(0, XMAX):
		i=s[x]
		d1=sum(s)-i
		D[d1][x][i]+=1

for x in range (0, XMAX):
	for c1 in range(1, N-1):
		Den=sum(D[c1][x])
		if (Den>0):
			OX=float(D[c1][x][1])/float(Den)
			print x, float(c1)/float(N-1), OX
