import math
import sys
import itertools

#if len(sys.argv)==2:
#	print sys.argv	
#	File=open(sys.argv[1])
#else:
File=sys.stdin



D=[]
this_run=[]
line=File.readline()

N=len(map(int, line.split('\t')[1:] ))/2

for x in range(0, N):
	this_run.append(0)
	D.append({})

for line in File:
	S=map(int, line.split('\t')[1:])
	s=[]
	for x in range(0, N*2, 2):
		s.append(int(S[x]==S[x+1]))
	for x in range(0, N):
		if (s[x]==0):
			try:
				D[x][this_run[x]]+=1;
				this_run[x]=0;
			except:
				D[x][this_run[x]]=1;
				this_run[x]=0;
		else:
			try:
				this_run[x]+=1;
			except:
				this_run[x]=1;

print "x length count"
for x in range (0, N):
	for size in D[x].keys():
		print x, size, D[x][size]
