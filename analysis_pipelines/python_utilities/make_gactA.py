import sys
import numpy

N=int(sys.argv[1])

if len(sys.argv)==3:
	n=int(sys.argv[2])
else:
	n=N

A=numpy.zeros( (n, n) )

names={}

for line in sys.stdin:
	l=line.strip('\n').split('\t')
	l=map(float, l)

	x=int(l[0]-1)
	y=int(l[1]-1)

	if x<n and y<n:
		A[x,y] = l[3]	
		A[y,x] = l[3]	

for x in range (0, n):
	for y in range (0, n-1):
		print str(A[x][y] )+", ",
	print str(A[x][n-1] )
