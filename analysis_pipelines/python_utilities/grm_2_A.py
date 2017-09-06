import sys
import numpy

N=int(sys.argv[1])

A=numpy.zeros( (N, N) )
f=numpy.zeros( (N) )

for line in sys.stdin:
	l=line.strip('\n').split('\t')

	x=int(l[0])-1
	y=int(l[1])-1
	try:
		A[x,y] = float(l[3])	
		A[y,x] = float(l[3])
	except:
		print l
		exit(0)

for x in range(0, N):
	A[x,x]=A[x,x]/N
#	f[x]=A[x,x]

#f_bar=numpy.mean(f)

for x in range (0, N):
	for y in range (0, N-1):
		print str(A[x][y] )+", ",
	print str(A[x][N-1] )
