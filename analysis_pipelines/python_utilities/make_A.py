import sys
import numpy

N=int(sys.argv[1])

A=numpy.zeros( (N, N) )
f=numpy.zeros( (N) )

names={}

for line in sys.stdin:
	l=line.split('\t')
	if l[0]=="@SAMPLE_X":
		break

for line in sys.stdin:
	l=line.strip('\n').split('\t')
	if len(l)<7:
		break
	l=map(float, l[:-1])
	try:
		x=names[l[0]]
	except:
		x=len(names.keys())
		names[l[0]]=x
	try:
		y=names[l[1]]
	except:
		y=len(names.keys() )
		names[l[1]]=y
	A[x,y] = l[6]	
	A[y,x] = l[6]	
	A[x,x] += (1.+l[2])/2.
	A[y,y] += (1.+l[4])/2.

for x in range(0, N):
	A[x,x]=A[x,x]/N
#	f[x]=A[x,x]

#f_bar=numpy.mean(f)

for x in range (0, N):
	for y in range (0, N-1):
		print str(A[x][y] )+", ",
	print str(A[x][N-1] )
