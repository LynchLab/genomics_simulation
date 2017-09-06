import sys
import numpy

N=int(sys.argv[1])

A=numpy.zeros( (N, N) )
f=numpy.zeros( (N) )

names={}

for line in sys.stdin:
	l=line.split()
	if l[0]=="FID1":
		break

for line in sys.stdin:
	l=line.strip('\n').split()
	if len(l)<7:
		break
	try:
		x=names[l[1]]
	except:
		x=len(names.keys())
		names[l[1]]=x
	try:
		y=names[l[3]]
	except:
		y=len(names.keys() )
		names[l[3]]=y
	A[x,y] = l[9]	
	A[y,x] = l[9]

for x in range (0, N):
	for y in range (0, N-1):
		print str(A[x][y] )+", ",
	print str(A[x][N-1] )
