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
	l=line.split('\t')
	if l[0]=="@SAMPLE_X":
		break

for line in sys.stdin:
	l=line.strip('\n').split('\t')
	if len(l)<5:
		break
	l=map(float, l)
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
	if x < n and y < n:
		f=float(l[9])
		fx=float(l[5])
		fy=float(l[6])
		A[x,y] = f*(-fx-fy)
		A[y,x] = f*(-fx-fy)

#f_bar=numpy.mean(f)

for x in range (0, n):
	for y in range (0, n-1):
		print str(A[x][y] )+", ",
	print str(A[x][n-1] )
