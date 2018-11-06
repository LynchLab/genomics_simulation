import numpy
import sys

x=0
lines=0

File=open(sys.argv[1])
File.readline()
File.readline()

for line in File:
	line=line.split()
	m=map(float, line[1:])
	if (x==0):
		x=(len(line)-1)/2
		A=numpy.zeros((x,x))
		J=numpy.ones((x,x))
		a=numpy.zeros((1,x))
	if max(m)>0:
		m=numpy.matrix(map(float, line[1:(x+1)]) )
		A+=m.transpose()*m
		a+=m
		lines+=1
a/=lines
a=a.transpose()

A=A-2.*a*a.transpose()

A=A/numpy.mean(numpy.diag(A))

A=A.tolist()

for line in A:
	for x in line[:-1]:
		print str(x)+",",
	print str(line[-1])
#print true_v
