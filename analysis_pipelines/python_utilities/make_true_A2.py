import numpy
import sys

x=10
M=[]
File=open(sys.argv[1])
File.readline()
File.readline()
N=0

for line in File:
	line=line.split()
	m=map(float, line[3:])
	a=float(line[1])
	d=float(line[2])
	if max(m)>0:
		if (N==0):
			N=int(len(m)/2)
		A=[i * a for i in map(float, line[3:(N+3)]) ]
		mu=float(reduce(lambda x, y: x + y, A) )/ float(N)
		A=[i - mu for i in A ]
		D=[i * d for i in map(float, line[(N+3):]) ] 
		mu=float(reduce(lambda x, y: x + y, D) ) / float(N)
		D=[i - mu for i in D ]
		M.append(A+D)


M=numpy.matrix(M)
m=M.mean(0)
M=M-m
v=numpy.dot(M.transpose(), M )

A=v[0:N, 0:N]

A=A/numpy.mean(numpy.diag(A))

A=A.tolist()

for line in A:
	for x in line[:-1]:
		print str(x)+",",
	print str(line[-1])
