import numpy
import sys

x=10
M=[]
File=open(sys.argv[1])
File.readline()
File.readline()

for line in File:
	line=line.split()
	m=map(float, line[1:])
	if max(m)>0:
		M.append(map(float, line[1:]) )

M=numpy.matrix(M)
v=numpy.cov(M.transpose() )

x=(len(line)-1)/2
A=v[0:x, 0:x]
D=v[x:,x:]
G=v[x:,0:x]+v[0:x,x:]

A=A/numpy.mean(numpy.diag(A))
D=D/numpy.mean(numpy.diag(D))
G=G/numpy.mean(numpy.diag(G))

G=G.tolist()

for line in G:
	for x in line[:-1]:
		print str(x)+",",
	print str(line[-1])
#print true_v
