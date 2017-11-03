import numpy
import sys

x=10
M=[]
File=open(sys.argv[1])
File.readline()
File.readline()

for line in File:
	line=line.split()
	M.append(map(float, line[1:]) )

M=numpy.matrix(M)
v=numpy.cov(M.transpose() )

x=(len(line)-1)/2
print v.shape
A=v[0:x, 0:x]
print A.shape
D=v[x:,x:]
print D.shape
G=v[0:x,x:]

A=A/numpy.mean(numpy.diag(A))
D=D/numpy.mean(numpy.diag(D))
G=G/numpy.mean(numpy.diag(G))

for line in A:
	for x in line:
		print str(x)+",",
	print
#print true_v
