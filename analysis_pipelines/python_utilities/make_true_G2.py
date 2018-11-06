import numpy
import sys

x=10
File=open(sys.argv[1])
File.readline()
File.readline()
N=0
lines=0

for line in File:
	line=line.split()
	m=map(float, line[3:])
	a=float(line[1])
	d=float(line[2])
	if max(m)>0:
		if (N==0):
			N=int(len(m)/2)
			G=numpy.zeros((N,N))
			mu_M=numpy.zeros((N,1))
			mu_H=numpy.zeros((N,1))
		A=[i * a for i in map(float, line[3:(N+3)]) ]
		mu=float(reduce(lambda x, y: x + y, A) )/ float(N)
		A=[i - mu for i in A ]
		D=[i * d for i in map(float, line[(N+3):]) ] 
		mu=float(reduce(lambda x, y: x + y, D) ) / float(N)
		D=[i - mu for i in D ]
		A=numpy.matrix(A).reshape((1,N))
		D=numpy.matrix(D).reshape((1,N))
		G+=A.transpose()*D+D.transpose()*A
		mu_M+=A.transpose()
		mu_H+=D.transpose()
		lines+=1

mu_M/=lines
mu_H/=lines

G=G-2*mu_M*mu_H.transpose()-2*mu_H*mu_M.transpose()

G=G/numpy.mean(numpy.diag(G))

G=G.tolist()

for line in G:
	for x in line[:-1]:
		print str(x)+",",
	print str(line[-1])
#print true_v
