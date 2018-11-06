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
		G=numpy.zeros((x,x))
		mu_m=numpy.zeros((x,1))
		mu_h=numpy.zeros((x,1))
	if max(m)>0:
		m=numpy.matrix(map(float, line[1:(x+1)]) ).reshape((1,x))
		h=numpy.matrix(map(float, line[(x+1):(x*2+1)]) ).reshape((1,x))
		G+=m.transpose()*h+h.transpose()*m;
		mu_m+=m.transpose()
		mu_h+=h.transpose()
		lines+=1
mu_m/=lines
mu_h/=lines

G=G-2*mu_m*mu_h.transpose()-2*mu_h*mu_m.transpose()

G=G/numpy.mean(numpy.diag(G))

G=G.tolist()

for line in G:
	for x in line[:-1]:
		print str(x)+",",
	print str(line[-1])
#print true_v
