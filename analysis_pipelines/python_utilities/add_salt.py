import sys
import numpy

PATH=sys.argv[1]

k=float(sys.argv[2])

name=[]
File=open(PATH+"name-file.txt")
for line in File:
	if line[0]!='@':
		line=line.split()
		name+=line[1:]

n = len(name)

#print n
#print name

Beta = numpy.asmatrix([0.0]*n).reshape((n,1))
Delta = numpy.asmatrix([0.0]*n).reshape((n,1))

A = numpy.asmatrix([0.0]*n).reshape((n,1))
D = numpy.asmatrix([0.0]*n).reshape((n,1))

File=open(PATH+"beta_delta.txt")

x=0
for line in File:
	if line[0]!='@':
		line=map(float, line.split('	')[1:] )
		Beta[x,0]=line[1]
		Delta[x,0]=line[2]
		A[x,0]=line[3]
		D[x,0]=line[4]
		x+=1

G=A+D
VarG=numpy.var(G)
sdE=numpy.sqrt(VarG/k-VarG)

E=numpy.asmatrix(numpy.random.normal(0, sdE, n)).reshape((n,1))
Z=G+E

for x in range(0, n):
	print name[x]+'\t'+name[x]+'\t'+str(Z[x,0])

File=open(PATH+"true_var.hsq", 'w')

File.write("Var(Beta)/Var(Z)\t"+str(numpy.var(Beta)/numpy.var(Z))+'\n')
File.write("Var(Delta)/Var(Z)\t"+str(numpy.var(Delta)/numpy.var(Z))+'\n')
File.write("2*Cov(Beta,Delta)/Var(Z)\t"+str(2*numpy.cov(Beta.transpose(),Delta.transpose())[0,1]/numpy.var(Z))+'\n')
File.write("\n")
File.write("Var(A)/Var(Z)\t"+str(numpy.var(A)/numpy.var(Z))+'\n')
File.write("Var(D)/Var(Z)\t"+str(numpy.var(D)/numpy.var(Z))+'\n')
File.write("2*Cov(A,D)/Var(Z)\t"+str(2*numpy.cov(A.transpose(),D.transpose())[0,1]/numpy.var(Z))+'\n')
File.write("\n")
File.write("2*Cov(G,E)/Var(Z)\t"+str(2*numpy.cov(E.transpose(),G.transpose())[0,1]/numpy.var(Z))+'\n')
File.write("Var(E)/Var(Z)\t"+str(numpy.var(E)/numpy.var(Z))+'\n')
