import sys
import numpy
from scipy.optimize import minimize

VA=[]
File=open(sys.argv[1]+"A.csv")
for line in File:
	line=map(float, line.split(',') )
 	VA.append(line)
n=len(VA)
VA=numpy.asmatrix(VA).reshape((n, n))

File=open(sys.argv[1]+"D.csv")
VD=[]
for line in File:
	line=map(float, line.split(',') )
 	VD.append(line)
VD=numpy.asmatrix(VD).reshape((n, n))

File=open(sys.argv[1]+"G.csv")
VG = []
for line in File:
	line=map(float, line.split(',') )
 	VG.append(line)
VG=numpy.asmatrix(VG).reshape((n, n))

dI=numpy.diag(numpy.repeat(1, n) )#-1/(n)


File=open(sys.argv[1]+"t_final.txt")
X = []
for line in File:
	if line[0]!='@':
		line=map(float, line.split('	') )
 		X.append(line[1])
X=numpy.asmatrix(X).reshape((n, 1))
#E=
def LL(params):
	a,d,g=params
	print a, d, g, 1-a-d-g
	global VA, VD, VG, X, n
	k=n-1
	V=a*VA+d*VD+g*VG+(1-a-d-g)*dI
	D, E= numpy.linalg.eig(V)

	idx = D.argsort()[::-1]   
	D = D[idx]
	E = E[:,idx]

	#print D
	D = numpy.diag(D)
	#print E
	#print D
	E = numpy.delete(E, -1, axis=0)
	E = E.reshape(n,k)
	D = numpy.delete(D, -1, axis=0)
	D = numpy.delete(D, -1, axis=1)
	D = D.reshape(k,k)
	#print numpy.linalg.det(D)
	#print E
	#print D
	C=numpy.matrix([[10]])
	giV=numpy.matmul(numpy.matmul(E, numpy.linalg.inv(D) ), numpy.transpose(E) ) 
	return ((-numpy.log(numpy.linalg.det(D))-numpy.matmul(numpy.matmul(numpy.transpose(X), giV), X ) ) /2)[0,0]
#print VD
#print VA

#res = minimize(LL, [0.5,0.5,0], method='nelder-mead', options={'xtol': 1e-8, 'disp': True})
for x in range(0, 100):
	for y in range(x, 100):
		print LL([x/100.,(y-x)/100.,0])
		
#print res

#print LL([1,0,0,0])
