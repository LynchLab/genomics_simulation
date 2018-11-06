import sys
import numpy
from scipy.optimize import minimize

PATH=sys.argv[1]

VA=[]
File=open(PATH+"A.csv")
for line in File:
	line=map(float, line.split(',') )
 	VA.append(line)
n=len(VA)
VA=numpy.asmatrix(VA).reshape((n, n))

File=open(PATH+"D.csv")
VD=[]
for line in File:
	line=map(float, line.split(',') )
 	VD.append(line)
VD=numpy.asmatrix(VD).reshape((n, n))

File=open(PATH+"G.csv")
VG = []
for line in File:
	line=map(float, line.split(',') )
 	VG.append(line)
VG=numpy.asmatrix(VG).reshape((n, n))


File=open(PATH+"mua.csv")
mua = []
for line in File:
	line=map(float, line.split(',') )
 	mua.append(line)
mua=numpy.asmatrix(mua).reshape((n, 1))

File=open(PATH+"mud.csv")
mud = []
for line in File:
	line=map(float, line.split(',') )
 	mud.append(line)
mud=numpy.asmatrix(mud).reshape((n, 1))

dI=numpy.full( (n,n), -1./(n-1.) )
numpy.fill_diagonal(dI, 1)
#dI=numpy.diag(numpy.repeat(1., n) )-1./(n-1)


File=open(PATH+"plink.pheno")
X = []
for line in File:
	if line[0]!='@':
		line=map(float, line.split('	') )
 		X.append(line[2])
X=numpy.asmatrix(X).reshape((n, 1))
#E=
def LL(params):
	a, d, g, e=params
	#print a, d, 1.-a-d-e, e
	global VA, VD, VG, Y, n
	k=n-1
	V=a*VA+d*VD+g*VG+e*dI
	D, E= numpy.linalg.eig(V)

	idx = D.argsort()[::-1]   
	D = D[idx]
	E = E[:,idx]

	#print min(D[0:k])

	P=0

	if min(a,d,e)<-0.1:
		P=min(a,d,e)**2
		#return sys.float_info.max
	if min(D[0:k])<10**-9:
		return sys.float_info.max


	D = numpy.diag(D)
	#E = numpy.delete(E, -1, axis=0)
	E = numpy.delete(E, -1, axis=1)
	D = numpy.delete(D, -1, axis=0)
	D = numpy.delete(D, -1, axis=1)

	C=numpy.matrix([[10]])
	giV=numpy.matmul(numpy.matmul(E, numpy.linalg.inv(D) ), numpy.transpose(E) ) 
	#giV=numpy.matmul(numpy.matmul(numpy.transpose(E), numpy.linalg.inv(D) ), E )  #Wrong way
	#print giV
	ret =-((-numpy.log(numpy.linalg.det(D))-numpy.matmul(numpy.matmul(numpy.transpose(Y), giV), Y) ) /2)[0,0]
	#print ret
	if (not(numpy.isfinite(ret))):
		return sys.float_info.max
	return ret+P
#print VD
#print VA

A = numpy.hstack([mua, mud])

print A.shape
print X.shape

Y=X-numpy.mean(X)
Y=Y/numpy.std(Y)

mua_bar, mud_bar=numpy.linalg.lstsq(A, Y)[0]

Z=Y[:]
Y=Y-mua_bar[0,0]*mua-mud_bar[0,0]*mud
v=numpy.std(Y)**2


print mua+mud
print
print Z

quit()
#print v
#regress Y on mua and mud

#for j in range(0, 10):
#	Z=numpy.asmatrix(numpy.random.normal(0, 0.2*j, n)).reshape((n,1))
#	Y=X+Z
#	Y=Y-numpy.mean(Y)
#	Y=Y/numpy.std(Y)

#res1 = minimize(LL, [1./3.*v,1./3.*v,0,1./3.*v], method='SLSQP', tol=10**-9)
#a,d,g,e=res1.x

A=numpy.var((mua_bar[0,0]*mua).T)+a
D=numpy.var((mud_bar[0,0]*mud).T)+d
G=numpy.cov((mua_bar[0,0]*mua).T, (mud_bar[0,0]*mud).T)+g
#E=e
#print
print  A, a, D, d, #, G, E

