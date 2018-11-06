import sys
import numpy
import scipy
from scipy.spatial.distance import mahalanobis
from scipy.optimize import minimize

#PATH=sys.argv[1]
PATH="../analysis_files/"

File=open(PATH+"mapgd_relatedness.out")
X = []
for line in File:
	line=map(float, line.split('\t') )
 	X.append(line)
n=len(X)

X=(numpy.asmatrix(X))

S1=numpy.asmatrix(X[0:n,0:n]).reshape((n,n))
S2=numpy.asmatrix(X[0:n, n:(2*n)]).reshape((n,n))
S3=numpy.asmatrix(X[0:n, (2*n):(3*n)]).reshape((n,n))
S4=numpy.asmatrix(X[0:n, (3*n):(4*n)]).reshape((n,n))
S5=numpy.asmatrix(X[0:n, (4*n):(5*n)]).reshape((n,n))

dI=numpy.full( (n,n), -1./(n-1.) )
numpy.fill_diagonal(dI, 1.)
#dI=numpy.identity(n)

File=open(PATH+"plink.pheno")

X = []

for line in File:
	if line[0]!='@':
		line=map(float, line.split('	') )
 		X.append(line[2])

X=numpy.asmatrix(X).reshape((n, 1))

def LL3(params):
	mua, mud, sigma_a, sigma_d, sigma_ad, sigma_e=params
	#mud, mua, sigma_a, sigma_d, sigma_e=params
	global S1, S2, S3, S4, S5, dI, Y, Z, n

	p=2./(1.+numpy.exp(-sigma_ad) )-1.

	sa=numpy.exp(sigma_a)
	sd=numpy.exp(sigma_d)
	se=numpy.exp(sigma_e)

	V=mua*S1+mud*S2+(mua*mud+sd*sa*p)*S3+(mua**2+sa)*S4+(mud**2+sd)*S5+se*dI
	sY=Y
#	V=mua*S1+mud*S2+sigma_ad*S3+sigma_a*S4+sigma_d*S5+sigma_e*dI


#	distval = mahalanobis(Y,Z,V)
#	          scipy.spatial.distance.mahalanobis
#mahalanobis(x, center = mean, cov = sigma)
#	D = numpy.linalg.eigvals(V)
#	D=sorted(D, reverse=True)
#	k=sum(i > 1e-10 for i in D)
#	D=D[0:k]

#	logdet = sum(numpy.log(D) )
#	logretval = -(n * numpy.log(2. * numpy.pi) + logdet + distval)/2.

#	print k
#	print distval
#	print logdet
#	print logretval
#
#	return(logretval)


	try:
		D = numpy.linalg.eigvals(V)
	except:
		print "Can't get eigen values..."
		return sys.float_info.max
	D=sorted(D, reverse=True)
	
	print min(D)
	k=sum(i > 1e-8 for i in D)

	print k
	#return D

	if k!=n-1:
		return sys.float_info.max

	D=numpy.array(D[0:k])

	#print "k ", len(D)

	print mua*k**2, mud*k**2, sigma_a**2*k, sigma_d**2*k, sigma_ad, sigma_e**2*k
	#print mua/k, mud/k, (sigma_a**2-mua**2)/k**2, (sigma_d**2-mud**2)/k**2, sigma_ad/k**2, sigma_e/k**2

	if max(abs(mua), abs(mud) )>10:
		return sys.float_info.max
		
	giV=numpy.linalg.pinv(V, rcond=1e-5)

	rank=numpy.linalg.matrix_rank(giV)
	print "k:", k, " rank:", numpy.linalg.matrix_rank(giV)
	if (rank!=k):
		return sys.float_info.max

	lnD=numpy.sum(numpy.log(2*numpy.pi*D))

	ret =-(-lnD-(numpy.matmul(numpy.matmul(numpy.transpose(sY), giV), sY)[0,0] ) /2)
	print "ret: "+str(numpy.real(ret) )
	print "lnD: "+str(lnD)
	#print "min(L): "+ str(numpy.ndarray.min(numpy.diag(D)))
	print "e^-: "+str( (numpy.matmul(numpy.matmul(numpy.transpose(sY), giV), sY)[0,0] ) )
	if (not(numpy.isfinite(ret))):
		print "Whoops!"
		return sys.float_info.max
	return numpy.real(ret)

def LL(params):
	mua, mud, sigma_a, sigma_d, sigma_ad, sigma_e=params
	#mud, mua, sigma_a, sigma_d, sigma_e=params
	global S1, S2, S3, S4, S5, dI, Y, Z, n

	p=2./(1.+numpy.exp(-sigma_ad) )-1.
	sa=numpy.exp(sigma_a)
	sd=numpy.exp(sigma_d)
	se=numpy.exp(sigma_e)

	V=mua*S1+mud*S2+(mua*mud+sd*sa*p)*S3+(mua**2+sa)*S4+(mud**2+sd)*S5+se*dI

	DEN=numpy.sum(numpy.diag(V))
	V=V/DEN
	sY=Y/numpy.sqrt(DEN)
#	print DEN
#	V=mua*S1+mud*S2+sigma_ad*S3+sigma_a*S4+sigma_d*S5+sigma_e*dI
#	distval = mahalanobis(Y,Z,V)
#	          scipy.spatial.distance.mahalanobis
#mahalanobis(x, center = mean, cov = sigma)
#	D = numpy.linalg.eigvals(V)
#	D=sorted(D, reverse=True)
#	k=sum(i > 1e-10 for i in D)
#	D=D[0:k]
#
#	logdet = sum(numpy.log(D) )
#	logretval = -(n * numpy.log(2. * numpy.pi) + logdet + distval)/2.
#
#	print k
#	print distval
#	print logdet
#	print logretval
#
#	return(logretval)


	try:
		D = numpy.linalg.eigvals(V)
	except:
		print "Can't get eigen values..."
		return sys.float_info.max
	D=sorted(D, reverse=True)
	
#	print min(D)
	k=sum(i > 1e-8 for i in D)
#	print k
# 	return D

	if k!=n-1:
		return sys.float_info.max

	D=numpy.array(D[0:k])

	#print "k ", len(D)

	#print mua/k, mud/k, (sigma_a**2-mua**2)/k**2, (sigma_d**2-mud**2)/k**2, sigma_ad/k**2, sigma_e/k**2

	if max(abs(mua), abs(mud) )>10:
		return sys.float_info.max
		
	giV=numpy.linalg.pinv(V, rcond=1e-5)

	rank=numpy.linalg.matrix_rank(giV)
	if (rank!=k):
		return sys.float_info.max

	lnD=numpy.sum(numpy.log(2*numpy.pi*D) )+numpy.log(DEN)*k

	ret =-(-lnD-(numpy.matmul(numpy.matmul(numpy.transpose(sY), giV), sY)[0,0] ) /2)
#	print mua*k**2, mud*k**2, sigma_a**2*k, sigma_d**2*k, sigma_ad, sigma_e**2*k
	print "ret:", numpy.real(ret), "vals:", sa*numpy.sqrt(k)/numpy.sqrt(2), sd*numpy.sqrt(k)/numpy.sqrt(2), 2./(1.+numpy.exp(sigma_ad))-1., se*numpy.sqrt(k)/numpy.sqrt(2)
#	print "lnD: "+str(lnD)
#	print "min(L): "+ str(numpy.ndarray.min(numpy.diag(D)))
#	print "e^-: "+str( (numpy.matmul(numpy.matmul(numpy.transpose(sY), giV), sY)[0,0] ) )
	if (not(numpy.isfinite(ret))):
		print "Whoops!"
		return sys.float_info.max
	return numpy.real(ret)

def LL2(params):
	sigma_a, sigma_d, sigma_ad, sigma_e=params
	return LL([0., 0., sigma_a, sigma_d, sigma_ad, sigma_e])

Y=X-numpy.mean(X)
Y=Y/numpy.std(Y)
Z=numpy.zeros((n,1), dtype=numpy.float )
v=numpy.std(Y)

print v

res1 = minimize(LL2, [-10, -10, 2.99443, 0.0], method='SLSQP')
#res2 = minimize(LL, [0.0, 0.0]+numpy.ndarray.tolist(res1.x), method='powell')

sigma_a, sigma_d, sigma_ad, sigma_e=res1.x
k=4800
print abs(sigma_a)*numpy.sqrt(k)/numpy.sqrt(2), abs(sigma_d)*numpy.sqrt(k)/numpy.sqrt(2), 2./(1.+numpy.exp(sigma_ad))-1., abs(sigma_e)*numpy.sqrt(k)/numpy.sqrt(2)
print LL2(res1.x)
print LL2([res1.x[0], res1.x[1], 2.94443, res1.x[3]])

#res1 = minimize(LL, [0, 0, v/numpy.sum(numpy.diag(S4))/3., v/numpy.sum(numpy.diag(S5))/3., v/3.], method='SLSQP', tol=10**-9)
#print dI
#print res1.x[0]/numpy.sum(numpy.diag(S4))
#print res1.x[1]/numpy.sum(numpy.diag(S5))
#print res1.x[2]/numpy.sum(numpy.diag(dI))

