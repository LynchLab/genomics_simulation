import sys
import numpy
from scipy.optimize import minimize

n=10

#dI=numpy.diag([1,2,3,4,5,6,7,8,9,10])

dI=numpy.full( (n,n), -1./(n-1.) )
numpy.fill_diagonal(dI, 1)

def LL(params):
	sigma_e=1
	global dI, Y, n
#	print "HI!", sigma_e
#	print dI
	k=n-1
	V=sigma_e**2*dI

	D, E= numpy.linalg.eig(V)


	idx = (D.argsort()[::-1])
	D = numpy.real(D[idx])
	E = numpy.real(E[:,idx])


#	print numpy.matmul(numpy.matmul(E, numpy.diag(D) ), numpy.linalg.inv(E) )

	P=0

	#if min(sigma_e, sigma_d, sigma_a)<0:
	#	return sys.float_info.max

	#if min(a,d,e)<-0.1:
	#	P=min(a,d,e)**2
		#return sys.float_info.max
#	if min(D[0:k])<10**-9:
#		return sys.float_info.max


	D = numpy.diag(D)
	#E = numpy.delete(E, -1, axis=0)
	E = numpy.delete(E, -1, axis=1)
	D = numpy.delete(D, -1, axis=0)
	D = numpy.delete(D, -1, axis=1)

#	print numpy.matmul(numpy.matmul(E, D), numpy.linalg.inv(E) )

#	print D
#	print E

#	print "Min E:"+str(numpy.ndarray.min(E))
#	print "Min D:"+str(numpy.ndarray.min(numpy.diag(D)))
#	print "Min pi:"+str(numpy.ndarray.min((2*numpy.pi*numpy.diag(D))) )
#	print "Min Log:"+str(numpy.ndarray.min(numpy.log(2*numpy.pi*numpy.diag(D))) )

#	giV=numpy.matmul(numpy.matmul(E, numpy.linalg.inv(D) ), numpy.transpose(E) ) 
#	print giV
	giV=numpy.linalg.pinv(dI)
	print giV

	#print giV
	#giV=numpy.matmul(numpy.matmul(numpy.transpose(E), numpy.linalg.inv(D) ), E )  #Wrong way
	#print "D", D
	#print "D", numpy.diag(D)
	#print "Det", numpy.linalg.det(D)

	lnD=numpy.sum(numpy.log(2*numpy.pi*numpy.diag(D)))
#	print "lnD: ", -lnD/2
#	print "lnD: ", -numpy.matmul(numpy.matmul(numpy.transpose(Y), giV), Y )[0,0] /2
#	ret =-( (-lnD-numpy.matmul(numpy.matmul(numpy.transpose(Y), giV), Y)[0,0] ) /2)
	ret = (numpy.matmul(numpy.transpose(Y), Y)/sigma_e**2)[0,0]/2.- n*(-numpy.log(2)-numpy.log(numpy.pi)-numpy.log(sigma_e) )/2.
	#print (numpy.matmul(numpy.transpose(Y), Y)/sigma_e**2)[0,0]
#	print "ret:"+str(numpy.real(ret) )
	if (not(numpy.isfinite(ret))):
#		print "Whoops!"
		return sys.float_info.max
	return numpy.real(ret)

Y=numpy.asmatrix([0,1,2,3,4,5,6,7,8,9]).reshape((n,1))
print LL([1])
Y=numpy.asmatrix([0,2,4,6,8,10,12,14,16,18]).reshape((n,1))
print LL([1])


#Y=numpy.asmatrix([1,2,3,4,5,6,7,8,9,10]).reshape((n,1))
#print LL([1])
