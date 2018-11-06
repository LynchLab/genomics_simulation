import sys
import numpy
import scipy
import sympy
import time

from scipy.spatial.distance import mahalanobis
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from sklearn.linear_model import LinearRegression

#PATH=sys.argv[1]
l=float(sys.argv[1])
PATH="../analysis_files/"

File=open(PATH+"mapgd_relatedness.out")
X = []
for line in File:
	line=map(float, line.split('\t') )
 	X.append(line)
n=len(X)

X=(numpy.asmatrix(X))

S1=numpy.asmatrix(X[0:n,0:1]).reshape((n,1))
S2=numpy.asmatrix(X[0:n, 1:2]).reshape((n,1))
S3=numpy.asmatrix(X[0:n, (2+0*n):(2+1*n)]).reshape((n,n))
S4=numpy.asmatrix(X[0:n, (2+1*n):(2+2*n)]).reshape((n,n))
S5=numpy.asmatrix(X[0:n, (2+2*n):(2+3*n)]).reshape((n,n))

def t(p):
	return numpy.transpose(p)

SS11=S1*t(S1)
SS12=S1*t(S2)+S2*t(S1)
SS22=S2*t(S2)

dI=numpy.full( (n,n), -1./(n-1.) )
numpy.fill_diagonal(dI, 1.)

print dI
#dI=numpy.identity(n)

File=open(PATH+"plink.pheno")

X = []

for line in File:
	if line[0]!='@':
		line=map(float, line.split('	') )
 		X.append(line[2])

X=numpy.asmatrix(X).reshape((n, 1))

def getOIM(params):
	mua, mud, sigma_a, sigma_d, sigma_ad, sigma_e=params
	global S1, S2, S3, S4, S5, dI, Y, Z, n
	MA, MD, SA, SD, SAD, SE=sympy.symbols("MA MD SA SD SAD SE")
	V = SAD*S3+SA*S4+SD*S5+SE*dI
	u, s, v = sympy.svd_r(V)
	K=n-1
	giV=sympy.transpose(u[0:n,0:K]*sympy.transpose(sympy.transpose(v[0:K, 0:n])/s[0:K]))
	print giV
	quit()	


class CIL:
	def __init__ (self, x, params):
		self.x=abs(x)-1
		self.params=params
		if x>0:
			self.s=1
		else:
			self.s=-1
		self.max=LL(params)
	def get (self, d):
		params=self.params[:]
		params[self.x]+=abs(d[0])*self.s
		ret=abs(LL(params)-self.max-5)
#		print "HI!", params, ret
		return ret
	def got (self, d):
		params=self.params[:]
		params[self.x]+=abs(d[0])*self.s
		return params[self.x]

def CI(x, params):
	cil=CIL(x, params)
	res3 = minimize(cil.get, [0], method='powell' )
	return cil.got([res3.x])

def CIfromR5(params):
	mua, mud, sigma_a, sigma_d, sigma_ad=params
	print mua, CI(-1, params), CI(1, params)
	print mud, CI(-2, params), CI(2, params)
	print numpy.sqrt(numpy.exp(sigma_a)),  numpy.sqrt(numpy.exp(CI(-3, params))), numpy.sqrt(numpy.exp(CI(3, params)))
	print numpy.sqrt(numpy.exp(sigma_d)),  numpy.sqrt(numpy.exp(CI(-4, params))), numpy.sqrt(numpy.exp(CI(4, params)))
	print 2./(1.+numpy.exp(-sigma_ad) )-1.,  2./(1.+numpy.exp(-CI(-5, params) ) )-1.,  2./(1.+numpy.exp(-CI(5, params) ) )-1.

def fromR5(params):
	mua, mud, sigma_a, sigma_d, sigma_ad=params
	return mua, mud, numpy.sqrt(numpy.exp(sigma_a)), numpy.sqrt(numpy.exp(sigma_d)), 2./(1.+numpy.exp(-sigma_ad) )-1. 


def fromR3(params):
	sigma_a, sigma_d, sigma_ad=params
	return fromR5([0,0,sigma_a, sigma_d, sigma_ad])

def CIfromR3(params):
	sigma_a, sigma_d, sigma_ad=params
	return CIfromR5([0,0,sigma_a, sigma_d, sigma_ad])


#def fromR5:
#	return

def print_V(params):
	mua, mud, sigma_a, sigma_d, sigma_ad, sigma_e=params
	#mud, mua, sigma_a, sigma_d, sigma_e=params
	global S1, S2, S3, S4, S5, dI, Y, Z, n

	try:
		p=2./(1.+numpy.exp(-sigma_ad) )-1.
		sa=numpy.exp(sigma_a) 
		sd=numpy.exp(sigma_d) 
		se=numpy.exp(sigma_e) 
	except:
		print "Whoops!"
		return sys.float_info.max

	V= (numpy.sqrt(sa*sd)*p)*S3+(sa)*S4+(sd)*S5+se*dI
	print V


def LL(params):
	mua, mud, sigma_a, sigma_d, sigma_ad=params
	s1=time.time()
	global SS11, SS12, SS22, S3, S4, S5, dI, Y, Z, n, l

	try:
		p=2./(1.+numpy.exp(-sigma_ad) )-1.
		sa=numpy.exp(sigma_a) 
		sd=numpy.exp(sigma_d) 
	except:
		print "Whoops!"
		return sys.float_info.max

	if max(abs(mua), abs(mud) )>10:
		return sys.float_info.max
	sY=(Y-mua*S1-mud*S2)
	
	s2=time.time()
	V  = ( (numpy.sqrt(sa*sd)*p+l*mua*mud/(l-1) )*S3-mua*mud/(l-1)*( SS12 )+(sa+l*mua**2/(l-1) )*S4-mua**2/(l-1)*SS11+(sd+mud**2*l/(l-1) )*S5-mud**2/(l-1)*SS22 )
	s3=time.time()

	se=numpy.var(sY)-numpy.trace(V)/n
	V=V+se*dI

	s4=time.time()
	try:
		DEN=scipy.linalg.eigh(V, eigvals=(0,3), eigvals_only=True) #(abs(numpy.sum(numpy.diag(V))))
	except:
		return sys.float_info.max
	s5=time.time()
		
	DEN=(DEN[1]+DEN[0])

	if (DEN<=0):
		return sys.float_info.max
		
	sY=sY/numpy.sqrt(numpy.abs(DEN))
	V=V/DEN

	s6=time.time()
	try:
		u, s, v = numpy.linalg.svd(V)
	except:
		print "Can't get eigen values..."
		return sys.float_info.max
	s7=time.time()

	k=n-1
	K=k

	giV=numpy.transpose(u[0:n,0:K]*numpy.transpose(numpy.transpose(v[0:K, 0:n])/s[0:K]))

	lnD=numpy.sum(numpy.log(2*numpy.pi*s[0:k]) )+numpy.log(abs(DEN))*k


	ret =-(-lnD/2-(numpy.matmul(numpy.matmul(numpy.transpose(sY), giV), sY)[0,0] )/2 )
	s8=time.time()
	print s8-s7, s7-s6, s6-s5, s5-s4, s4-s3, s3-s2, s2-s1
	print "ret:", numpy.real(ret), "vals:", mua*l, mud*l, numpy.sqrt(sa)*l, numpy.sqrt(sd)*l, p, se
	if (not(numpy.isfinite(ret))):
		print "Whoops!"
		return sys.float_info.max
	return numpy.real(ret)

def LL2(params):
	sigma_a, sigma_d, sigma_ad=params
	return LL([0., 0., sigma_a, sigma_d, sigma_ad])


def LLE(params):
	sigma_e=params
	return LL([0., 0., -100., -100., 0., sigma_e])

class LLS:
	def __init__ (self, params):
		self.params=params
	def LL (self, params):
		return LL(list(params)+self.params)



Y=X-numpy.mean(X)
Y=Y/numpy.std(Y)
Z=numpy.zeros((n,1), dtype=numpy.float )
v=numpy.std(Y)

#getOIM([0,0,-7,-7,0, -7])

res2 = minimize(LL2, [-11,-11,0], method='powell')
res3 = minimize(LL, [0,0,-11,-11,0], method='powell' )

print res2
a,b,c,d,e=fromR3(res2.x)
print 0/c, 0/c,c/c,d/c,e
print 0,0,c,d,e
CIfromR3(res2.x)
print res3
a,b,c,d,e=fromR5(res3.x)
print a/c,b/c,c/c,d/c,e
print a,b,c,d,e
CIfromR5(res3.x)

#print_V(res3.x)

#X=res3.x
#X[4]=2.94443
#print LL(X)
#X=res3.x
#X[4]=0.51082
#print LL(X)
#X=res3.x
#X[4]=0.00000
#print LL(X)
#X=res3.x
#X[4]=-0.51082
#print LL(X)
#X[4]=-2.9443
#print LL(X)


#res1 = minimize(LL, [0, 0, v/numpy.sum(numpy.diag(S4))/3., v/numpy.sum(numpy.diag(S5))/3., v/3.], method='SLSQP', tol=10**-9)
#print dI
#print res1.x[0]/numpy.sum(numpy.diag(S4))
#print res1.x[1]/numpy.sum(numpy.diag(S5))
#print res1.x[2]/numpy.sum(numpy.diag(dI))

