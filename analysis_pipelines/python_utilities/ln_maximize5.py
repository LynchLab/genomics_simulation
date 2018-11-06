#!/bin/python
# -*- coding: utf-8 -*-

import sys
import numpy
import scipy
import sympy
import time

from numpy import exp
from numpy import log
from numpy import sqrt
from scipy.optimize import minimize


#PATH=sys.argv[1]
PATH="../analysis_files/"

File=open(PATH+"mapgd_relatedness.out")
X = []
for line in File:
	if line[0]!='@':
		line=map(float, line.split('\t')[1:] )
		X.append(line)
n=len(X)

X=(numpy.asmatrix(X))

print X.shape
S1=numpy.asmatrix(X[0:n, 0:1]).reshape((n,1))
print S1.shape
S2=numpy.asmatrix(X[0:n, 1:2]).reshape((n,1))
print S2.shape
S3=numpy.asmatrix(X[0:n, (2+0*n):(2+1*n)]).reshape((n,n))
print S3.shape
S4=numpy.asmatrix(X[0:n, (2+1*n):(2+2*n)]).reshape((n,n))
print S4.shape
S5=numpy.asmatrix(X[0:n, (2+2*n):(2+3*n)]).reshape((n,n))
print S5.shape


def t(p):
	return numpy.transpose(p)

SS11=S1.dot(t(S1))
SS12=S1.dot(t(S2))+S2.dot(t(S1))
SS22=S2.dot(t(S2))

dI=numpy.full( (n,n), -1./(n-1.) )
numpy.fill_diagonal(dI, 1.)

print(dI)
#dI=numpy.identity(n)

File=open(PATH+"plink.pheno")

X = []

for line in File:
	if line[0]!='@':
		line=map(float, line.split('	') )
		X.append(line[2])

X=numpy.asmatrix(X).reshape((n, 1))


def lm():
	global SS11, SS12, SS22, S1, S2, S3, S4, S5, dI, Y, Z, n	
	A=S4
	B=SS11
	C=S5
	D=SS22
	E=S3
	F=SS12
	G=dI
	r=n**2
	X=S1
	X=S2
	X=numpy.concatenate( (S1, S2), axis=1)
	V=numpy.linalg.inv( t(X)*X )*t(X)*Y
	b=V[0,0]
	d=V[1,0]
	yyt=(Y-b*S1-d*S2)*t(Y-b*S1-d*S2)
	X=numpy.concatenate( (A.reshape(r,1), C.reshape(r,1), E.reshape(r,1),  G.reshape(r,1) ), axis=1)
	V=numpy.linalg.inv( t(X)*X )*t(X)*yyt.reshape(r,1)
	a=max(V[0,0], 10**-6)
	c=max(V[1,0], 10**-6)
	e=V[2,0]
	f=max(V[3,0], 10**-6)
	return numpy.asmatrix(bus([a,b,c,d,e,f])).reshape((6,1))

def opls():
	global SS11, SS12, SS22, S1, S2, S3, S4, S5, dI, Y, Z, n	
	A=S4
	B=SS11
	C=S5
	D=SS22
	E=S3
	F=SS12
	G=dI
	H=[A, C, E, G]

	X=S1
	X=S2
	X=numpy.concatenate( (S1, S2), axis=1)
	V=numpy.linalg.inv( numpy.transpose(X)*X )*numpy.transpose(X)*Y

	b=V[0,0]
	d=V[1,0]

	var=numpy.var(Y-b*S1-d*S2)
	bb=numpy.ones((4,1))
	L=numpy.ones((4,4))
	print var
	for i in range(0, 4):
		bb[i,0]=numpy.trace( var * H[i])
		for j in range(0, 4):
			L[i,j]=numpy.trace( H[i]*H[j])

	print L.shape, bb.shape
	t=numpy.linalg.inv(L).dot(bb)
	print t
	a=t[0,0]
	c=t[1,0]
	e=t[2,0]
	f=t[3,0]
	
	return numpy.asmatrix(bus([a,b,c,d,e,f])).reshape((6,1))
	
	


# sqrt(exp(a)*exp(c) )*( 2/(1+exp(e) ) -1)
def sub(params):
	a, b, c, d, e, f=map(float, params)
	return [exp(a), b, exp(c), d, e, exp(f) ]
#	return [exp(a), b, exp(c), d, sqrt(exp(a+c) )*(2./(1.+exp(e) )-1.), exp(f) ]
def bus(params):
	if len(params)==6:
		a, b, c, d, e, f=map(float, params)
	elif len(params)==1:
		a = params[0,0]
		b = params[0,1]
		c = params[0,2]
		d = params[0,3]
		e = params[0,4]
		f = params[0,5]
	else:
		print len(params)
		print params
	return  [log(a),b,log(c),d,e,log(f)]

def get_Sigma(params):
	a, b, c, d, e, f=sub(params)
	global SS11, SS12, SS22, S1, S2, S3, S4, S5, dI, Y, Z, n	
	A=S4
	B=SS11
	C=S5
	D=SS22
	E=S3
	F=SS12
	G=dI
	Sigma=a*A+c*C+e*E+f*G
	if ( (a+b**2) !=0 ):
		Sigma -= 2*a*b**2/(a+b**2)*B
	if ( (c+d**2) !=0 ):
		Sigma -= 2*c*d**2/(c+d**2)*D
	if ( (e+b*d) !=0 ):
		Sigma -= 2*e*b*d/(e+b*d)*F
	return Sigma

last_ret=[]
last_params=[]
def get_pseudo(Sigma, params):
	global n, last_ret, last_params
	K=n-1
#	print "Try: ", params
#	print params, last_params
#	print params != last_params
#	print params == last_params
#	if( list(params)!=list(last_params) ):
	if( True ):
#		print "calc!"
		try:
			u, s, v = numpy.linalg.svd(Sigma)
		except:
			print("Can't get eigen values...")
			exit(0)
		Sigma_inv = numpy.transpose(u[0:n,0:K]*numpy.transpose(numpy.transpose(v[0:K, 0:n])/s[0:K]))
		last_ret = [Sigma_inv, s]
		last_params = list(params[:])
	else:
#		print "Return last"
		Sigma_inv, s=last_ret
	return Sigma_inv, s

def nabla(params):
	if len(params)==6:
		a, b, c, d, e, f=map(float, params)
	elif len(params)==1:
		a = params[0,0]
		b = params[0,1]
		c = params[0,2]
		d = params[0,3]
		e = params[0,4]
		f = params[0,5]
		print params.shape
	global SS11, SS12, SS22, S1, S2, S3, S4, S5, dI, Y, Z, n
	A=S4
	B=SS11
	C=S5
	D=SS22
	E=S3
	F=SS12
	G=dI
	k=n-1
	K=k
	dSigma, dmu, ddSigma=get_partials(params)
	S=get_Sigma(params)
	P, s=get_pseudo(S, params)
	I=numpy.identity(n)
	J=numpy.zeros( (6) )
	for x in range(0, 6):
		dS=dSigma[x]
		dm=dmu[x]
		dP=(-P*dS*P+P*P*dS*(I-S*P)+(I-P*S)*dS*P*P)
		J[x]=-(-t(Y-b*S1-d*S2)*dP*(Y-b*S1-d*S2)/2.-(2.*numpy.pi*P*dS).trace()/(4*numpy.pi)-t(dm)*P*(Y-b*S1-d*S2)/2.-t(Y-b*S1-d*S2)*P*(dm)/2.)[0,0]
	print J
	return J

def nabla2(params):
	print "hi", params
	if len(params)==6:
		a, b, c, d, e, f=map(float, params)
	elif len(params)==1:
		a = params[0,0]
		b = params[0,1]
		c = params[0,2]
		d = params[0,3]
		e = params[0,4]
		f = params[0,5]
		print params.shape
	global SS11, SS12, SS22, S1, S2, S3, S4, S5, dI, Y, Z, n
#	A=S4
#	B=SS11
#	C=S5
#	D=SS22
#	E=S3
#	F=SS12
#	G=dI
	k=n-1
	K=k
	dSigma, dmu, ddSigma=get_partials(params)
	S=get_Sigma(params)
	P, s=get_pseudo(S, params)
	I=numpy.identity(n)
	H=numpy.zeros( (6,6) )
	PP=P*P
	IPS=I-P*S
	ISP=t(IPS)
	sY=Y-b*S1-d*S2
	dP=[0,0,0,0,0,0]
	PP=P*P
	ISP=I-S*P
	for x in range(0, 6):
		dS=dSigma[x]
		PPdSISP=PP*dS*ISP
		dP[x]=(-P*dS*P)+PPdSISP+t(PPdSISP)
	for x in range(0, 6):
		for y in range (x, 6):
			dmx=dmu[x]
			dmy=dmu[y]
			dPx=(-P*dSigma[x]*P+PP*dSigma[x]*ISP+IPS*dSigma[x]*PP)
			dPy=(-P*dSigma[y]*P+PP*dSigma[y]*ISP+IPS*dSigma[y]*PP)
			ddPxy=get_ddP(dSigma[x], dSigma[y], ddSigma[x][y], dP[x], dP[y], P, S, I)
			a1=t(dmx)*P*(dmy)/2.
			b1=t(dmx)*dPy*(sY)/2.
			c1=t(dmy)*dPx*(sY)/2.
			d1=((dPx*dSigma[y]+P*ddSigma[x][y]).trace()/(2.))
			H[x,y]=(d1+t(sY)*ddPxy*(sY)/2.+a1+t(a1)+b1+t(b1)+c1+t(c1)
					)[0,0] 
			if(x!=y):
				H[y,x]=H[x,y]
	return H

	
def ll(params):
	print params
	if len(params)==6:
		a, b, c, d, e, f=map(float, params)
	elif len(params)==1:
		a = params[0,0]
		b = params[0,1]
		c = params[0,2]
		d = params[0,3]
		e = params[0,4]
		f = params[0,5]
		print params.shape
	global Y, S1, S2, n
	K=n-1
	sY = (Y-b*S1-d*S2)
	try:
		S=get_Sigma(params)
	except:
		return sys.float_info.max

	try:
		P, s = get_pseudo(S, params)
	except: 
		return sys.float_info.max

	return t(sY)*P*sY/2. + numpy.sum(numpy.log(2*numpy.pi*s[0:K]) )/2.

def get_ddP(D1, D2, D12, dX1, dX2, P, S, I):
	PP=P*P #1
	SP=S*P #2
	ISP=I-SP
	D1ISP=D1*ISP #3
	D1P=D1*P
	D2P=D2*P
	dX2D1P=dX2*D1P
	dX2P=dX2*P

	RR= (dX2P*(D1ISP) 
	    +t(dX2P)*(D1ISP) 
	    +PP*D12*ISP    
	    +PP*D1*(-D2P-S*dX2) ) 
	return (
	      -dX2D1P
              -P*D12*P
              -t(dX2D1P)
		+RR
		+t(RR) )




def get_partials(params):
	a, b, c, d, e, f=map(float, params)
	global SS11, SS12, SS22, S1, S2, S3, S4, S5, dI, Y, Z, n
	A=S4
	B=SS11
	C=S5
	D=SS22
	E=S3
	F=SS12
	G=dI

	dSigma=[]
	dmu=[]
	ddSigma=[]
	for x in range (0, 6):
		dSigma.append(numpy.zeros((n,n)) )
		dmu.append(numpy.zeros((n,1)))
		ddSigma.append([])
		for y in range (0, 6):
			ddSigma[x].append(numpy.zeros((n,n)))
#	    0
	dmu[1]=-S1
#	    2
	dmu[3]=-S2

	dSigma[0] = A*exp(a) + B*(-2*b**2*exp(a)/(b**2 + exp(a)) + 2*b**2*exp(2*a)/(b**2 + exp(a))**2)
	dSigma[1] = B*(4*b**3*exp(a)/(b**2 + exp(a))**2 - 4*b*exp(a)/(b**2 + exp(a))) + F*(2*b*d**2*e/(b*d + e)**2 - 2*d*e/(b*d + e))
	dSigma[2] = C*exp(c) + D*(-2*d**2*exp(c)/(d**2 + exp(c)) + 2*d**2*exp(2*c)/(d**2 + exp(c))**2)
	dSigma[3] = D*(4*d**3*exp(c)/(d**2 + exp(c))**2 - 4*d*exp(c)/(d**2 + exp(c))) + F*(2*b**2*d*e/(b*d + e)**2 - 2*b*e/(b*d + e))
	dSigma[4] = E + F*(2*b*d*e/(b*d + e)**2 - 2*b*d/(b*d + e))
	dSigma[5] = G*exp(f)

	ddSigma[0][0] = (A + B*(-2*b**2/(b**2 + exp(a)) + 6*b**2*exp(a)/(b**2 + exp(a))**2 - 4*b**2*exp(2*a)/(b**2 + exp(a))**3))*exp(a)
	ddSigma[0][1] = 4*B*b*(b**2/(b**2 + exp(a)) - 2*b**2*exp(a)/(b**2 + exp(a))**2 - 1 + exp(a)/(b**2 + exp(a)))*exp(a)/(b**2 + exp(a))
	ddSigma[0][2] = 0
	ddSigma[0][3] = 0
	ddSigma[0][4] = 0
	ddSigma[0][5] = 0
	ddSigma[1][1] = 4*B*(-4*b**4*exp(a)/(b**2 + exp(a))**3 + 5*b**2*exp(a)/(b**2 + exp(a))**2 - exp(a)/(b**2 + exp(a))) + 4*F*(-b*d**3*e/(b*d + e)**3 + d**2*e/(b*d + e)**2)
	ddSigma[1][2] = 0
	ddSigma[1][3] = 2*F*e*(-2*b**2*d**2/(b*d + e)**2 + 3*b*d/(b*d + e) - 1)/(b*d + e)
	ddSigma[1][4] = 2*F*d*(-2*b*d*e/(b*d + e)**2 + b*d/(b*d + e) + e/(b*d + e) - 1)/(b*d + e)
	ddSigma[1][5] = 0
	ddSigma[2][2] = (C + D*(-2*d**2/(d**2 + exp(c)) + 6*d**2*exp(c)/(d**2 + exp(c))**2 - 4*d**2*exp(2*c)/(d**2 + exp(c))**3))*exp(c)
	ddSigma[2][3] = 4*D*d*(d**2/(d**2 + exp(c)) - 2*d**2*exp(c)/(d**2 + exp(c))**2 - 1 + exp(c)/(d**2 + exp(c)))*exp(c)/(d**2 + exp(c))
	ddSigma[2][4] = 0
	ddSigma[2][5] = 0
	ddSigma[3][3] = 4*D*(-4*d**4*exp(c)/(d**2 + exp(c))**3 + 5*d**2*exp(c)/(d**2 + exp(c))**2 - exp(c)/(d**2 + exp(c))) + 4*F*(-b**3*d*e/(b*d + e)**3 + b**2*e/(b*d + e)**2)
	ddSigma[3][4] = 2*F*b*(-2*b*d*e/(b*d + e)**2 + b*d/(b*d + e) + e/(b*d + e) - 1)/(b*d + e)
	ddSigma[3][5] = 0
	ddSigma[4][4] = 4*F*b*d*(-e/(b*d + e) + 1)/(b*d + e)**2
	ddSigma[4][5] = 0
	ddSigma[5][5] = G*exp(f)


	return [dSigma, dmu, ddSigma]

def getOIM(params):
	a,b,c,d,e,f=map(float, params)
	global SS11, SS12, SS22, S1, S2, S3, S4, S5, dI, Y, Z, n

	H=numpy.full( (6,6), 0. )
	g=numpy.full( (6,1), 0. )

	#for delta in [10, 1, 0.1, 0.01, 0.001]:#, 0.0000001, 0.000000001]:
	ds=[1., 10., 1., 10., 1., 1.]
	for delta in [1, 0.1, 0.01, 0.001]:
		Theta=[a,b,c,d,e,f]
		ret0 = ll(Theta) 
		for x in range (0,6):
			dxy0=Theta[:]
			dx=ds[x]*delta
			dxy0[x]-= dx
			for y in range (0,6):
				dy=ds[y]*delta
				dxdy=Theta[:]
				dxdy[x]-= dx
				dxdy[y]-= dy
				x0dy=Theta[:]
				x0dy[y]-= dy
				H[x,y]=( ( ll(dxdy)-ll(x0dy) )/dx-( ll(dxy0)-ret0 )/dx )/dy
			g[x,0]=( ll(dxy0)-ret0 )/dx
		print H[0,0], H[1,1], H[2,2], H[3,3], H[4,4], H[5,5], H[1,3], H[1,4], H[3,4]

	V=get_Sigma(params)
	P, s=get_pseudo(V, params)
	print
	I=numpy.identity(n)

	dSigma, dmu, ddSigma=get_partials(params)
	dP=[0,0,0,0,0,0]
	eg=nabla(params)
	for x in range(0, 6):
		dS=dSigma[x]
		dm=dmu[x]
		dP[x]=(-P*dS*P+P*P*dS*(I-V*P)+(I-P*V)*dS*P*P)
		print "dl/dTheta_{"+str(x)+"}", eg[x,0], g[x,0], eg[x,0]/g[x,0]
#(-t(Y-b*S1-d*S2)*dP[x]*(Y-b*S1-d*S2)/2.-(2.*numpy.pi*P*dS).trace()/(4*numpy.pi)-t(dm)*P*(Y-b*S1-d*S2)/2.-t(Y-b*S1-d*S2)*P*(dm)/2.)[0,0], g[x,0]

	eH=numpy.full( (6,6), 0. )
	eH2=numpy.full( (6,6), 0. )
	eH=nabla2(params)
	for x in range(0, 6):
		for y in range (x, 6):
#			t0=time.time()
			dmx=dmu[x]
			dmy=dmu[y]
			dPx=(-P*dSigma[x]*P+P*P*dSigma[x]*(I-V*P)+(I-P*V)*dSigma[x]*P*P)
			dPy=(-P*dSigma[y]*P+P*P*dSigma[y]*(I-V*P)+(I-P*V)*dSigma[y]*P*P)
			ddPxy=get_ddP(dSigma[x], dSigma[y], ddSigma[x][y], dP[x], dP[y], P, V, I)
#			t1=time.time()
#			print t1-t0
			print "d^2 l/dTheta_{"+str(x)+" "+str(y)+"}", ((dPx*dSigma[y]+P*ddSigma[x][y]).trace()/(2.))[0,0], 
			print  (t(Y-b*S1-d*S2)*ddPxy*(Y-b*S1-d*S2)/2.
					+t(dmx)*dPy*(Y-b*S1-d*S2)/2.
					+t(dmx)*P*(dmy)/2.
					+t(dmy)*P*(dmx)/2.
					+t(Y-b*S1-d*S2)*dPy*(dmx)/2.
					+t(dmy)*dPx*(Y-b*S1-d*S2)/2.
					+t(Y-b*S1-d*S2)*dPx*(dmy)/2.
					)[0,0], 
			print  H[x,y]/( (dPx*dSigma[y]+P*ddSigma[x][y]).trace()/(2.)
					+t(Y-b*S1-d*S2)*ddPxy*(Y-b*S1-d*S2)/2.
					+t(dmx)*dPy*(Y-b*S1-d*S2)/2.
					+t(dmx)*P*(dmy)/2.
					+t(dmy)*P*(dmx)/2.
					+t(Y-b*S1-d*S2)*dPy*(dmx)/2.
					+t(dmy)*dPx*(Y-b*S1-d*S2)/2.
					+t(Y-b*S1-d*S2)*dPx*(dmy)/2.
					)[0,0], 
			eH2[x,y]=(t(Y-b*S1-d*S2)*ddPxy*(Y-b*S1-d*S2)/2.+t(dmx)*dPy*(Y-b*S1-d*S2)/2.+t(dmx)*P*(dmy)/2.+t(dmy)*P*(dmx)/2.+t(Y-b*S1-d*S2)*dPy*(dmx)/2.+t(dmy)*dPx*(Y-b*S1-d*S2)/2.+t(Y-b*S1-d*S2)*dPx*(dmy)/2.)[0,0]+((dPx*dSigma[y]+P*ddSigma[x][y]).trace()/(2.))[0,0]
			print H[x,y]
#			print eH[x,y], eH[x,y]/H[x,y], H[x,y], eH2[x,y]
	try:
		print numpy.diag(numpy.linalg.inv(H))
	except:
		print " A thing"
	print numpy.diag(numpy.linalg.inv(eH))

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

theta=lm()
#theta=opls()
print "least-sq fit", sub(theta)
print ll(theta)
#res3 = minimize(ll, theta, method='SLSQP', jac=nabla, tol=1e-8 )
#res3 = minimize(ll, theta, method='SLSQP', jac=nabla )
res3 = minimize(ll, theta, method='BFGS', jac=nabla, hess=nabla2 )
print res3

#theta=numpy.asmatrix(res3.x).reshape(6,1)

for z in range(0,10):
	print z
	h=nabla2(theta)
	j=nabla(theta)
	print "h:", h
	print "t:", theta
	print "j:", j
	print "ll:", ll(theta)
	h=numpy.linalg.inv(h)
	D=numpy.diag(h)
	if abs(max([3, 7, -10], key=abs)) < 10**-5:
		break
	for x in range(0, 6):
		print x, D[x]
	print h.dot(j).reshape(6,1)
	theta_2=theta+h.dot(j).reshape(6,1)
	theta=theta_2

h=nabla2(theta)
h=numpy.linalg.inv(h)
T=numpy.asmatrix([sqrt(abs(h[0,0])), sqrt(abs(h[1,1])), sqrt(abs(h[2,2])), sqrt(abs(h[3,3])), sqrt(abs(h[4,4])), sqrt(abs(h[5,5])) ]).reshape((6,1))

print sub(theta-T)
print sub(theta)
print sub(theta+T)

print "start ", ll(theta)
for x in range(0, 6):
	this=theta.copy()
	this[x]+=T[x]
	print this
	print x, "!!", ll(this) 
quit()

#getOIM(res3.x)

