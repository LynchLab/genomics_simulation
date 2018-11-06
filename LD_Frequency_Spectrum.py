from scipy.stats import binom
import sys
import scipy
import math
import numpy

N=int(sys.argv[1])
MU=float(sys.argv[2])
RHO=float(sys.argv[3])
F=float(sys.argv[4])
TIME=int(sys.argv[5])

n=N			# population size
mu=MU			# mutation rate
theta=2.*mu*n		# heterozygosity
r=RHO			# recombintion rate
f=F			# inbreeding

# Declaring the probability vectors has 2*n rows and 1 column.
	
# Stochastic models typically use row vectors to represent probabilities, 
# but I will be using column vectors, so everything is transposed.
				
# The probability that a  randomly selected pair of alleles in the focal  
# generation have ? is p[x,y,z]:
p=numpy.arange( float(8.*n*n*n) ).reshape( (2*n, 2*n, 2*n) )	

# The probability that a randomly selected pair of allele in the 
# generation  preceding the focal generation has ? is p_prime[x,y,z]:
p_prime=numpy.arange( float(8.*n*n*n) ).reshape( (2*n, 2*n, 2*n) ) 

h_n=[0]
ln_fact=[0]
binom_pmf={}

#this gives the n^th harmonic number, as defined on wikipedia at https://en.wikipedia.org/wiki/Harmonic_number
def H(n):
	global h_n
	try :
		return h_n[n-1]
	except:
		s=len(h_n)
		for i in range(s, n+1):
			h_n.append(h_n[i-1]+1./float(i) )
	return h_n[n-1]

#defined a lookup binomial pmf for faster iteration
def this_binom (s, n, p):
	key=str(s)+","+str(n)+","+str(p)
	try :
		return binom_pmf[key]
	except:

		binom_pmf[key]=binom.pmf(s, n, p)
		return binom_pmf[key]

#used for the binomial pmf lookup
def lnfact(x):
	global ln_fact
	try:
		return ln_fact[x]

	except:
		s=len(ln_fact)
		for i in range(s, x+1):
			ln_fact.append(ln_fact[i-1]+math.log(i) )
		
	return ln_fact[x]
	
#used for the binomial pmf lookup
def lnchoose(n,k):
	return lnfact(n)-lnfact(k)-lnfact(n-k)

#used for the binomial pmf lookup
def lookup_binom_pmf(s, n, p):
	f=n-s
								#if (x0+x2 != 0):
	lnp=lnchoose(n,s)+math.log(1.-p)*f+math.log(p)*s
	return math.exp(lnp)	


#initializing p_prime to the equations shown in slide ?

print "t, E[p], E[p'], E[p]/( (2*n)*4*n*mu )"

R=numpy.arange( float(2**6) ).reshape( (2, 2, 2, 2, 2, 2) )
T=numpy.arange( float( (2*n)**6) ).reshape( (2*n, 2*n, 2*n, 2*n, 2*n, 2*n) )
TX=numpy.arange( float( (2*n)**3) ).reshape( (2*n, 2*n, 2*n) )
TXI=numpy.arange( float( (2*n)**3) ).reshape( (2*n, 2*n, 2*n) )

for A in range (0, 2):
	for B in range (0, 2):
		for A1 in range (0, 2):
			for B1 in range (0, 2):
				for A2 in range (0, 2):
					for B2 in range (0, 2):
						if A1==A2:
							if B1==B2:
								if (A==A1 and B==B1):
									R[A,B,A1,B1,A2,B2] = 1.0  #this_binom(0, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)+this_binom(2*n, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)
								else:
									R[A,B,A1,B1,A2,B2] = 0.0  #this_binom(0, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)+this_binom(2*n, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)
							else:
								if (A==A1):
									R[A,B,A1,B1,A2,B2] = 0.5  #this_binom(0, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)+this_binom(2*n, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)
								else:
									R[A,B,A1,B1,A2,B2] = 0.0  #this_binom(0, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)+this_binom(2*n, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)
						else:
							if B1==B2:
								if (B==B1):
									R[A,B,A1,B1,A2,B2] = 0.5  #this_binom(0, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)+this_binom(2*n, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)
								else:
									R[A,B,A1,B1,A2,B2] = 0.0  #this_binom(0, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)+this_binom(2*n, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)
							else:
								if ( (A==A1 and B==B1) or (A==A2 and B==B2) ) :
									R[A,B,A1,B1,A2,B2] = (1.-r)/2. #this_binom(0, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)+this_binom(2*n, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)
								else:
									R[A,B,A1,B1,A2,B2] = r/2.    #this_binom(0, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)+this_binom(2*n, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)

G=numpy.arange(float(2**4) ).reshape(2,2,2,2)
H=numpy.arange(float(2**2) ).reshape(2,2)

def multinomial_pmf (X, N, Y):
	tX=[]
	tY=[]
	Y=[Y[0,0], Y[0,1], Y[1,0], Y[1,1] ]
	for x in range(0, 4):
		if Y[x] <= 0:
			if Y[x] == 0:
				if X[x] != 0:
					return 0.0
			else:
				return 0.0
		else:
			tX.append(X[x])
			tY.append(Y[x])
		if X[x] < 0:
			return 0.0
#	x1, x2, x3, x4=X
#	p1, p2, p3, p4=Y
#	return math.exp(lnfact(int(N))-lnfact(x1)-lnfact(x2)-lnfact(x3)-lnfact(x4)+x1*math.log(p1)+x2*math.log(p2)+x3*math.log(p3)+x4*math.log(p4) )
	return scipy.stats.multinomial.pmf(tX, N, tY) 
	

for x0 in range (0, 2*n):
	for x1 in range (0, 2*n):
		for x2 in range (0, 2*n):
			if x0==0 and x1==2 and x2 == 2 :
				p_prime[x0, x1, x2] = 1.0
			else :
				p_prime[x0, x1, x2] = 0.0


for x0_p in range (0, 2*n):
	for x1_p in range (0, 2*n):
		print x0_p, x1_p
		for x2_p in range (0, 2*n):
			x3_p=int(2.*n-x0_p-x1_p-x2_p)

			#simple mutation model where 
			#        1
			#   x0 <---> x1
			#   ^  \   /  ^	
			#   |   \ /   |
			# 2 | 3  X  4 | 5
			#   |   / \   | 
			#   v  /   \  v
			#   x2 <---> x3
			#        6
			#

			x0_r=float(x0_p)*(1-mu)**2/(2.*n)+(x2_p+x1_p)*mu*(1-mu)/(2.*n)+x3_p*mu**2/(2.*n)
			x1_r=float(x1_p)*(1-mu)**2/(2.*n)+(x0_p+x3_p)*mu*(1-mu)/(2.*n)+x2_p*mu**2/(2.*n)
			x2_r=float(x2_p)*(1-mu)**2/(2.*n)+(x0_p+x3_p)*mu*(1-mu)/(2.*n)+x1_p*mu**2/(2.*n)
			x3_r=float(x3_p)*(1-mu)**2/(2.*n)+(x2_p+x1_p)*mu*(1-mu)/(2.*n)+x0_p*mu**2/(2.*n)

			#x0_r=float(x0_p)
			#x1_r=float(x1_p)*(1-mu)/(2.*n)+(x0_p+x3_p)*mu/(2.*n)
			#x2_r=float(x2_p)*(1-mu)/(2.*n)+(x0_p+x3_p)*mu/(2.*n)
			#x3_r=float(x3_p)*(1-mu)/(2.*n)+(x2_p+x1_p)*mu/(2.*n)

			f=F

			G[0,0,0,0] = x0_r*x0_r*(1.-f)+f*x0_r
			G[1,0,0,0] = x2_r*x0_r*(1.-f)
			G[0,1,0,0] = x1_r*x0_r*(1.-f)
			G[1,1,0,0] = x3_r*x0_r*(1.-f)

			G[0,0,1,0] = x0_r*x2_r*(1.-f)
			G[1,0,1,0] = x2_r*x2_r*(1.-f)+f*x2_r
			G[0,1,1,0] = x1_r*x2_r*(1.-f)
			G[1,1,1,0] = x3_r*x2_r*(1.-f)

			G[0,0,0,1] = x0_r*x1_r*(1.-f)
			G[1,0,0,1] = x2_r*x1_r*(1.-f)
			G[0,1,0,1] = x1_r*x1_r*(1.-f)+f*x1_r
			G[1,1,0,1] = x3_r*x1_r*(1.-f)

			G[0,0,1,1] = x0_r*x3_r*(1.-f)
			G[1,0,1,1] = x2_r*x3_r*(1.-f)
			G[0,1,1,1] = x1_r*x3_r*(1.-f)
			G[1,1,1,1] = x3_r*x3_r*(1.-f)+f*x3_r
		
			H=numpy.tensordot(R,G,4)

			if (x0_p==8 and x1_p==0 and x2_p==0):
				print x0_p, x1_p, x2_p, x3_p, H[0,0]*2*n, H[0,1]*2*n, H[1,0]*2*n, H[1,1]*2*n
			if (x0_p==4 and x1_p==4 and x2_p==4):
				print x0_p, x1_p, x2_p, x3_p, H[0,0]*2*n, H[0,1]*2*n, H[1,0]*2*n, H[1,1]*2*n
			if (x0_p==4 and x1_p==4 and x2_p==0):
				print x0_p, x1_p, x2_p, x3_p, H[0,0]*2*n, H[0,1]*2*n, H[1,0]*2*n, H[1,1]*2*n

			for x0 in range (0, 2*n):
				for x1 in range (0, 2*n):
					for x2 in range (0, 2*n):
						x3=int(2.*n-x0-x1-x2)

						#if (x3_p >= 0):
						#	print x0_p/(2.*n), x1_p/(2.*n), x2_p/(2.*n), x3_p/(2.*n)
						#	print x0_r, x1_r, x2_r, x3_r
						#	print H[0,0], H[0,1], H[1,0], H[1,1]

						#if x3 >= 0 and x3_p >=0:
						#	if x0==x0_p and x1==x1_p and x2==x2_p:
						#		T[x0,x1,x2,x0_p,x1_p,x2_p]=1.
						#	else:
						#		T[x0,x1,x2,x0_p,x1_p,x2_p]=0
						#else:
						#	T[x0,x1,x2,x0_p,x1_p,x2_p]=0

						if ( x3_p >= 0 and x3 >= 0 ):
							if ( x3 == 2*n ):
								T[x0,x1,x2,x0_p,x1_p,x2_p] = multinomial_pmf( [0, 0, 0, x3], 2.*n, H ) +  multinomial_pmf( [0,0,x3,0], 2.*n, H ) + multinomial_pmf( [0, x3, 0, 0], 2.*n, H ) +  multinomial_pmf( [x3,0,0,0], 2.*n, H )
							elif ( x3 + x1 == 0 ):
								T[x0,x1,x2,x0_p,x1_p,x2_p] = 0
							elif ( x3 + x2 == 0 ):
								T[x0,x1,x2,x0_p,x1_p,x2_p] = 0
							elif ( x3 + x1 == 2*n ):
								T[x0,x1,x2,x0_p,x1_p,x2_p] = multinomial_pmf( [x1,  0, x3,  0], 2.*n, H ) + multinomial_pmf( [  0, x1,  0,  x3], 2.*n, H ) 
							elif ( x3 + x2 == 2*n ):
								T[x0,x1,x2,x0_p,x1_p,x2_p] = multinomial_pmf( [x2, x3,  0,  0], 2.*n, H ) + multinomial_pmf( [  0,  0, x2,  x3], 2.*n, H )
							else:
								T[x0,x1,x2,x0_p,x1_p,x2_p] = multinomial_pmf( [x0, x1, x2, x3], 2.*n, H )
						else:
							T[x0,x1,x2,x0_p,x1_p,x2_p] = 0

#for x in range(0, 2*n):
#	for y in range(0, 2*n):
#		for z in range(0, 2*n):
#			TX[x,y,z]=T[10,0,0,x,y,z]
#			TXI[x,y,z]=T[6,4,4,x,y,z]

print "----====0====----"

for t in range(0, TIME):
	p=numpy.tensordot(T, p_prime, 3)
	print t, numpy.sum(numpy.absolute(p-p_prime))
#	for x in range(0, 2*n):
#		p[0,x,0]=0
#		p[0,0,x]=0
	p_prime=p[:]/numpy.sum(p)
		
#print  p_prime
#print "---===1===---"
#for x in range (0, 2*n):
#	for y in range (0, 2*n) :
#		for z in range (0, 2*n):
#			for X in range (0, 2*n):
#				for Y in range (0, 2*n) :
#					for Z in range (0, 2*n) :
#						if(T[x,y,z,X,Y,Z]>0):
#							print x, y, z, X, Y, Z, T[x,y,z,X,Y,Z]
#print "---===2===---"
#for x in range (0, 2*n):
#	for y in range (0, 2*n) :
#		for z in range (0, 2*n):
#			if(p[x,y,z]>0):
#				print x,y,z,p[x,y,z]


#if x0==1 and x1==2 and x2 == 2 :

D_hat=numpy.arange( float( (2.*n)**2) ).reshape( (2*n, 2*n) )
r_hat=numpy.arange( float( (2.*n)**2) ).reshape( (2*n, 2*n) )
D2_hat=numpy.arange( float( (2.*n)**2) ).reshape( (2*n, 2*n) )
r2_hat=numpy.arange( float( (2.*n)**2) ).reshape( (2*n, 2*n) )
p_xy =numpy.arange( float( (2.*n)**2) ).reshape( (2*n, 2*n) )
p_x  =numpy.arange( float( (2.*n) ) ).reshape( (2*n ) )

NORM=0

for x in range(0, 2*n):
	p_x[x]=0
	for y in range(0, 2*n):
		D2_hat[x,y] = 0
		r2_hat[x,y] = 0
		D_hat[x,y] = 0
		r_hat[x,y] = 0
		p_xy[x,y]  = 0
		for z in range(0, min(x+1,y+1) ):
			X=x-z
			Y=y-z
			Z=z
			x0=Z/(2.*n)
			x1=X/(2.*n)
			x2=Y/(2.*n)
			x3=1-x0-x1-x2
			p1=float(x0+x1)
			p2=float(x0+x2)
			if (p1*p2*(1-p1)*(1-p2) != 0 and p1 < 1.0 and p2 < 1.0):
				r_hat[x,y]+= p[Z,X,Y]*(x0*x3-x1*x2)/math.sqrt(p1*(1.-p1)*p2*(1-p2) )
				r2_hat[x,y]+= p[Z,X,Y]*( (x0*x3-x1*x2)/math.sqrt(p1*(1.-p1)*p2*(1-p2) ) )**2
				D_hat[x,y]+= p[Z,X,Y]*(x0*x3-x1*x2)
				D2_hat[x,y]+= p[Z,X,Y]*(x0*x3-x1*x2)**2
				print x0, x1, x2, x3, p1, p2, (x0*x3-x1*x2)/math.sqrt(p1*(1.-p1)*p2*(1-p2) ), p[Z,X,Y]
				NORM+=p[Z,X,Y]
			p_xy[x,y] += p[Z,X,Y]
			p_x[x]    += p[Z,X,Y]

print "----====3====----"

print "x, y, P[A=x,B=y], E[D | A=x, B=y], V[D | A=x,  B=y], E[r | A=x, B=y], V[r | A=x, B=y]", numpy.sum(r_hat)/NORM, 
for x in range(0, 2*n):
	for y in range(0, 2*n):
		print x, y, p_xy[x,y], D_hat[x,y]/p_xy[x,y], (D2_hat[x,y]-D_hat[x,y]**2)/p_xy[x,y], r_hat[x,y]/p_xy[x,y],  (r2_hat[x,y]-r_hat[x,y]**2)/p_xy[x,y]

print "----====4====----"

for x in range(0, 2*n):
	print x, p_x[x]

print NORM
print "E[r]=", numpy.sum(r_hat)/NORM, 
print "sqrt(V[r])=", math.sqrt(numpy.sum(r2_hat)/NORM-numpy.sum(r_hat)**2/NORM)
print "E[D]=", numpy.sum(D_hat)/NORM,
print "sqrt(V[D])=", math.sqrt(numpy.sum(D2_hat)/NORM-numpy.sum(D_hat)**2/NORM)

quit()
#print R

for t in range(0, TIME):

	p=numpy.matmul(T, p_prime)

	E_p=0
	E_p_prime=0

	for x in range (0, 2*n):
		E_p+=p[x,0]*x
		E_p_prime+=p_prime[x,0]*x

	p_prime=p[:]

	SUM=numpy.sum(p_prime)
	for y in range (0, 2*n):
		p_prime[y,0]=p_prime[y,0]/SUM

	print t, E_p, E_p_prime, E_p/(8*n*n*mu)


pmp_prime=0
cdf=[0]*(2*n+1)		
c=[0]*(2*n)
SUM=0

for x in range (0, 2*n):
	cdf[x+1]=cdf[x]+p_prime[x,0]
	for y in range (0, 2*n):
		if(x!=0):
			pmp_prime+=(x-y)*T[x,y]*p_prime[y,0]
			SUM+=T[x,y]*p_prime[y,0]
			c[x]+=(x-y)*T[x,y]*p_prime[y,0]

for x in range (0, 2*n):
	c[x]=c[x]/SUM*2*n


pmp_prime=pmp_prime/SUM

abs_delta_p=0
abs_delta_p_i=0

print "x, P(p=x), P( p=x )-P( p'=x ), E [ (x-y)*P( p=x & p'=y ) for x ], P( p<x | p != 0 )"

for x in range (0, 2*n):
	print x, p[x,0], p[x,0]-p_prime[x,0], c[x], 
	if (x!=0):
		print (cdf[x]-cdf[1])/(1.-cdf[1]) 
	else :
		print "NaN"
	abs_delta_p+=abs(p[x,0]-p_prime[x,0] )
	abs_delta_p_i+=abs(p[x,0]-p_prime_i[x,0] )

#print SUM, (1-l0), (1-cdf[1])

print "Initially the value of Sum | p_x-p'_x | is ", abs_delta_p_i, ", after forward simulation this becomes ", abs_delta_p,"(should be ~0)" 

print "Total probability should sum to 1, it sums to ", numpy.sum(p)
print "Total: ", cdf[2*n]

print "l0 =", l0
print "E[ p-p' | 0<p<2*n ] =", pmp_prime, "(Method 1 =", n*mu/(1.-cdf[1])/16.-n*mu/16., ", Method 2 =", c[n]/(2*n)*SUM*(2*n+1.5)/(1.-p[0,0]), ")"
print "E[p] =", E_p, "(should be ~", (2*n)*4*n*mu, ")"
print "please be aware that it takes some time for p to reach equilibrium, and |p_x-p'_x| does not appear to be a good indicator of equilibrium."
