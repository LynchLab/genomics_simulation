from scipy.stats import binom
import sys
import math
import numpy

N=int(sys.argv[1])
MU=float(sys.argv[2])
TIME=int(sys.argv[3])

n=N			# population size
mu=MU			# mutation rate
theta=2.*mu*n		# heterozygosity

# Declaring the probability vectors has 2*n rows and 1 column.
	
# Stochastic models typically use row vectors to represent probabilities, 
# but I will be using column vectors, so everything is transposed.
				
# The probability that a  randomly selected allele in the focal  generation 
# has a count of x is p[x]:
p=numpy.arange( float(2.*n) ).reshape( (2*n, 1) )	

# The probability that a  randomly selected allele in the  generation 
# preceding the focal generation has a count of x is p_prime[x]:
p_prime=numpy.arange( float(2.*n) ).reshape( (2*n, 1) ) 

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
	lnp=lnchoose(n,s)+math.log(1.-p)*f+math.log(p)*s
	return math.exp(lnp)	


#initializing p_prime to the equations shown in slide ?

p_prime[0,0]=1-4.*n*theta*(H(2*n-1)/(n*2.-2.) )

print p_prime[0,0]
print "D2=", n
print "G1=", theta
print "G2=", 1./H(2*n-1)
print 4.*float(n)*theta/(H(2*n-1)*(n*2-2) )
l0=1.-4.*float(n)*theta/(H(2*n-1)*(n*2-2) )

if (p_prime[0,0] < 0.0):
	p_prime[0,0]=0.05

for y in range (1, 2*n):
	p_prime[y,0]=(1-p_prime[0,0])/( H(2*n-1)*y )

SUM=numpy.sum(p_prime)
for y in range (0, 2*n):
	p_prime[y,0]=p_prime[y,0]/SUM

p_prime_i=p_prime[:]

#running the simulation out to time TIME

print "t, E[p], E[p'], E[p]/( (2*n)*4*n*mu )"

#T=numpy.matrix(ncol=2*n, nrow=2*n)
T=numpy.arange( float(4.*n*n) ).reshape( (2*n, 2*n) )

for x in range (0, 1):
	for y in range (0,2*n):
		T[x,y]=this_binom(0, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)+this_binom(2*n, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)
#		p[0]+=this_binom(0, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)*p_prime[y]+this_binom(int(2*n), int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)*p_prime[y]

for x in range (1,2*n):
	for y in range (0,2*n):
		T[x,y]=this_binom(x, int(2*n), float(y)/float(2.*n)*(1.-2.*mu)+mu)


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
#The expectation on slide 11
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
		print (cdf[x]-cdf[1])/(1.-cdf[1]) #,
	else :
		print "NaN"
	abs_delta_p+=abs(p[x,0]-p_prime[x,0] )
	abs_delta_p_i+=abs(p[x,0]-p_prime_i[x,0] )

print SUM, (1-l0), (1-cdf[1])

print "Initially the value of Sum | p_x-p'_x | is ", abs_delta_p_i, ", after forward simulation this becomes ", abs_delta_p,"(should be ~0)" 

print "Total probability should sum to 1, it sums to ", numpy.sum(p)
print "Total: ", cdf[2*n]

print "l0 =", l0
print "E[ p-p' | 0<p<2*n ] =", pmp_prime, "(Method 1 =", (n*mu/(1.-cdf[1])/16.-n*mu/16.)*n/8., ", Method 2 =", c[n]/(2*n)*(2*n+1.5), ")"
print "E[p] =", E_p, "(should be ~", (2*n)*4*n*mu, ")"
print "please be aware that it takes some time for p to reach equilibrium, and |p_x-p'_x| does not appear to be a good indicator of equilibrium."
