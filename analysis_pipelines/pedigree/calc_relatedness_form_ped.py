#THIS WILL ONLY WORK FOR VERY SMALL PEDIGREES!!!

import math
import sympy 
import numpy

mu=sympy.Symbol("mu")
#P=[[1-mu/2, mu/2], [0.5, 0.5], [mu/2, 1-mu/2] ]
P=[ [1, 0], [0.5, 0.5], [0, 1] ]

def f0(pa,pb,i,j,k,l):
	s=[0,0,0,0]
	m=[i,j,k,l]
	sign=(-1)**m.count(-1)
	for n in 0, 1:
		if m[n]<0:
			s[n]=1.
		elif m[n]==0:
			s[n]=1.-pb
		else :
			s[n]=pb
	for n in 2, 3:
		if m[n]<0:
			s[n]=1.
		elif m[n]==0:
			s[n]=1.-pa
		else :
			s[n]=pa
#	print sign
	return sign*s[0]*s[1]*s[2]*s[3]

def f2(pa, pb, a, b, c, d):
	ret=[[0,0,0],[0,0,0],[0,0,0]]
	for i in 1,2:
		for j in 1,2:
			for k in 1,2:
				for l in 1,2:
					ret[i+j-2][k+l-2]+=(f0(pa,pb, min(-a*i,i-1), min(-b*j,j-1), min(-c*k,k-1), min(-d*l, l-1) ) )
	#Do[Do[Do[Do[AppendTo[ret, f0[pa,pb, Min[-a*i,i-1],Min[-b*j,j-1],Min[-c*k,k-1], Min[-d*l, l-1] ] ], {i,2}], {j,2}], {k,2} ], {l,2} ]; 
	ret=numpy.matrix(ret)
#	print ret.shape
	return ret


#def contract:
#	[X3][Xn]+=t
		
	
O=[0,0,0]

g1=0
g2=0
P1=0
P1s=0
P2=0
P2s=0
P12=0

#let T be a 
#T=

#contract(T, 
pa=sympy.Symbol("p1")
qa=1-pa

pb=sympy.Symbol("p2")
qb=1-pb

#p2=sympy.Symbol("p2")
#q2=1-p

s2a=pa*qa**2+qa*(-pa)**2
s3a=pa*qa**3+qa*(-pa)**3
s4a=pa*qa**4+qa*(-pa)**4
 
s2b=pb*qb**2+qb*(-pb)**2
s3b=pb*qb**3+qb*(-pb)**3
s4b=pb*qb**4+qb*(-pb)**4

fa=sympy.Symbol("f1")
fb=sympy.Symbol("f2")
tab=sympy.Symbol("t")
Gab=sympy.Symbol("Gab")
Gba=sympy.Symbol("Gba")
Dab=sympy.Symbol("Dab")
dab=sympy.Symbol("dab")

fbdev=fb*s2b*(f2(pb,pb, 1, 1, -1, -1))
fadev=fa*s2a*(f2(pa,pa, -1, -1, 1, 1))
tdev=tab*sympy.sqrt(s2a*s2b)*(f2(pa,pb, 1, -1, 1, -1)+f2(pa,pb, 1, -1, -1, 1)+f2(pa,pb, -1, 1, 1, -1)+f2(pa,pb, -1, 1, -1, 1) )/2
gabdev=Gab*(s3a*s3a*s3b)**(1/3)*(f2(pa,pb, 1, 1, 1, -1)+f2(pa,pb, 1, 1, -1, 1) )
gbadev=Gba*(s3b*s3b*s3a)**(1/3)*(f2(pa,pb, -1, 1, 1, 1)+f2(pa,pb, 1, -1, 1, 1) )
ddev=dab*sympy.sqrt(s4a*s4b)*(f2(pa,pb, 1, 1, 1, 1) )
Ddev=Dab*s2a*s2b*(f2(pa,pb, 1, 1, 1, 1) )
 
Q=f2(pa, pb, -1,-1,-1,-1)+fadev+fbdev+tdev+gabdev+gbadev+ddev+Ddev
print Ddev
#Q=f2(pa, pb, -1,-1,-1,-1)+Ddev+ddev

Qs=0

#Q=numpy.matrix([[qa**2+pa*qa*fa,0,0],[0,2*(1-fa)*pa*qa,0],[0,0,pa**2+pa*qa*fa]])
for i in range(0, 3):
	for j in range(0, 3):
		g1+=P[i][1]*Q[i,j]
		g2+=P[j][1]*Q[i,j]
		Qs+=Q[i,j]
		P1+=i*Q[i,j]/2.
		P1s+=(i/2.)**2*Q[i,j]
		P2+=j*Q[i,j]/2.
		P2s+=(j/2.)**2*Q[i,j]
		P12+=(i/2.)*(j/2.)*Q[i,j]
		for k in range (0, 2):
			for l in range (0, 2):
				O[k+l]+=P[i][k]*P[j][l]*Q[i,j]

#print sympy.simplify(O)
#print sympy.simplify(g1)
#print sympy.simplify(g2)
print "1:"
print sympy.simplify( Qs )
print "2:"
print sympy.simplify( (O[2]-g1*g2)/sympy.sqrt(g1*(1-g1)*g2*(1-g2) ) )
print "3:"
print sympy.simplify( (P12-P1*P2)/sympy.sqrt( (P1s-P1**2)*(P2s-P2**2) ) )
