#!/usr/bin/python 

from math import sqrt

def print_hap(G):
	print "Haplotypes:"
	global L
	H=[0]*2**L
	for x in range(0, 2**(2*L) ):
		H[l(x)]+=G[x]/2.		
		H[r(x)]+=G[x]/2.	
	for x in range(2**L-1, -1, -1):
		print '{0:03b}'.format(x), H[x]

def print_gen(G):
	print "Genotypes:"
	for x in range(2**(2*L)-1, -1, -1):
		print '{0:03b}'.format(l(x)), '{0:03b}'.format(r(x)), G[x]
		

def comp(x):
	global L
	if min(x)<L:
		return [x[0]+L, x[1]+L]
	else :
		return [x[0]-L, x[1]-L]
def l(y):
	global L
	return y >> L
def r(y):
	global L
	m = 2**L-1
	return y & m

def is_set(G, x):
	global L
	return (int)(G & (1 << x) != 0)

def get_r(DP, Cm):
	dist=abs(DP[0]-DP[1])
	if dist==1:
		return Cm[L-min(DP)-1]
	else :
		state1=1.0
		state2=0.0
#		print "min/max: ", min(DP), max(DP)
		for x in range (L-max(DP), L-min(DP)):
#			print "x", x
			t_state1=state1*(1-Cm[x])+state2*Cm[x]
			t_state2=state1*Cm[x]+state2*(1-Cm[x])
			state1=t_state1
			state2=t_state2	
#		print state2
		return state2
	

def get_D(Ga, Gb, DP1):
	global L
	DP2=[DP1[0]+L,DP1[1]+L]
	op=(get_moment(Ga, DP1)+get_moment(Ga, DP2 ) )/2.
	o =(get_moment(Gb, DP1)+get_moment(Gb, DP2 ) )/2.
#	print  op, o
	r=get_r(DP1)
	if r!=0:
		return (o/(r*op)-1./r+1.)
	else :
		return 0

def get_S(Ga, Gb, DT1):
	global L
	DP2=[DT1[0]+L, DT1[1]+L, PT1[2]+L]
	op=(get_moment(Ga, DT1)+get_moment(Ga, DT2 ) )/2.
	o =(get_moment(Gb, DT1)+get_moment(Gb, DT2 ) )/2.
#	print  op, o
	r1=get_r([DT1[0], DT1[1]])
	r2=get_r([DT1[1], DT1[2]])
	r3=get_r([DT1[0], DT1[2]])
#	C
#	G
#	I
	if r!=0:
#		=(C+I-op)/(I-C*G)
		return (o/(r*op)-1./r+1.)
	else :
		return 0

def get_moment(g, loci):
	global L
	order=len(loci)
	if order==1:
		E=0
		for x in range(0, 2**(2*L) ):
			if ( is_set(x, loci[0] ) ):
				E+=g[x]
		return E
	elif order==2:
		E0=get_moment(g, [loci[0]] )
		E1=get_moment(g, [loci[1]] )
		E01=0
	#	print E0, E1
		for x in range(0, 2**(2*L) ):
			E01+=g[x]*(is_set(x, loci[0])-E0)*(is_set(x, loci[1])-E1)
		return E01	
	elif order==3:
		E0=get_moment(g, [loci[0]] )
		E1=get_moment(g, [loci[1]] )
		E2=get_moment(g, [loci[2]] )
		E012=0
		for x in range(0, 2**(2*L) ):
			E012+=g[x]*(is_set(x, loci[0])-E0)*(is_set(x, loci[1])-E1)*(is_set(x, loci[2])-E2)
		return E012

Rtable={}

def R(x,y, Cm):
	global Rlookup

	h1=l(y)
	h2=r(y)

	try:
		return Rtable[x][y]

	except:
		if x not in Rtable.keys():
			Rtable[x]={}

	state1=0.5
	state2=0.5

	for s in range (0, L):
		m = 2**(L-s-1)
		t_state1=state1*(1-Cm[s])+state2*Cm[s]
		t_state2=state1*Cm[s]+state2*(1-Cm[s])
		state1=t_state1
		state2=t_state2

		if( (h1 & m ) != (x & m )  ):
			state1=0
		if( (h2 & m ) != (x & m )  ):
			state2=0
	Rtable[x][y]=state1+state2
	return Rtable[x][y]
	
		
L=3

for j in range(1, 33):
	
	Cm=[ j/(32.*2.), j/(32.*2.), j/(32.*2.)]
	Rtable={}
	#Cm=[0.005, 0.25, 0.25]
	#Cm=[0.5, 0.5, 0.5]
	
	G0=[0]*2**(2*L)
	G1=[0]*2**(2*L)
	G2=[0]*2**(2*L)
	G3=[0]*2**(2*L)
	G2q=[0]*2**(2*L)
	G3q=[0]*2**(2*L)
	G1p=[0]*2**(2*L)
	G2p=[0]*2**(2*L)
	G3p=[0]*2**(2*L)
	
	H1=[0]*2**(L)
	
	A=0.5
	B=0.5
	C=0.5
	
	a=1-A
	b=1-B
	c=1-C
	
	r_ab=0.25
	r_ac=0.25
	r_bc=0.25
	s_abc=0.0
	
	H0=    [a*b*c+( r_ab+r_ac+r_bc)/2.-s_abc,
		a*b*C+( r_ab-r_ac-r_bc)/2.+s_abc,
		a*B*c+(-r_ab+r_ac-r_bc)/2.+s_abc,
		a*B*C+(-r_ab-r_ac+r_bc)/2.-s_abc,
		A*b*c+(-r_ab-r_ac+r_bc)/2.+s_abc,
		A*b*C+(-r_ab+r_ac-r_bc)/2.-s_abc,
		A*B*c+( r_ab-r_ac-r_bc)/2.-s_abc,
		A*B*C+( r_ab+r_ac+r_bc)/2.+s_abc]
	
	#H0=    [a*b+( r_ab),
	#	a*B+(-r_ab),
	#	A*b+(-r_ab),
	#	A*B+( r_ab)]
	
	
	#for x in range (0, 2**(L) ):
	#	print '{0:03b}'.format(x), H0[x]
	
	#print 
	
	for x in range (0, 2**(2*L) ):
		G0[x]=H0[r(x)]*H0[l(x)]
	#	print '{0:03b}'.format(l(x) ), '{0:03b}'.format(r(x)), G0[x]
	
		
	#for x in range(0, 2**(2*L) ):
	#	for y in range(0, 2**(2*L) ):
	#		for z in range(0, 2**(2*L) ):
	#			G1[x]+=R(l(x),y, Cm)*R(r(x),z, Cm)*G0[y]*G0[z]
	#G0[:]=G1[:]
	#G1=[0]*2**(2*L)

	#for x in range(0, 2**(2*L) ):
	#	for y in range(0, 2**(2*L) ):
	#		for z in range(0, 2**(2*L) ):
	#			G1[x]+=R(l(x),y, Cm)*R(r(x),z, Cm)*G0[y]*G0[z]
	#G0[:]=G1[:]
	#G1=[0]*2**(2*L)

	#for x in range(0, 2**(2*L) ):
	#	for y in range(0, 2**(2*L) ):
	#		for z in range(0, 2**(2*L) ):
	#			G1[x]+=R(l(x),y, Cm)*R(r(x),z, Cm)*G0[y]*G0[z]
	#G0[:]=G1[:]
	#G1=[0]*2**(2*L)
	
	for x in range(0, 2**(2*L) ):
		for y in range(0, 2**(2*L) ):
			G1p[x]+=R(l(x),y, Cm)*R(r(x),y, Cm)*G0[y]
			for z in range(0, 2**(2*L) ):
				G1[x]+=R(l(x),y, Cm)*R(r(x),z, Cm)*G0[y]*G0[z]
	
	G2vw=[0]*2**(2*L)
	
	#for v in range(0, 2**(2*L) ):
	#	lv=l(v)
	#	rv=r(v)
	#	G2vw[v]=[0]*2**(2*L)
	#	print v, 2**(2*L)
	#	for w in range(0, 2**(2*L) ):
	#		lw=l(w)
	#		rw=r(w)
	#		for x in range(0, 2**(2*L) ):
	#			Rrvx_G0x=R(rv,x)*G0[x]
	#			for y in range(0, 2**(2*L) ):
	#				RrvxRlvyRrwy_G0xG0y=R(lv,y)*R(rw,y)*Rrvx_G0x*G0[y]
	#				for z in range(0, 2**(2*L) ):
	#					G2vw[v][w]+=RrvxRlvyRrwy_G0xG0y*R(lw,z)*G0[z]
	
	#for u in range(0, 2**(2*L) ):
	#	G2q[u]=0
	#	lu=l(u)
	#	ru=r(u)
	#	for v in range(0, 2**(2*L) ):
	#		print u, '\t', v
	#		lv=l(v)
	#		rv=r(v)
	#		Rluv=R(lu,v)
	#		for w in range(0, 2**(2*L) ):
	#			G2q[u]+=R(l(u),v)*R(r(u),w)*G2vw[v][w]
	
	for x in range(0, 2**(2*L) ):
		for y in range(0, 2**(2*L) ):
			for z in range(0, 2**(2*L) ):
				G2[x]+=R(l(x),y, Cm)*R(r(x),z, Cm)*G1[y]*G1[z]
				G2p[x]+=R(l(x),y, Cm)*R(r(x),z, Cm)*G1p[y]*G1[z]
	
	for x in range(0, 2**(2*L) ):
		for y in range(0, 2**(2*L) ):
			for z in range(0, 2**(2*L) ):
				G3p[x]+=R(l(x),y, Cm)*R(r(x),z, Cm)*G2p[y]*G2[z]
	#			G3q[x]+=R(l(x),y)*R(r(x),z)*G2q[y]*G2[z]
				G3[x]+=R(l(x),y, Cm)*R(r(x),z, Cm)*G2[y]*G2[z]
	
	S0=0
	S1=0
	S1p=0
	S2=0
	S2p=0
	S3=0
	S3p=0
	S3q=0
	
	for x in range (0, 2**(2*L) ):
		S0+=G0[x]
		S1+=G1[x]
		S1p+=G1p[x]
		S2+=G2[x]
		S2p+=G2p[x]
		S3+=G3p[x]
		S3p+=G3p[x]
		S3q+=G3q[x]
	#	print '{0:03b}'.format(l(x) ), '{0:03b}'.format(r(x)), G3p[x], G3p[x]-G3[x]
	
	#print S0, S1, S1p, S2, S2p, S3, S3p, S3q
	
	#A is 2
	#B is 1
	#C is 0
	
	#print "A:", A, get_moment(G3p, [2] )
	#print "B:", B, get_moment(G3p, [1] )
	#print "C:", C, get_moment(G3p, [0] )
	
	#2 1 0
	#5 4 3
	
	DP1=[0,1]
	DP2=[1,2]
	DP3=[0,2]
	
	FP=[0,2]
	
	#print "D1 ", get_D(G0, G1, DP1)
	#print "D2 ", get_D(G1, G2, DP1)
	#print "D3 ", get_D(G2, G3, DP1)
	#print "f3 ", get_moment(G3, FP )
	
	#print_gen(G1p)
	#print_gen(G0)
	#print_hap(G1)
	#print_hap(G1p)
	#print_gen(G1p)
	#print_hap(G1p)
	
	#print "D1p 1 ", get_D(G0, G1p, DP1)
	
	#print "D2p 1 ", get_D(G1, G2p, DP1)
	#print "D2  1 ", get_D(G1, G2, DP1)
	
	#print "D3p 1 ", get_D(G2, G3p, DP1)
	#print "D3q 1 ", get_D(G2, G3q, DP1)
	
	#print "D3p 2", get_D(G2, G3p, DP2)
	#print "D3q 2", get_D(G2, G3q, DP2)
	
	#print "D3p 3", get_D(G2, G3p, DP3)
	#print "D3q 3", get_D(G2, G3q, DP3)
	
	#print "S3q 3", get_S(G2, G3q, DP3)
	
	print Cm[0], (get_moment(G3p, [0,1] )+get_moment(G3p, [0+L,1+L]) )/ 2., (get_moment(G3, [0,1] )+get_moment(G3, [0+L,1+L]) )/ 2.
