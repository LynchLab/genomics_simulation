import numpy
import sys
import itertools

M=[]
File=open(sys.argv[1])
File.readline()
File.readline()

a=0.01
d=-0.01
N=12

for line in File:
	line=line.split()
	m=map(int, line[:])
	if len(m)/2==0:
		break
	N=len(m)/2
	p=sum(m)/float(2*N)
#	print m
#	quit()
	#print p
	if (p>0 and p<1):
		q=1.-p
		dm=[]
		for x in range(0, N):
			dm.append(m[x*2]+m[x*2+1])
		#print dm
		H = sum(tm == 1 for tm in dm)
		D=2.*p*q-float(H)/float(N);
		alpha=a+d*(q-p)
		#this_beta=[alpha*(0.-2*p), alpha*(1.-2.*p), alpha*(2.-2.*p) ]
		this_beta=[-p*2.*alpha,(q-p)*alpha,2.*q*alpha]
		this_delta=[d*(-2*p**2+D),d*(2*p*q+D), d*(-2*q**2+D)]
		t_b=[]
		t_d=[]

		for x in range(0, N):
			t_b.append(this_beta[dm[x]])
			t_d.append(this_delta[dm[x]])
#		print t_b+t_d
#		quit()

		M.append(t_b+t_d)

M=numpy.matrix(M)
v=numpy.cov(M.transpose() )
#print N
A=v[0:N, 0:N]
D=v[N:,N:]
G=v[0:N,N:]+v[N:,0:N]

#A=A/numpy.mean(numpy.diag(A))
#D=D/numpy.mean(numpy.diag(D))
#G=G/numpy.mean(numpy.diag(G))

#A=A.tolist()

print "@NAME:CORREL	VERSION:0.0	FORMAT:TEXT	CONCATENATED"
print "@SAMPLE_X	SAMPLE_Y	BETA_XY	DELTA_XY	GAMMA_XY"
for x in range(0, N):
	for y in range(x, N):
		print str(x)+"\t"+str(y)+"\t"+str(A[x,y])+"\t"+str(D[x,y])+"\t"+str(G[x, y])
print "junk"
