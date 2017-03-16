from scipy import optimize
from scipy import linalg
from scipy import stats
import numpy
import math

#def sqrtm(M):
#    (U,S,VT) = LinearAlgebra.singular_value_decomposition(M)
#    D = MLab.diag(sqrt(S))
#    return matrixmultiply(matrixmultiply(U,D),VT

def vmodel (x):
	global N
	global irA
	global irD
	a=numpy.matrix(x).reshape(N,1)
	d=z-a
	lnsum=0
	u1=numpy.array(irA * a)
	u2=numpy.array(irD * d)
	for x in range(0, N):
		print u1[x], u2[x]
		p=stats.norm(0, 1).pdf(u1[x])
		if p!=0:
			lnsum+=math.log(p)
		else:
			lnsum-=350
		p=stats.norm(0, 1).pdf(u2[x])
		if p!=0:
			lnsum+=math.log(p)
		else:
			lnsum-=350
	return (-lnsum)
def model (x):
	global N
	global irA
	global irD
	lnsum=0
	a=numpy.matrix(x).reshape(N,1)
	d=z-a
	u1=numpy.array(irA * a)
	u2=numpy.array(irD * d)
	for x in range(0, N):
		p=stats.norm(0, 1).pdf(u1[x])
		if p!=0:
			lnsum+=math.log(p)
		else:
			lnsum-=9223372
		p=stats.norm(0, 1).pdf(u2[x])
		if p!=0:
			lnsum+=math.log(p)
		else:
			lnsum-=9223372
	return (-lnsum)
#+sum( numpy.array(a+d-z)**2)/10000

def func(x):
	#print x
	return model(x)

def norm(x):
	return math.sqrt(numpy.mean(numpy.array(x)**2))

zfile=open("zfile.txt")
Afile=open("Afile.txt")
Dfile=open("Dfile.txt")

z=[]
for line in zfile:
	z.append(float(line) )	
N=len(z)
A=numpy.zeros(shape=(N,N))
D=numpy.zeros(shape=(N,N))
x=0
for line in Afile:
	A[x,] = map(float, line.split('\t'))
	x+=1
x=0
for line in Dfile:
	D[x,] = map(float, line.split('\t'))
	x+=1

rA=linalg.sqrtm(A)
rD=linalg.sqrtm(D)
irA=linalg.inv(rA)
irD=linalg.inv(rD)

sigma_a=numpy.sqrt(numpy.mean(numpy.diagonal(A)))
sigma_d=numpy.sqrt(numpy.mean(numpy.diagonal(D)))
rescale=numpy.sqrt(norm(z)**2/ (sigma_a**2+sigma_d**2) )
sigma_a=sigma_a*rescale
sigma_d=sigma_d*rescale
sigma_z=norm(z)

z=numpy.matrix(z).reshape(N, 1)
#Q using a solution for a as the starting point
iA=linalg.inv(A)
Q=linalg.inv(sigma_d**2/sigma_a**2*D*iA)*z

z_perp=numpy.matrix(Q).reshape(N,1)
z_perp=z_perp-numpy.mean(z_perp)
z_perp=z_perp- numpy.asscalar((z_perp.transpose() * z)/(z.transpose() * z))*z
z_perp=z_perp/norm(z_perp)*sigma_z

a_theta=sigma_a**2/(sigma_a**2+sigma_d**2)*z+(sigma_a*sigma_d)/(sigma_a**2+sigma_d**2)*z_perp
d_theta=sigma_d**2/(sigma_a**2+sigma_d**2)*z-(sigma_a*sigma_d)/(sigma_a**2+sigma_d**2)*z_perp

cons = ({'type': 'eq','fun' : lambda x: sum(x) } 
	,{'type': 'eq','fun' : lambda x: norm(x)-sigma_a }
	,{'type': 'eq','fun' : lambda x: ((z-numpy.matrix(x).reshape(N,1) ).transpose()*numpy.matrix(x).reshape(N,1) )[0,0] }
	)
B=10
my_bounds=[[-B, B]]*N
res = optimize.minimize(func, a_theta,  options={'disp':True, 'maxiter':1000}, bounds=my_bounds, constraints=cons, method='SLSQP')
#res = optimize.minimize(func, init, constraints=cons, method='SLSQP', options={'disp': True})
a=numpy.matrix(res['x']).reshape(N,1)
d=z-a

f = open('solfile.txt', 'w')

f.write("a\td\tu1\tu2\n")

u1=numpy.array(irA * a)
u2=numpy.array(irD * d)


for x in range(0, N):
	f.write(str(a[x][0,0])+'\t'+str(d[x][0,0])+'\t'+str(u1[x][0])+'\t'+str(u2[x][0])+"\n")
print norm(a), sigma_a
print norm(d), sigma_d

