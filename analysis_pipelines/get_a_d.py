import sys
import scipy.optimize
import numpy
import operator

def fun(o):
	a=o[0]
	d=o[1]
	global X
	global s
	#return (s[0]-X[0]*a**2-X[1]*a*d-X[2]*d**2)**2+(s[1]-X[3]*d**2)**2+(s[2]-X[4]*a*d-X[5]*d**2)**2+(s[3]-X[6]*d)**2
	return (s[0]-X[0]*a**2-X[1]*a*d-X[2]*d**2)**2+(s[2]-X[4]*a*d-X[5]*d**2)**2

X=map(float, sys.argv[1:8])
s=map(float, sys.argv[8:12])

print X
print s

results={}

for x in range(0, 100):
	c=[numpy.random.uniform(-50,50), numpy.random.uniform(-50, 50) ]
	f=scipy.optimize.minimize(fun, c, method="SLSQP")
	try:
		results[str(round(f.x[0],3))+","+str(round(f.x[1],3))][1]+=1
	except:
		results[str(round(f.x[0],3))+","+str(round(f.x[1],3))]=[ f.fun, 1]

sorted_results = sorted(results.items(), key=operator.itemgetter(1))
for key in sorted_results:
	print key[0], key[1][0], key[1][1]
