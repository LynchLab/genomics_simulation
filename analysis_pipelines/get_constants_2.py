import sys
N=0
Aa2=0
Ad2=0
Aad=0
Dd2=0
ADad=0
ADad2=0
Dm=0
K=0
a=float(sys.argv[2])
d=float(sys.argv[3])
File=open(sys.argv[1])
for line in File:
	l=line.split(' ')
	if l[0]!="POS":
		p=float(l[1])
		q=1-p
		if p==0 or q==0:
			continue
		f=float(l[2])
		#f=10
		Aa2+=2*p*q*(1+f)
		Aad+=2*p*q*(1+f)*(2*q-2*p)
		Ad2+=2*p*q*(1+f)*(p**2-2*p*q+q**2)
		Dd2+=4*p*q*(p*q + f*(1-2*p)**2-f**2*p*q )
		ADad+=4*(p**3*q - q**3*p )*f
		ADad2+=4*(p**3*q - q**3*p )*f*(1-2*p)
		Dm-=2*f*p*q
		K+=1 / (1-p) + 1/p -3.
		N+=1

T=Aa2*a**2/N+Aad*a*d/N+Ad2*d**2/N+Dd2*d**2/N+Dd2*d**2/N
print Aa2*a**2/N+Aad*a*d/N+Ad2*d**2/N, (Aa2*a**2/N+Aad*a*d/N+Ad2*d**2/N) / T
print Dd2*d**2/N, Dd2*d**2/N / T*2
print -(ADad*a*d)/N-(ADad2*d**2)/N, -((ADad*a*d/N)+(ADad2*d**2/N)) / T
print d*Dm/N
print K/N
