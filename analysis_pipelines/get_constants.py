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
File=open(sys.argv[1])
for line in File:
	l=line.split(' ')
	if l[0]!="POS":
		p=float(l[1])
		q=1-p
		if p==0 or q==0:
			continue
		f=float(l[2])
		Aa2+=2*p*q*(1+f)
		Aad+=2*p*q*(1+f)*(2*q-2*p)
		Ad2+=2*p*q*(1+f)*(p**2-2*p*q+q**2)
		Dd2+=4*p*q*(p*q + f*(1-2*p)**2-f**2*p*q )
		ADad+=4*(p**3*q - q**3*p )*f
		ADad2+=4*(p**3*q - q**3*p )*f*(1-2*p)
		Dm-=2*p*q*f
		K+=1 / (1-p) + 1/p -3.
		N+=1.
print "var_A =", Aa2/N, "a^2 +", Aad/N, "a d +", Ad2/N, "d^2"
print "var_D =", Dd2/N, "d^2"
print "cov_A,D =", ADad/N, "a d +", ADad2/N, "d^2"
print "mu_D =", Dm/N
print "kappa =", K/N
