from sympy import exp 
from sympy import symbols 
from sympy import sqrt
from sympy import collect
from sympy import ccode

from sympy.utilities.codegen import codegen
from sympy.printing import print_ccode

#DATA
								# A	B   C	  D   E     F   G
trA, trB, trC, trD, trE, trF, trG = symbols("trMtM trMJM trHtH trHJH trMtH trMJH trE")
#parameters
a, b, c, d, e, f, l = symbols("a b c d e f l") 

#obj = A*exp(a) - 2*(exp(a)*b**2)/(exp(a) + b**2) * B + C * exp(c) - 2*(exp(c)*d**2)/(exp(c) + d**2) * D+ sqrt( exp(a) * exp(c) )*( 2 / (1+exp(e) ) -1 ) * E - 2*b*d*sqrt( exp(a)*exp(c) )*( 2/(1+exp(e) ) - 1 )/( sqrt( exp(a)*exp(c) )*( 2/(1+exp(e) ) -1)+b*d )*F+exp(f)*G
#obj = A*exp(a) - 2*(exp(a)*b**2)/(exp(a) + b**2) * B + C * exp(c) - 2*(exp(c)*d**2)/(exp(c) + d**2) * D+e*exp(e**2)* E - 2*b*d*e*exp(e**2)/(e*exp(e**2)+b*d )*F+exp(f)*G
#obj = A*exp(a) - 2*(exp(a)*b**2)/(exp(a) + b**2) * B + C * exp(c) - 2*(exp(c)*d**2)/(exp(c) + d**2) * D+e*E - 2*b*d*e/(e+b*d )*F+exp(f)*G

#(sa^2)*p(M,M)-2*( (sa^2*ma^2)/(sa^2+ma^2) )*q(M) %*% t(q(M) ) ) 
#(sd^2)*p(H,H)-2*( (sd^2*md^2)/(sd^2+md^2) )*q(H) %*% t(q(H) ) ) 
#sad/(sa*sd)*2*( 
#sad*p(H,M)-
#(sad*ma*md)/(sad+ma*md)*( q(H) %*% t(q(M))+q(M) %*% t(q(H)) ) ) 


obj = l*(1.-a*trA-c*trC-e*trE-f*trG+2*(a*b**2)/(a + b**2)*trB+2*(c*d**2)/(c + d**2)*trD+2*e*b*d/(e+b*d)*trF )

all_as=[a,b,c,d,e,f,l]

print "switch(x)"
print "{"
for x in range (0,7):
	s1=all_as[x]
	R=obj.diff(s1) 
	if(R!=0):
		print "\tcase("+str(x)+"):"
		print "\t{"
		K = ccode(R)
		print "\t\treturn", K, "; //", x
		print "\t}"
print "\tdefault:"
print "\t\treturn empty;"
print "}"
print
print "MJM=r.Mt1*r.Mt1.transpose();"
print "HJH=r.Ht1*r.Ht1.transpose();"
print "MJH=(r.Mt1*r.Ht1.transpose()+r.Ht1*r.Mt1.transpose() )/2;"
print "switch(x)"
print "{"
for x in range (0,7):
	print "\tcase("+str(x)+"):"
	print "\t\tswitch(y)"
	print "\t\t{"
	for y in range (x,7):
		s1=all_as[x]
		s2=all_as[y]

		R=obj.diff(s1,s2)

		if(R!=0):
			print "\t\t\tcase("+str(y)+"):"
			print "\t\t\t{"
			K = ccode(R)
			print "\t\t\t\treturn", K, "; //", x,y
			print "\t\t\t}"
	print "\t\t\tdefault:"
	print "\t\t\t\treturn empty;"
	print "\t\t}"
print "}"
