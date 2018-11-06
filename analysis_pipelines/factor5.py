from sympy import exp 
from sympy import symbols 
from sympy import sqrt
from sympy import collect
from sympy import ccode

from sympy.utilities.codegen import codegen
from sympy.printing import print_ccode


A, B, C, D, E, F, G, trA, trB, trC, trD, trE, trF, trG = symbols("r.MtM MJM r.HtH HJH r.MtH MJH v trMtM trMJM trHtH trHJH trMtH trMJH trF")

a, b, c, d, e = symbols("a b c d e") 

#expr = A*exp(a) - 2*(exp(a)*b**2)/(exp(a) + b**2) * B + C * exp(c) - 2*(exp(c)*d**2)/(exp(c) + d**2) * D+ sqrt( exp(a) * exp(c) )*( 2 / (1+exp(e) ) -1 ) * E - 2*b*d*sqrt( exp(a)*exp(c) )*( 2/(1+exp(e) ) - 1 )/( sqrt( exp(a)*exp(c) )*( 2/(1+exp(e) ) -1)+b*d )*F+exp(f)*G
#expr = A*exp(a) - 2*(exp(a)*b**2)/(exp(a) + b**2) * B + C * exp(c) - 2*(exp(c)*d**2)/(exp(c) + d**2) * D+e*exp(e**2)* E - 2*b*d*e*exp(e**2)/(e*exp(e**2)+b*d )*F+exp(f)*G
#expr = A*exp(a) - 2*(exp(a)*b**2)/(exp(a) + b**2) * B + C * exp(c) - 2*(exp(c)*d**2)/(exp(c) + d**2) * D+e*E - 2*b*d*e/(e+b*d )*F+exp(f)*G

expr = A*a -2* (a*b**2)/(a + b**2)*B +C*c - 2 * (c*d**2)/(c + d**2) * D+(e * E - 2*e*b*d/(e+b*d)*F )+(1.-a*trA-c*trC-e*trE+2*(a*b**2)/(a + b**2)*trB+2*(c*d**2)/(c + d**2)*trD+2*e*b*d/(e+b*d)*trF )*G

all_as=[a,b,c,d,e]

print "MJM=r.Mt1*r.Mt1.transpose();"
print "HJH=r.Ht1*r.Ht1.transpose();"
print "MJH=(r.Mt1*r.Ht1.transpose()+r.Ht1*r.Mt1.transpose() )/2;"
print "switch(x)"
print "{"
for x in range (0,5):
	s1=all_as[x]
	R=collect(expr.diff(s1), [A, B, C, D, E, F])
	if(R!=0):
		print "\tcase("+str(x)+"):"
		K = ccode(R)
		print "\t\treturn", K, "; //", x
print "\tdefault:"
print "\t\treturn empty;"
print "}"
print
print "MJM=r.Mt1*r.Mt1.transpose();"
print "HJH=r.Ht1*r.Ht1.transpose();"
print "MJH=(r.Mt1*r.Ht1.transpose()+r.Ht1*r.Mt1.transpose() )/2;"
print "switch(x)"
print "{"
for x in range (0,5):
	print "\tcase("+str(x)+"):"
	print "\t\tswitch(y)"
	print "\t\t{"
	for y in range (x,5):
		s1=all_as[x]
		s2=all_as[y]

		R=collect(expr.diff(s1,s2), [A, B, C, D, E, F])

		if(R!=0):
			print "\t\t\tcase("+str(y)+"):"
			K = ccode(R)
			print "\t\t\t\treturn", K, "; //", x,y
	print "\t\t\tdefault:"
	print "\t\t\t\treturn empty;"
	print "\t\t}"
print "}"
