import sys
import math
var=[]
state_file=open(sys.argv[1])
var_file=open(sys.argv[2])
line=[]

DIST=int(sys.argv[3])

def getld(line1, line2):

	N=float(len(line1))
	p_i=0
	p_j=0
	f_i=0
	f_j=0
	mu=0
	s22=0
	count=[0,1,1,2]
	#genotype encoding 0=0|0 1=1|0 2=0|1 3=1|1

	for x in range(0, len(line1) ):
		s22+=count[line1[x] & line2[x] ]

		p_i+=count[line1[x] ]
		p_j+=count[line2[x] ]
	p_i=p_i/float(N)/2.
	p_j=p_j/float(N)/2.
	q_i=(1-p_i)
	q_j=(1-p_j)
	s22=s22/float(N)/2.
	count_i=[p_i**2,-p_i*q_i,-p_i*q_i,q_i**2]
	count_j=[p_j**2,-p_j*q_j,-p_j*q_j,q_j**2]
	for x in range(0, len(line1) ):
                mu+=(count_i[line1[x]]*count_j[line2[x]])
		f_i+=count_i[line1[x] ]
		f_j+=count_j[line2[x] ]

	f_i=f_i/float(N)
	f_j=f_j/float(N)
	mu=mu/float(N)

	k1i=(p_i*q_i)**2
	k2i=p_i*q_i*(3*p_i**2-3*p_i+1)
	k1j=(p_j*q_j)**2
	k2j=p_j*q_j*(3*p_j*p_j-3*p_j+1)
	k1=math.sqrt(k1i*k1j);
	k2=math.sqrt(k2i*k2j);

	s2_i=p_i*q_i
	s2_j=p_j*q_j
#	den=math.sqrt(s2_i*s2_j)
	if p_i>0 and p_j>0:
		return [2*(s22-p_i*p_j), 4*mu, 4*(f_j*s2_j*(p_i-p_j)**2+f_i*s2_i*(p_i-p_j)**2+f_i*f_j*s2_i*s2_j), k1, k2, math.sqrt(p_i*p_j)]
	else:
		return ["NaN", "NaN", "NaN"]

for thisline in state_file:
	Ns=map(int, thisline.strip('\n').split('\t')[1:]) 
	sums=[]
	#genotype encoding 0=0|0 1=1|0 2=0|1 3=1|1
	for x in range(0, len(Ns), 2):
		sums.append(Ns[x]+Ns[x+1]*2)
	line.append(sums[:])
#	print line
		
for thisline in var_file:
	thisline=thisline.strip('\n').split(' ')
	var.append(int(thisline[1]))

print "X Y POS1 POS2 RSQ MU F K1 K2 P"
for x in range(0, len(var) ):
	posx=var[x]
	for y in range(x, min(len(var), x+DIST) ):
		posy=var[y]
		if (posy-posx)<DIST:
			print x, y, var[x], var[y], ' '.join(map(str, getld(line[x], line[y])))
