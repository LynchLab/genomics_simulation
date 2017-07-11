import sys
var=[]
state_file=open(sys.argv[1])
var_file=open(sys.argv[2])
line=[]

DIST=int(sys.argv[3])

def getld(line1, line2):

	N=float(len(line1))
	p1=0
	p2=0
	s22=0
	count=[0,1,1,2]
	#genotype encoding 0=0|0 1=1|0 2=0|1 3=1|1

	for x in range(0, len(line1) ):
		s22+=count[line1[x] & line2[x] ]
		p1+=count[line1[x] ]
		p2+=count[line2[x] ]
	p1=p1/float(N)/2.
	p2=p2/float(N)/2.
	s22=s22/float(N)/2.
	q1=(1-p1)
	q2=(1-p2)
	if p1>0 and p2>0:
		return (s22-p1*p2)**2/(p1*p2*q1*q2)
	else:
		return "NaN"

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

print "POS1", "POS2", "RSQ"
for x in range(0, len(var) ):
	posx=var[x]
	for y in range(x+1, min(len(var), x+DIST) ):
		posy=var[y]
		if (posy-posx)<DIST:
			print var[x], var[y], getld(line[x], line[y])
