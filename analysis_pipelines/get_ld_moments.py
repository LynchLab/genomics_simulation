import sys
import math
var=[]
state_file=open(sys.argv[1])
var_file=open(sys.argv[2])
line=[]

DIST=int(sys.argv[3])

def getld(line1, line2):

	N=float(len(line1))
	count=[0]*16
	#genotype encoding 0=0|0 1=1|0 2=0|1 3=1|1

	for x in range(0, len(line1) ):
		count[line1[x]+line2[x]*4]+=1
	A=count[1]+count[2]+count[3]*2+count[5]+count[6]+count[7]*2+count[9]+count[10]+count[11]*2+count[13]+count[14]+count[15]*2
	B=count[4]+count[5]+count[6]+count[7]+count[8]+count[9]+count[10]+count[11]+(count[12]+count[13]+count[14]+count[15])*2
	AA=count[3]+count[7]+count[11]+count[15]
	BB=count[12]+count[13]+count[14]+count[15]
	ABc=count[5]+count[7]+count[10]+count[11]+count[13]+count[14]+count[15]*2
	ABt=count[6]+count[7]+count[9]+count[11]+count[13]+count[14]+count[15]*2
	AAB=count[7]+count[11]+count[15]*2
	ABB=count[13]+count[14]+count[15]*2
	AABB=count[15]*2
	SS=float(sum(count)*2)
	print count
	if SS==0:
		return "NaN"
	A/=SS
	B/=SS
	AA/=SS
	BB/=SS
	ABc/=SS
	ABt/=SS
	AAB/=SS
	ABB/=SS
	AABB/=SS
	return A, B, AA, BB, ABc, ABt, AAB, ABB, AABB

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

print "POS1", "POS2", "A", "B", "AA", "BB", "ABc", "ABt", "AAB", "ABB", "AABB"
for x in range(0, len(var) ):
	posx=var[x]
	for y in range(x+1, min(len(var), x+DIST) ):
		posy=var[y]
		if (posy-posx)<DIST:
			J=getld(line[x], line[y])
			print var[x], var[y], J, (J[4]-J[0]*J[1])/math.sqrt(J[0]*(1-J[0])*J[1]*(1-J[1]) )
