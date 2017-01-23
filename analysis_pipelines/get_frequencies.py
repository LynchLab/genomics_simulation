import sys
freq=[]
var=[]
fstat=[]
state_file=open(sys.argv[1])
var_file=open(sys.argv[2])
for line in state_file:
	line=line.strip('\n').split('\t')
	N=float(len(line)-1)
	s=map(float, line[1:])
	p=(sum(s)/N )
	H=0
	for x in range(0, int(N), 2):
		H+=int(s[x]!=s[x+1])
#		print s[x], s[x+1], "|", 
#	print H
	H/=(N/2.)
	q=1-p
	freq.append(p)
	if p>0 and p<1:
		fstat.append(1-H/(2.*p*q) )
	else:
		fstat.append( 0 )
		

for line in var_file:
	line=line.strip('\n').split(' ')
	var.append(line[1])

print "POS", "VR_FREQ", "F_STAT"
for x in range(0, len(var) ):
	print var[x], freq[x], fstat[x]
