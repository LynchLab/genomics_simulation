import sys
freq=[]
var=[]
state_file=open(sys.argv[1])
var_file=open(sys.argv[2])
for line in state_file:
	line=line.strip('\n').split('\t')
	N=float(len(line)-1)
	freq.append(sum(map(float, line[1:]) )/N )
for line in var_file:
	line=line.strip('\n').split(' ')
	var.append(line[1])
print "POS", "VR_FREQ"
for x in range(0, len(var) ):
	print var[x], freq[x]
