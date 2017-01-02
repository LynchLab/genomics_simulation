#!/usr/bin/python
import sys
var_file=open(sys.argv[1])
print "POS", "VR_FREQ"
ACI=9
ANI=10
line="#"

for x in range(0, 6):
	line=var_file.readline()

for line in var_file:
	line=line.strip('\n').split('\t')
	POS=line[1]
	Ref=line[2]
	Major=line[3]
	Minor=line[4]

	if len(line[4])>1:
		continue
	keys={}
	freq=float(line[7].split('/')[0])
	if Ref==Major:
		VR_FREQ=1-freq
	elif Ref==Minor:
		VR_FREQ=freq
	else:
		VR_FREQ="NaN"
	print POS, VR_FREQ
