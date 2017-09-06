#!/usr/bin/python
import sys
var_file=open(sys.argv[1])
print "POS", "VR_FREQ", "F_STAT", "QUAL"
ACI=9
ANI=10
line="#"

while line[0]=="#":
	line=var_file.readline()

for line in var_file:
	line=line.strip('\n').split('\t')
	POS=line[1]
	if len(line[4])>1:
		continue
	keys={}
	for l in line[7].split(';'):
		try:
			key, value=l.split('=')
			keys[key]=value
		except:
			K=0
	try:
		if "AF" not in keys.keys():
			AC=float(keys["AC"])
			AN=float(keys["AN"])
			VR_FREQ=AC/AN
			Q=line[5]
		else:
			VR_FREQ=float(keys["AF"])
			Q=line[5]
	except:
		print line[7]
		continue
	try:
		H=float(keys["InbreedingCoeff"])
	except:
		H=0
	print POS, VR_FREQ, H, Q
