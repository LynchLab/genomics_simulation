
import sys

MAP=open(sys.argv[1])
loc={}
MAP.readline()
for line in MAP:
	line=line.strip('\n').split('\t')
	try:
		name=line[0]
		bp=int(line[1])
		cm=float(line[2])
	except:
		continue
	try:
		loc[name][bp]=cm+int(name)*1000
	except:
		loc[name]={}
		loc[name][bp]=cm+int(name)*1000
MAP.close()

VCF=open(sys.argv[2])
i=0
last_name=""
for line in VCF:
	if line[0]=="#":
		continue
	line=line.strip('\n').split('\t')
	try:
		name=line[0]
		bp=int(line[1])
	except:
		continue
	if (name!=last_name):
		last_name=name
		i=0
		keyList=sorted(loc[name].keys() )
		next_cm=int(name)*1000
		next_bp=0
	if (bp>next_bp):
		if (i!=len(keyList)-1 ):
			last_cm=next_cm
			last_bp=next_bp
			i+=1
			next_cm=loc[name][keyList[i]]
			next_bp=keyList[i]
	if (next_bp!=last_bp):
		print name, bp, last_cm+(next_cm-last_cm)*(bp-last_bp)/(next_bp-last_bp)
	else:
		print name, bp, "NaN"
VCF.close()
