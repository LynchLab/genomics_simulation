import sys
File=open(sys.argv[1])
out={}
for line in File:
	line=line.strip('\n').split(' ')
	while(True):
		try:
			get=out[line[1]]
			line[1]=str(int(line[1])+1)
		except:
			out[line[1]]=line
			break
keys=map(int, out.keys() )
keys.sort()
for k in keys:
	print ' '.join(out[str(k)])

