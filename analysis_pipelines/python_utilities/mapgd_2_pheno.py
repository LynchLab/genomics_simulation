import fileinput

for line in fileinput.input():
	if line[0]!='@':
		name, value=line.split('\t') 
		print name+' '+name+' '+value,
