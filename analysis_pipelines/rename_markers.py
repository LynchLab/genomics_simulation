import fileinput
x=0
for line in fileinput.input():
	line=line.split('\t')
	print line[0]+"\tSNP"+str(x)+'\t'+'\t'.join(line[2:]),
	x+=1
