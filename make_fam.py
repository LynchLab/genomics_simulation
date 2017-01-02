import sys
File=open(sys.argv[1])
x=0
for line in File:
        line=line.split()
	x+=1
        l=[str(x)]+line[1:]
        print '\t'.join(l)
