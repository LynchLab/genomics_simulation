import sys
import numpy

N=int(sys.argv[1])

D=numpy.zeros( (N, N) )
f=numpy.zeros( (N) )
names={}

for line in sys.stdin:
	l=line.split('\t')
	if l[0]=="@SAMPLE_X":
		break
#print "HI"

for line in sys.stdin:
	l=line.strip('\n').split('\t')
	if len(l)<7:
	#	print "broken!\n"
		break
	#print l
	#exit
	try:
		x=names[l[0]]
	except:
		x=len(names.keys())
		names[l[0]]=x
	try:
		y=names[l[1]]
	except:
		y=len(names.keys() )
		names[l[1]]=y
	l=map(float, l[:-1])
	D[x,y] = l[14]	
	D[y,x] = l[14]	
	D[x,x] += 1.-l[2]
	D[y,y] += 1.-l[4]
#print "Done"
#quit()

for x in range(0, N):
	D[x,x]=D[x,x]/N
#	f[x]=D[x,x]

#f_bar=numpy.mean(f)

for x in range (0, N):
	for y in range (0, N-1):
		print str(D[x][y] )+", ",
	print str(D[x][N-1] )
