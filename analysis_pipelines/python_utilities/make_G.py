import sys
import numpy

N=int(sys.argv[1])

G=numpy.zeros( (N, N) )
f=numpy.zeros( (N) )

for line in sys.stdin:
	l=line.split('\t')
	if l[0]=="@SAMPLE_X":
		break
#print "HI"
names={}

for line in sys.stdin:
	l=line.strip('\n').split('\t')
	if len(l)<7:
	#	print "broken!\n"
		break
	l=map(float, l[:-1])
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
	#print l
	#exit
	G[x,y] = l[8]+l[10]	
	G[y,x] = l[8]+l[10]
	G[x,x] += l[2]
	G[y,y] += l[4]
#print "Gone"
#quit()

for x in range(0, N):
	G[x,x]=G[x,x]/N
#	f[x]=G[x,x]
 
#f_bar=numpy.mean(f)

for x in range (0, N):
	for y in range (0, N-1):
		print str(G[x][y])+", ",
	print str(G[x][N-1])
