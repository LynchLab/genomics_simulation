import sys
import numpy

N=int(sys.argv[1])

d=numpy.zeros( (N, N) )
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

	d[x,y] = l[12] 
	d[y,x] = l[12] 
	d[x,x] += l[2] 
	d[y,y] += l[4] 
#print "done"
#quit()

for x in range(0, N):
	d[x,x]=d[x,x]/N
#	f[x]=d[x,x]

#f_bar=numpy.mean(f)

for x in range (0, N):
	for y in range (0, N-1):
		print str(d[x][y] )+", ",
	print str(d[x][N-1]  )
