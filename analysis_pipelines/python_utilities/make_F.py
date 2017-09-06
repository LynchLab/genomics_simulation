import sys
import numpy

N=int(sys.argv[1])

F=numpy.zeros( (N, N) )
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
	f[x]+=l[2]
	f[y]+=l[4]

for x in range(0, N):
	f[x]=f[x]/N

fbar=numpy.mean(f)

for x in range (0, N):
	for y in range(0, N):
		F[x,y] = fbar*(fbar-f[x]-f[y])	
		F[y,x] = fbar*(fbar-f[x]-f[y])
		F[x,x] = fbar**2
		F[y,y] = fbar**2
#print "Fone"
#quit()


for x in range (0, N):
	for y in range (0, N-1):
		print str(F[x][y])+", ",
	print str(F[x][N-1])
