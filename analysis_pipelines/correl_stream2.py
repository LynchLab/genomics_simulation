import math
import sys
from scipy import stats
import numpy
import itertools

if len(sys.argv)==2:
	print sys.argv	
	File=open(sys.argv[1])
else:
	File=sys.stdin


print "run freq rho1 rho2 rho3 P OX EX OY EY"

D=[]
Sbuffer=[]
line=File.readline()
DIST=2

XMAX=1

N=len(map(int, line.split('\t')[1:] ))/2
print N

for d1 in range(0, N*2):
	D.append([])
	print len(D)
	for d2 in range(0, N*2):
		D[d1].append([])
		for X in range (0, XMAX):
			D[d1][d2].append([])
			for Y in range (0, DIST):
				D[d1][d2][X].append([])
				for i in range (0, 2):
					D[d1][d2][X][Y].append([])
					for j in range (0, 2):
						D[d1][d2][X][Y][i].append(0)
for line in File:
	S=map(int, line.split('\t')[1:])
	s=[]
	for x in range(0, N*2, 2):
		s.append(int(S[x]==S[x+1]))
	Sbuffer.append(s)
	if( len(Sbuffer) > DIST):
		Sbuffer.pop(0)
	for x in range(0, XMAX):
		for y in range (0, len(Sbuffer) ):
			t=Sbuffer[-y-1]
			i=s[x]
			j=t[x]
			d1=sum(s)-i
			d2=sum(t)-j
#			print d1, d2, 
			D[d1][d2][x][y][i][j]+=1
			#print d, x, y, i, j, D[d][x][y][i][j]
print "done reading"
for x in range (0, XMAX):
	for y in range(1, DIST):
		for c1 in range(1, N-1):
			for c2 in range(1, N-1):
#				print x, y, c1, c2, D[c1][c2][x][y]
				d1=float(c1)/float(N-1)
				d2=float(c2)/float(N-1)
				Den=D[c1][c2][x][y][0][0]+D[c1][c2][x][y][1][0]+D[c1][c2][x][y][0][1]+D[c1][c2][x][y][1][1]
				if Den>10:
					OX=float(D[c1][c2][x][y][1][0]+D[c1][c2][x][y][1][1])/float(Den)
					OY=float(D[c1][c2][x][y][0][1]+D[c1][c2][x][y][1][1])/float(Den)
					Num=D[c1][c2][x][y][0][0]*float(0-d1)*float(0-d2)+D[c1][c2][x][y][1][0]*float(1.0-d1)*float(0-d2)+D[c1][c2][x][y][0][1]*float(0-d1)*float(1.-d2)+D[c1][c2][x][y][1][1]*float(1.-d1)*float(1.-d2)
					Num2=D[c1][c2][x][y][0][0]*float(0-OX)*float(0-OY)+D[c1][c2][x][y][1][0]*float(1.0-OX)*float(0-OY)+D[c1][c2][x][y][0][1]*float(0-OX)*float(1.-OY)+D[c1][c2][x][y][1][1]*float(1.-OX)*float(1.-OY)
					try:
						print D[c1][c2][x][y]
						print x, y, c1, c2, float(Num)/float(Den)/math.sqrt((1-d1)*d1*d2*(1-d2) ), OX, OY, float(Num2)/float(Den)/math.sqrt((1-OX)*OX*OY*(1-OY) )

					except:
						print "whoops!"					
	quit()
