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

D=[]
line=File.readline()
N=len(map(int, line.split('\t')[1:] ))/2

XMAX=1
for d in range(0, N*2):
	D.append([])
	for X in range (0, XMAX):
		D[d].append([])
		for Y in range (0, N):
			D[d][X].append([])
			for i in range (0, 3):
				D[d][X][Y].append([])
				for j in range (0, 3):
					D[d][X][Y][i].append(0)
for line in File:
	S=map(int, line.split('\t')[1:] )
#	for x in range(0, N ):
	x=0
	for y in range (x+1, N ):
		i=S[x*2]+S[x*2+1]
		j=S[y*2]+S[y*2+1]
		d=sum(S)-i-j
		D[d][x][y][i][j]+=1
		#print d, x, y, i, j, D[d][x][y][i][j]

print "run freq rho1 rho2 rho3 P OX EX OY EY"
for x in range (0, XMAX):
	for y in range(x+1, N):
		for c in range(0, N*2-4):
			d=float(c)/float(N*2-4)
			Den=D[c][x][y][0][0]+D[c][x][y][1][0]+D[c][x][y][2][0]+D[c][x][y][0][1]+D[c][x][y][1][1]+D[c][x][y][2][1]+D[c][x][y][0][2]+D[c][x][y][1][2]+D[c][x][y][2][2]
			if (Den>10 and d!=0 and d!=1):
				OX=float(D[c][x][y][1][0]+D[c][x][y][1][1]+D[c][x][y][1][2]+2*(D[c][x][y][2][0]+D[c][x][y][2][1]+D[c][x][y][2][2]) )/float(Den*2)
				OY=float(D[c][x][y][0][1]+D[c][x][y][1][1]+D[c][x][y][2][1]+2*(D[c][x][y][0][2]+D[c][x][y][1][2]+D[c][x][y][2][2]) )/float(Den*2)
				EX=d
				EY=d
				Num1=D[c][x][y][0][0]*float(0-d)*float(0-d)+D[c][x][y][1][0]*float(0.5-d)*float(0-d)+D[c][x][y][2][0]*float(1.-d)*float(0-d)+D[c][x][y][0][1]*float(0-d)*float(0.5-d)+D[c][x][y][1][1]*float(0.5-d)*float(0.5-d)+D[c][x][y][2][1]*float(1.-d)*float(0.5-d)+D[c][x][y][0][2]*float(0-d)*float(1.-d)+D[c][x][y][1][2]*float(0.5-d)*float(1.-d)+D[c][x][y][2][2]*float(1.-d)*float(1.-d)
				Num2=D[c][x][y][0][0]*float(0-EX)*float(0-EY)+D[c][x][y][1][0]*float(0.5-EX)*float(0-EY)+D[c][x][y][2][0]*float(1.-EX)*float(0-EY)+D[c][x][y][0][1]*float(0-EX)*float(0.5-EY)+D[c][x][y][1][1]*float(0.5-EX)*float(0.5-EY)+D[c][x][y][2][1]*float(1.-EX)*float(0.5-EY)+D[c][x][y][0][2]*float(0-EX)*float(1.-EY)+D[c][x][y][1][2]*float(0.5-EX)*float(1.-EY)+D[c][x][y][2][2]*float(1.-EX)*float(1.-EY)
				Num3=D[c][x][y][0][0]*float(0-OX)*float(0-OY)+D[c][x][y][1][0]*float(0.5-OX)*float(0-OY)+D[c][x][y][2][0]*float(1.-OX)*float(0-OY)+D[c][x][y][0][1]*float(0-OX)*float(0.5-OY)+D[c][x][y][1][1]*float(0.5-OX)*float(0.5-OY)+D[c][x][y][2][1]*float(1.-OX)*float(0.5-OY)+D[c][x][y][0][2]*float(0-OX)*float(1.-OY)+D[c][x][y][1][2]*float(0.5-OX)*float(1.-OY)+D[c][x][y][2][2]*float(1.-OX)*float(1.-OY)
				try:
					print x*N+(y-x), c, float(Num1)/float(Den)/(1-d)/d, float(Num2)/float(Den)/(1-d)/d, float(Num3)/float(Den)/math.sqrt(OX*(1-OX)*OY*(1-OY) ), d, OX, EX, OY, EY
				except:
					print "whoops!"
	quit()
