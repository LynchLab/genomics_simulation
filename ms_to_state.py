import math
import sys
from bitarray import bitarray
import numpy
import fileinput
WORD=64

class locus:
	def __init__ (self, pos, N):
		global WORD
		self.pos=pos
		self.state=bitarray('0'*N)
		self.N=N

def readms (FILENAME):
	loci=[]
	FILE=""
	if FILENAME!='-':
		FILE=open(FILENAME)
	else:
		FILE=sys.stdin
	line=FILE.readline().strip('\n')
	N=int(line.strip('\n').split()[1] )
	#print N
	while(line !="//"):
		line=FILE.readline().strip('\n')
	#	print line
	segsites=int(FILE.readline().split(':')[1])
	#print segsites
	line=FILE.readline().strip('\n').split(':')[1].split()
	#print line
	for x in line:
		loci.append(locus(float(x), N))
	for x in range(0, N):
		line=FILE.readline().strip('\n')
		for y in range(0, segsites) :
			loci[y].state[x]=(line[y]=='1')
	return loci

def print_state (loci):
	GS=int(len(loci)/32)
	#del loci[GS*32:]

	NS=loci[0].N/2
	print "@NAME:STATE	VERSION:0.0	FORMAT:TEXT	CONCATENATED"
	print "@NS:"+str(NS)+"\tGS:"+str(GS)+"\tBS:"+str(2)
	for x in range(0, GS*32):
		print str(loci[x].state.to01()[0])+str(loci[x].state.to01()[1]),
		for y in range(2, NS*2, 2):
			print '\t'+str(loci[x].state.to01()[y] )+str(loci[x].state.to01()[y+1] ),
		print 
	print "@END_TABLE"

loci=readms(sys.argv[1])
print_state(loci)
