import sys
import math
from copy import deepcopy

def get_color(z):
	global min_z, max_z
	if max_z!=min_z:
		value=float(z-min_z)/float(max_z-min_z)
	else :
		value=1.
	G="{:02x}".format(int(value*255))
	return "\"#ff"+G+G+"\""

	return "blue"

def print_parent (ind, depth, td):
	global individuals, joins, File, linesize, pad, start, max_z, min_z
	if depth >= td:
		File.seek((int(ind.pid)-start)*linesize+pad)
		line=File.readline()
		time, id, pid, mid, sex, z=line.strip('\n').split()
		if id not in individuals.keys():
			individuals[id]=Individual(id, pid, mid, sex, float(z), depth, ind.ID2+2**(td-1) )
			joins[id]=set([ind.ID2+2**(td-1) ])
		else :
			joins[id].add( ind.ID2+2**(td-1) )
		max_z=max(max_z, float(z) )
		min_z=min(min_z, float(z) )

		File.seek((int(ind.mid)-start)*linesize+pad)
		line=File.readline()
		time, id, pid, mid, sex, z=line.strip('\n').split()
		if id not in individuals.keys():
			individuals[id]=Individual(id, pid, mid, sex, float(z), depth, ind.ID2+2**(td) )
			joins[id]=set([ind.ID2+2**(td) ])
		else :
			joins[id].add( ind.ID2+2**(td) )

		max_z=max(max_z, float(z) )
		min_z=min(min_z, float(z) )

		print_parent(individuals[ind.pid], depth, td+1)
		print_parent(individuals[ind.mid], depth, td+1)

def node_sum(H):
	rsum=0
	for x in range(0, len(H) ):
		for y in range(0, len(H[x])):
			for z in range(0, len(H[x][y]) ):
				rsum+=H[x][y][z]
	return rsum

class Individual:
	def __init__ (self, id, pid, mid, sex, z, depth, ID2):
		self.id=id
		self.pid=pid
		self.mid=mid
		self.sex=sex
		self.z=z
		self.draw=False
		self.depth=depth
		self.ID2=ID2
	def node(self):
		return self.id+":"+str(self.depth)
	def edges(self):
		return self.pid+","+self.id+":"+self.mid+","+self.id

File=open(sys.argv[1])

NameFile=open(sys.argv[2])

for line in NameFile:
	if line[0]!='@':
		names=map(int, line.strip('\n').split('\t')[1:] )
NameFile.close()

DEAPTH=int(sys.argv[3])

line  = File.readline()
pad   = len(line)

start = int(File.readline().strip('\n').split()[1])

min_z=10000
max_z=-10000

linesize=61

count={}

#is x the ancestor of y?
def ancestor (x, y):
	if (x <= y):
		return False
	m=2**(math.floor(math.log(y+1)/math.log(2) ) )
	return x%m==y%m

def rotate(H, n):
	rH=deepcopy(H)
	m=int(2**(math.floor(math.log(n+1)/math.log(2) ) ) )
	l=n+m
	r=n+2*m
	for x in range(0, len(H) ):
		for y in range(0, len(H[x])):
			for z in range(0, len(H[x][y]) ):
				if (ancestor(H[x][y][z], n) ):
					if( ancestor(H[x][y][z], l) or H[x][y][z]==l ):
						rH[x][y][z]+=m
					else:
						rH[x][y][z]-=m
	return rH

#not defined correctly for intergeneration
def swap(H, X, Y):
	#print "Attempt swap ", X, Y
	rH=deepcopy(H)
	slide=Y-X
	for x in range(0, len(H) ):
		for y in range(0, len(H[x])):
			for z in range(0, len(H[x][y]) ):
				if ( ancestor(H[x][y][z], X) or H[x][y][z]==X ):
					rH[x][y][z] += slide
				if H[x][y][z]==Y:
					rH[x][y][z] = X
				if ancestor(H[x][y][z], Y):
					print "Invalid swap"
	#print H, rH
	return rH

#H=[[], [[4,3]], [[9, 8]]]
#print swap(H,4,3)
#quit()
def min_swap(H):
	tH=deepcopy(H)
	for x in range(0, len(tH) ):
		for y in range(0, len(tH[x])):
			if len(H[x][y])!=0:
				if min(tH[x][y])!=tH[x][y][0]:
					tH=deepcopy(swap(tH,tH[x][y][0], min(tH[x][y]) ) )
	return tH

def canonize(H):
	tH=min_swap(H)
	if (node_sum(tH)!=0):
		for x in range(0, 2**len(tH) ):
			rH=rotate(tH, x)
			rH=min_swap(rH)
			#print x, rH
			#print tH
			if (node_sum(rH)<node_sum(tH) ):
				tH=deepcopy(rH)
	return tH

#N = ?

for Ind in names:
	individuals={}
	joins={}
	Ind=str(Ind)
#	print int(Ind)-start, start
	File.seek((int(Ind)-start)*linesize+pad)
	line=File.readline()
	time, id, pid, mid, sex, z=line.strip('\n').split()
	individuals[id]=Individual(id, pid, mid, sex, float(z), DEAPTH, 0)
	max_z=max(max_z, float(z) )
	min_z=min(min_z, float(z) )

	HASH=[]
	for x in range(0, DEAPTH):
		HASH.append([])

	print_parent(individuals[Ind], DEAPTH, 1)

	for key in joins.keys():
		this_join=list(joins[key])
		fmin=max([int(math.floor(math.log(x+1)/math.log(2)-1) ) for x in this_join ])
		if len(this_join)>1:
			HASH[fmin].append(this_join)
	
	THASH=""

	HASH=canonize(HASH)

	for x in range(0, DEAPTH):
		THASH=THASH+"{"
		HASH[x].sort()
		#J=[]
		#if (len(HASH[x])>1):
		#	J=sorted(HASH[x])
		#	HASH[x]=HASH[x].sort()
		for y in range(0, len(HASH[x])-1):
			THASH=THASH+'='.join(map(str,HASH[x][y]))+","
		if len(HASH[x])>0:
			rotate(HASH, 3)
			THASH=THASH+'='.join(map(str,HASH[x][-1]))
		THASH=THASH+"}"
	#if THASH != '{}'*DEAPTH:
	print names.index(int(Ind)), THASH

	if THASH not in count.keys():
		count[THASH]=1
	else :
		count[THASH]+=1
#	HASH=?

for key in count.keys():
	print key, count[key]
