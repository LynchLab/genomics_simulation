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
		#print "seeking ", ind.pid, " and ", ind.mid, " parents of ", ind.id
		if ind.pid=="NaN":
			print "Error, out of depth"
		else:
			if int(ind.pid)>=start:
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
			else :
				id=ind.pid
				z=0
				if id not in individuals.keys():
					individuals[id]=Individual(id, "NaN", "NaN", "NaN", float(z), depth, ind.ID2+2**(td-1) )
					joins[id]=set([ind.ID2+2**(td-1) ])
				else :
					joins[id].add( ind.ID2+2**(td-1) )
		#print "seeking ", ind.mid, start
		if ind.mid=="NaN":
			print "Error, out of depth"
		else:
			if int(ind.mid)>=start:
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
			else:
				id=ind.mid
				z=0
				if id not in individuals.keys():
					individuals[id]=Individual(id, "NaN", "NaN", "NaN", float(z), depth, ind.ID2+2**(td-1) )
					joins[id]=set([ind.ID2+2**(td) ])
				else :
					joins[id].add( ind.ID2+2**(td) )

		if ind.pid!="NaN":
			print_parent(individuals[ind.pid], depth, td+1)
		if ind.mid!="NaN":
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

Ind1=(sys.argv[2])

DEAPTH=int(sys.argv[3])

N = int(sys.argv[4])

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
	rH=deepcopy(H)
	slide=Y-X
	for x in range(0, len(H) ):
		for y in range(0, len(H[x])):
			for z in range(0, len(H[x][y]) ):
				if (ancestor(H[x][y][z], X) or H[x][y][z]==X ):
					rH[x][y][z] += slide
				if ( H[x][y][z]==Y ):
					rH[x][y][z] = X
				if (ancestor(H[x][y][z], Y) ):
					print "Invalid swap"
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

#H=[[],[],[[3,12,14]],[[17,24,26]]]
#print canonize(H)
#quit()


for t in range(N-1, -1, -1):
	individuals={}
	joins={}
	Ind=str(int(Ind1)-t)
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
		#print this_join
		if len(this_join)>1:
			HASH[fmin].append(this_join)
	
	THASH=""

	#print "uncanonized", HASH
	HASH=canonize(HASH)

	for x in range(0, DEAPTH):
		THASH=THASH+"{"
		for y in range(0, len(HASH[x])-1):
			THASH=THASH+'='.join(map(str,HASH[x][y]))+","
		if len(HASH[x])>0:
			rotate(HASH, 3)
			THASH=THASH+'='.join(map(str,HASH[x][-1]))
		THASH=THASH+"}"
	#if THASH != '{}'*DEAPTH:
	print str(int(Ind) ), THASH

	if THASH not in count.keys():
		count[THASH]=1
	else :
		count[THASH]+=1
#	HASH=?

for key in count.keys():
	print key, count[key]
