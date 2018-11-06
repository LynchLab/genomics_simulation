import sys
import numpy
import bisect

c=0.003125
u=0.00000

recombination_rate=[0,c] 	#male/female
ploidy=[1,2] 			#male/female
mutation_rate=[u,u] 		#male/female

class chromosome:
	def __init__ (self, size):
		self.size=size
		self.poly=[]
	def insert (self, poly):
		for p in poly:
			bisect.insort(self.poly, p)
	def splice (self, chrm, loc):
		self.poly=self.poly[:bisect.bisect_right(self.poly, loc)]+chrm.poly[bisect.bisect_left(chrm.poly, loc):]

def read_ms (File):
	line=""
	line=File.readline().strip('\n')
	size=int(line.strip('\n').split('-')[-1].split(' ')[2])
	while (line!="//"):
		line=File.readline().strip('\n')
	segsites=int(File.readline().strip('\n').split(':')[1])
	locs=map(float, File.readline().strip('\n').split(':')[1].split(' ')[1:-1] )
	chrm=[]
	for x in range (0, len(locs) ):
		locs[x]=int(round(locs[x]*size))
	for line in File.readlines():
		chrm.append( chromosome(size) )
		for x in range (0, len(line) ):
			if(line[x]=="1"):
				chrm[-1].insert([locs[x]])
	return chrm

def read_ms_fake (File):
	line=""
	line=File.readline().strip('\n')
	size=int(line.strip('\n').split('-')[-1].split(' ')[2])
	while (line!="//"):
		line=File.readline().strip('\n')
	segsites=int(File.readline().strip('\n').split(':')[1])
	locs=map(float, File.readline().strip('\n').split(':')[1].split(' ')[1:-1] )
	chrm=[]
	for x in range (0, len(locs) ):
		locs[x]=int(round(locs[x]*size))

	for line in File.readlines():
		chrm.append( chromosome(size) )
		for x in range (0, len(locs) ):
			#if(line[x]=="1"):
			if (len(chrm)%2==0):
				chrm[-1].insert([locs[x]])
	return chrm

def print_ms (chrms, header):
	locs=[]
	#print chrms[0].poly
	size=chrms[0].size
	for chrm in chrms:
		locs = sorted(set(locs+chrm.poly))	
	print "none "+str(len(chrms))+" "+str(header)+"\n\n\n//"
	print "segsites: "+str(len(locs) )

	segsites=locs[:]
	for x in range (0, len(segsites) ):
		segsites[x]=float(segsites[x])/float(chrms[0].size)
	print "positions: "+' '.join(map(str, segsites) )
	for y in range (0, len(chrms) ):
		out=[]
	#	print chrms[y].poly
		for p in locs:
			if p not in chrms[y].poly:
				out.append('0')
			else:
				out.append('1')
		print ''.join(out)
def mutate (chrm, rate):
	x=numpy.random.binomial(chrm.size, rate)
	switch=numpy.random.randint(0,chrm.size,x)
	chrm.insert(switch)
	return chrm

def recombine (chrm1, chrm2, rate):
	if chrm1.size!=chrm2.size:
		print "Cannot handle different sized chromosomes..."
	x=numpy.random.binomial(chrm1.size, rate)	
	switch=numpy.random.randint(0, chrm1.size, x)
	switch.sort()
	switch=numpy.append(switch, chrm1.size)
	chrm3=chromosome(chrm1.size)
	even=False
	chrm3.poly=chrm1.poly[:]
	for s in switch:
		if (even):
			chrm3.splice(chrm1, s)
			even= ~even
		else:
			chrm3.splice(chrm2, s)
			even= ~even
	return chrm3

File=open(sys.argv[1])

chrm=read_ms(File)
mu = u
c  = 0.003125

#pre=read_ms_fake(File)
#chrm=[]
#N=len(pre)
#for x in range(0, 3200, 2):
#	chrm.append( recombine(pre[0], pre[1], c) )

N=len(chrm)
if (N%2==1):
	N=N-1
#print N
END=0 #N/8
O1=[]

for x in range(0, END*7, 2):
	O1.append( recombine(chrm[x], chrm[x+1], c) )
	O1[-1]=mutate(O1[-1], mu)
#	O1.append( recombine(chrm[x], chrm[x+1], c) )
#	O1[-1]=mutate(O1[-1], mu)

I=[]

for x in range(END*7, N, 2):
	I.append( recombine(chrm[x], chrm[x+1], c) )
	I[-1]=mutate(I[-1], mu)
	I.append( recombine(chrm[x], chrm[x+1], c) )
	I[-1]=mutate(I[-1], mu)

print_ms(I, len(I) )

quit()

O2=[]

for x in range(0, len(I), 2):
	O2.append( recombine( I[x], I[x+1], c) )
	O2[-1]=mutate(O2[-1], mu)
	O2.append( recombine( O1[x], O1[x+1], c) )
	O2[-1]=mutate(O2[-1], mu)

if len(O1)%2==1:
	O1.pop()

if len(I)%2==0:
	for x in range(len(I), len(O1), 2):
		O2.append( recombine(O1[x], O1[x+1], c) )
		O2[-1]=mutate(O2[-1], mu)
else :
	for x in range(len(I)+1, len(O1), 2):
		O2.append( recombine(O1[x], O1[x+1], c) )
		O2[-1]=mutate(O2[-1], mu)

print_ms(O2, len(I) )
