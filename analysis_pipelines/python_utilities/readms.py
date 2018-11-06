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

def get_segsites (File):
	line=""
	line=File.readline().strip('\n')
	#size=int(line.strip('\n').split('-')[-1].split(' ')[2])
	while (line!="//"):
		line=File.readline().strip('\n')
	line=File.readline().strip('\n')
	locs=map(float, File.readline().strip('\n').split(':')[1].split(' ')[1:-1] )
	for x in range (0, len(locs) ):
		locs[x]=(locs[x]*100)
	return locs
