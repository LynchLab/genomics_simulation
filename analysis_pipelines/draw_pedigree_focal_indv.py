import sys

def get_color(z):
	global min_z, max_z
	if max_z!=min_z:
		value=float(z-min_z)/float(max_z-min_z)
	else :
		value=1.
	G="{:02x}".format(int(value*255))
	return "\"#ff"+G+G+"\""

	return "blue"

def set_parent (ind, depth):
	global individuals, File, linesize, pad, start, max_z, min_z
	if not (ind.draw) :
		ind.draw=True
		if depth > 0:
			File.seek((int(ind.pid)-start)*linesize+pad)
			line=File.readline()
			time, id, pid, mid, sex, z=line.strip('\n').split()
			if id not in individuals.keys():
				individuals[id]=Individual(id, pid, mid, sex, float(z))
			max_z=max(max_z, float(z) )
			min_z=min(min_z, float(z) )
			set_parent(individuals[ind.pid], depth-1)

			File.seek((int(ind.mid)-start)*linesize+pad)
			line=File.readline()
			time, id, pid, mid, sex, z=line.strip('\n').split()
			if id not in individuals.keys():
				individuals[id]=Individual(id, pid, mid, sex, float(z))
			max_z=max(max_z, float(z) )
			min_z=min(min_z, float(z) )
			set_parent(individuals[ind.mid], depth-1)

class Individual:
	def __init__ (self, id, pid, mid, sex, z):
		self.id=id
		self.pid=pid
		self.mid=mid
		self.sex=sex
		self.z=z
		self.draw=False
	def node(self):
		return "\t\t"+self.id+" [shape=circle, fillcolor="+get_color(self.z)+", fixedsize=true, label=\"\", style=filled]\n"
	def edges(self):
		return "\t"+self.pid+" -> "+self.id+"\n\t"+self.mid+" -> "+self.id+'\n'

File=open(sys.argv[1])

Ind1=(sys.argv[2])

DEAPTH=int(sys.argv[3])

individuals={}

line=File.readline()
pad=len(line)
start=int(File.readline().strip('\n').split()[1])

min_z=10000
max_z=-10000

linesize=61
File.seek((int(Ind1)-start)*linesize+pad)
line=File.readline()
time, id, pid, mid, sex, z=line.strip('\n').split()
individuals[id]=Individual(id, pid, mid, sex, float(z))
max_z=max(max_z, float(z) )
min_z=min(min_z, float(z) )

set_parent(individuals[Ind1], DEAPTH)

print "digraph G {\n",
print "\tnode [shape=circle, label=\"\"]\n\t{"
#, label=\"\"]"
for ind in individuals.keys():
	this_ind=individuals[ind]
	if (this_ind.draw):
		print this_ind.node(),
print "\t}" 
for ind in individuals.keys():
	this_ind=individuals[ind]
	if (this_ind.draw):
		try:
			if (individuals[this_ind.pid].draw and individuals[this_ind.mid].draw ):
				print this_ind.edges(),
		except:
			E=0
print "}"
