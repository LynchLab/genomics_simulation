import sys

def get_color(z):
	global min_z, max_z
	value=float(z-min_z)/float(max_z-min_z)
	G="{:02x}".format(int(value*255))
	return "\"#ff"+G+G+"\""

	return "blue"

def set_parent (ind, depth):
	global individuals
	if not (ind.draw) :
		ind.draw=True
		if depth > 0:
			try:
				set_parent(individuals[ind.pid], depth-1)
			except:
				E=0
			try:
				set_parent(individuals[ind.mid], depth-1)
			except:
				E=0

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

Ind2=(sys.argv[3])

DEAPTH=int(sys.argv[4])

individuals={}
File.readline()

min_z=10000
max_z=-10000

for line in File:
	time, id, pid, mid, sex, z=line.strip('\n').split('\t')
	individuals[id]=Individual(id, pid, mid, sex, int(z))
	max_z=max(max_z, int(z))
	min_z=min(min_z, int(z))



set_parent(individuals[Ind1], DEAPTH)
set_parent(individuals[Ind2], DEAPTH)

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
		if (individuals[this_ind.pid].draw and individuals[this_ind.mid].draw ):
			print this_ind.edges(),
print "}"
