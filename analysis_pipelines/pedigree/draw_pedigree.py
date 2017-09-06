import sys

def get_color(z):
	global min_z, max_z
	value=float(z-min_z)/float(max_z-min_z)
	G="{:02x}".format(int(value*255))
	return "\"#ff"+G+G+"\""

	return "blue"

class Individual:
	def __init__ (self, id, pid, mid, sex, z):
		self.id=id
		self.pid=pid
		self.mid=mid
		self.sex=sex
		self.z=z
	def node(self):
		return "\t\t"+self.id+" [shape=circle, fillcolor="+get_color(self.z)+", fixedsize=true, label=\"\", style=filled]\n"
	def edges(self):
		ret=""
		if self.pid!="0":
			ret="\t"+self.pid+" -> "+self.id+"\n"
		if self.mid!="0":
			ret=ret+"\t"+self.mid+" -> "+self.id+'\n'
		return ret

File=open(sys.argv[1])
individuals={}
File.readline()

min_z=10000
max_z=-10000

for line in File:
	time, id, pid, mid, sex, z=line.strip('\n').split('\t')
	individuals[id]=Individual(id, pid, mid, sex, float(z))
	max_z=max(max_z, float(z))
	min_z=min(min_z, float(z))
print "digraph G {\n",
print "\tnode [shape=circle, label=\"\"]\n\t{"
#, label=\"\"]"
for ind in individuals.keys():
	this_ind=individuals[ind]
	print this_ind.node(),
print "\t}" 
for ind in individuals.keys():
	this_ind=individuals[ind]
	print this_ind.edges(),
print "}"
