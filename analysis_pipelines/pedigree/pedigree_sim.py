import sys
import random

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
		self.z=-1
		if self.pid=="0" and self.mid=="0":
			self.z=random.randint(0,2)
	def node(self):
		return "\t\t"+self.id+" [shape=circle, fillcolor="+get_color(self.z)+", fixedsize=true, label=\"\", style=filled]\n"
	def update(self):
		global individuals
		if self.z==-1:
			self.z=0
			if self.pid!="0":
				self.z+=individuals[self.pid].get()
			else:
				self.z+=random.randint(0,1)
			if self.mid=="0":
				self.z+=individuals[self.pid].get()
			else:
				self.z+=random.randint(0,1)
	def get(self):
		if self.z==-1:
			self.update()
		if self.z==0:
			return 0
		elif self.z==2:
			return 1
		else :
			return random.randint(0,1)
			
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

min_z=0
max_z=2

for line in File:
	time, id, pid, mid, sex, z=line.strip('\n').split('\t')
	individuals[id]=Individual(id, pid, mid, sex, float(z))

for ind in individuals:
	individuals[ind].update()

print "digraph G {\n",
print "bgcolor=\"transparent\""
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
