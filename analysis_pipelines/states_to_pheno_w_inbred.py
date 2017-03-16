import sys
import numpy

def samples_to_string(samples):
	string=[]
	for x in range(0, len(samples)/2):
		string.append('|'.join(map(str, [samples[x*2], samples[x*2+1]])))
	return '\t'.join(string)

class pop_file:
	def __init__ (self):
		self.N=0
		self.pos=0
		self.samples=[]
		self.e=[]
		self.a=[]
		self.d=[]
		self.g=[]
	def set_sample(self, a, d, i, z, ind, gen):
		self.a[ind]+=a[gen]
		self.d[ind]+=d[gen]
		self.i[ind]+=i[gen]
		self.g[ind]+=z[gen]
	def set_sample_size(self, N):
		self.N=N
		self.a=[0]*N
		self.d=[0]*N
		self.i=[0]*N
		self.g=[0]*N
		self.e=[0]*N#numpy.random.normal(0,800, N)


def read_state(field, vcf, N):
	global F
	global d
	mask=[0x0000001, 0x00000002, 0x00000004, 0x00000008,
	      0x0000010, 0x00000020, 0x00000040, 0x00000080,
	      0x0000100, 0x00000200, 0x00000400, 0x00000800,
	      0x0001000, 0x00002000, 0x00004000, 0x00008000,
	      0x0010000, 0x00020000, 0x00040000, 0x00080000,
	      0x0010000, 0x00200000, 0x00400000, 0x00800000,
	      0x0100000, 0x02000000, 0x04000000, 0x08000000,
	      0x1000000, 0x20000000, 0x40000000, 0x80000000]
	a=1
	g=[-a,d,a]
	alpha=[-a,0,a]
	delta=[0,d,0]
	iota=[0,d,0]
	genotype=[0]*N
	for b in range(0, 32):
		count=[0,0,0]
		for x in range(0, N):
			G=int(field[x*2] & mask[b] !=0)+int(field[x*2+1] & mask[b]!=0) 
			genotype[x]=G
			count[G]+=1
		count=[float(x) / float(N) for x in count]
		p=count[2]+count[1]/2.
		q=1-p
		if q==0 or p==0:
			continue
		f=1-count[1]/(2.*p*q)
		if f<0:
			continue
		alpha_=a+d*(q-p)
		alpha[0]=-alpha_*2*p
		alpha[1]=alpha_*(q-p)
		alpha[2]=alpha_*2*q

		for x in range(0, N):
			fi=f+F[x]#1.#f#F[x]
			dhwe=0#2*p*q*f*d
			dhwe2=2*p*q*f*d
			dibd=-dhwe2+2*p*q*(f-F[x])*d
			try:
				c=(fi+1)*p*q/((fi*p+q)*(fi*q+p) )
			except:
				print fi, p, q, ((fi*p+q)*(fi*q+p) )
				quit()
			c2=2*fi/(fi+1)
			c3=(q-fi*q) / (fi*p+q)
			c4=(p-fi*p) / (fi*q+p)
			delta[0]=-2*p**2 * d * c3 / c + dhwe
			delta[1]=2*p*q * d / c + dhwe
			delta[2]=-2*q**2 * d * c4 / c + dhwe
			iota[0]=2*p*d*(q-p) * c2 + dibd
			iota[1]=-d*(q-p)**2 * c2 + dibd
			iota[2]=-2*q*d*(q-p) * c2 + dibd
			vcf.set_sample(alpha, delta, iota, g, x, genotype[x] )
N=int(sys.argv[1])
name=int(sys.argv[2])
SIZE=int(sys.argv[3])
LOCI=int(sys.argv[4])
inb=open(sys.argv[5])
d=float(sys.argv[6])
F=[]
for line in inb:
	F.append(float(line) )
loci=sorted(map(int, numpy.random.uniform(0,SIZE, LOCI)))
field=[0]*N*2
pop=pop_file()
pop.set_sample_size(N)
string=[]
locus=loci.pop(0)
L=0
while (True):
	bfield=sys.stdin.read(N*8)
	if bfield == '':
        	break
	if L>locus:
		for x in range(0, N*2):
			field[x]=(ord(bfield[x*4+0])<<24)+(ord(bfield[x*4+1])<<16)+(ord(bfield[x*4+2])<<8)+ord(bfield[x*4+3])
		read_state(field, pop, N)
		if len(loci)>0:
			locus=loci.pop(0)
		else:
			break
	L+=1

a_hat=float(sum(pop.a))/float(N)
g_hat=float(sum(pop.g))/float(N)
d_hat=float(sum(pop.d))/float(N)
i_hat=float(sum(pop.i))/float(N)
e_hat=float(sum(pop.e))/float(N)
z_hat=0#a_hat+d_hat+e_hat+i_hat
for x in range(0, N):
	print '\t'.join(map(str, [x+name, x+name, pop.a[x]+pop.d[x]+pop.e[x]+pop.i[x]-z_hat, pop.a[x], pop.d[x], pop.i[x], pop.e[x]-e_hat, pop.g[x]-g_hat, F[x]]))
