import sys

def samples_to_string(samples):
	string=[]
	for x in range(0, len(samples)/2):
		string.append('|'.join(map(str, [samples[x*2], samples[x*2+1]])))
	return '\t'.join(string)

class vcf_file:
	def __init__ (self):
		self.samples=""
		self.N=0
	def print_line(self):
		for x in self.samples:
			sys.stdout.write(chr(x))
	def set_sample(self, b, N, A, B):
		rr=b%8
		rl=int(b/8)
		mask=[0x0000001, 0x00000002, 0x00000004, 0x00000008,
		      0x0000010, 0x00000020, 0x00000040, 0x00000080]

#		      0x0000100, 0x00000200, 0x00000400, 0x00000800,
#		      0x0001000, 0x00002000, 0x00004000, 0x00008000,
#		      0x0010000, 0x00020000, 0x00040000, 0x00080000,
#		      0x0010000, 0x00200000, 0x00400000, 0x00800000,
#		      0x0100000, 0x02000000, 0x04000000, 0x08000000,
#		      0x1000000, 0x20000000, 0x40000000, 0x80000000]
		if A:
			self.samples[N*4+rl]+=mask[rr] 
		if B:
			self.samples[self.N*4+N*4+rl]+=mask[rr]
	def set_sample_size(self, N):
		self.zero=[0]*4*N*2
		self.samples=self.zero[:]
		self.N=N
	def clear(self):
		self.samples=self.zero[:]


def write_state(lines, vcf, N):
	vcf.clear()
	for b in range(0, 32):
		field=lines[b]
		for x in range(0, N):
			vcf.set_sample(b, x, field[x*2],  field[x*2+1])
	vcf.print_line()

vcf=vcf_file()
K=[""]*32
N=0
while (True):
	L=[]
	while len(L)<32:
		line=sys.stdin.readline()
		if line[0]=='#':
			continue
		else:
			L.append(line)
	for x in range(0, 32):
		line=L[x].split('\t')
		if N==0:
			N=len(line)-9
			sys.stderr.write(str(N)+'\n' )
			vcf.set_sample_size(N)
		K[x]=[0]*N*2
		for y in range(0, N):
			K[x][y*2]=int(line[9+y].split(':')[0].split('/')[0]=='1')
			K[x][y*2+1]=int(line[9+y].split(':')[0].split('/')[1]=='1')
		#	K[x][y*2]=int(line[9+y].split(':')[0].split('|')[0]=='1')
		#	K[x][y*2+1]=int(line[9+y].split(':')[0].split('|')[1]=='1')
	write_state(K, vcf, N)
