import sys

def samples_to_string(samples):
	string=[]
	for x in range(0, len(samples)/2):
		string.append('|'.join(map(str, [samples[x*2], samples[x*2+1]])))
	return '\t'.join(string)

class state_file:
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


def write_state(lines, state, N):
	state.clear()
	for b in range(0, 32):
		field=lines[b]
		for x in range(0, N):
			state.set_sample(b, x, field[x*2],  field[x*2+1])
	state.print_line()

state=state_file()
K=[""]*32
N=0
while (True):
	L=[]
	while len(L)<32:
		line=sys.stdin.readline()
		if line=='':
			quit()
		if line[0]=='#':
			continue
		else:
			L.append(line)
	for x in range(0, 32):
		line=L[x].split('\t')
		if N==0:
			N=(len(line)-1)/2
			state.set_sample_size(N)
		K[x]=[0]*N*2
		K[x]=map(int, line[1:])
	write_state(K, state, N)
