import sys
import mappy

print "@NAME:SCAFFOLDS	VERSION:0.4.27-2016-12-20	FORMAT:TEXT	CONCATENATED"
print "@NAME          	LENGTH"
print "1	230218"
print "@END_TABLE"
print "@NAME:GENOTYPES	VERSION:0.4.27-2016-12-20	FORMAT:TEXT	CONCATENATED	INDEXED"
print "@SCFNAME	POS	MN_FREQ	F_STAT	",

lookup={ 0:{0:"1/0/0/30", 1:"0/1/0/30"}, 1:{0:"0/1/0/30", 1:"0/0/1/30"} }
def samples_to_string(samples):
	string=[]
	N=len(samples)/2
	for x in range(0, len(samples)/2):
		string.append(lookup[samples[x]][samples[x+N]])
	return '\t'.join(string)

class vcf_file:
	def __init__ (self):
		self.N=0
		self.chrom="1"
		self.pos=0
		self.AF=0
		self.info=str(self.AF)+"\t0"
		self.samples=[]
	def print_line(self):
		self.pos+=1
		print '\t'.join(map(str, [self.chrom, self.pos, self.info, samples_to_string(self.samples)]))
	def set_sample(self, N, A, B):
		self.samples[N*2]=int(A)
		self.samples[N*2+1]=int(B)
	def update_info(self):
		self.info=str(self.AF)+"\t0"
	def update_af(self):
		self.AF=float(sum(self.samples))/float(self.N*2)
		self.update_info()
	def set_sample_size(self, N):
		self.samples=[0]*N*2
		self.N=N

def read_state(field, vcf, N):
	mask=[0x0000001, 0x00000002, 0x00000004, 0x00000008,
	      0x0000010, 0x00000020, 0x00000040, 0x00000080,
	      0x0000100, 0x00000200, 0x00000400, 0x00000800,
	      0x0001000, 0x00002000, 0x00004000, 0x00008000,
	      0x0010000, 0x00020000, 0x00040000, 0x00080000,
	      0x0010000, 0x00200000, 0x00400000, 0x00800000,
	      0x0100000, 0x02000000, 0x04000000, 0x08000000,
	      0x1000000, 0x20000000, 0x40000000, 0x80000000]
	for b in range(0, 32):
		for x in range(0, N):
			vcf.set_sample(x, field[x*2] & mask[b] !=0, field[x*2+1] & mask[b]!=0)
		vcf.update_af()
		vcf.print_line()

namefile=mappy.open(sys.argv[1])
name=mappy.read(namefile).split('\t')[1:]

N=len(name)

field=[0]*N*2
vcf=vcf_file()
vcf.set_sample_size(N)

print '\t'.join(name)

while (True):
	bfield=sys.stdin.read(N*8)
	if bfield == '':
        	break
	for x in range(0, N*2):
		field[x]=(ord(bfield[x*4+0])<<24)+(ord(bfield[x*4+1])<<16)+(ord(bfield[x*4+2])<<8)+ord(bfield[x*4+3])
	read_state(field, vcf, N)
print "@END_TABLE"
