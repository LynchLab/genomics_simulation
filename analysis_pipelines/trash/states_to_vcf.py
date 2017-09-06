import sys
import mappy

print "##fileformat=VCFv4.2\n\
##fileDate=20090805\n\
##source=myImputationProgramV3.1\n\
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta\n\
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x> ##phasing=partial\n\
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n\
##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n\
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">\n\
##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\"> ##FILTER=<ID=q10,Description=\"Quality below 10\">\n\
##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n\
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">\n\
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	",

def samples_to_string(samples):
	string=[]
	N=len(samples)/2
	for x in range(0, len(samples)/2):
		string.append('|'.join(map(str, [samples[x], samples[x+N]])))
	return '\t'.join(string)

class vcf_file:
	def __init__ (self):
		self.N=0
		self.AA=0
		self.chrom="1"
		self.pos=0
		self.id='.'
		self.ref='T'
		self.alt='A'
		self.qual='0'
		self.filter="pass"
		self.AF=0
		self.info="NS="+str(N)+";DP=30;AA=T;AF="+str(self.AF)
		self.format="GT"
		self.samples=[]
	def print_line(self):
		self.pos+=1
		print '\t'.join(map(str, [self.chrom, self.pos, self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format, samples_to_string(self.samples)]))
	def set_sample(self, N, A, B):
		self.samples[N*2]=int(A)
		self.samples[N*2+1]=int(B)
	def update_info(self):
		self.info="NS="+str(self.N)+";DP=30;AA=T;AF="+str(self.AF)
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
