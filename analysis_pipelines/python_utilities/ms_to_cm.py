import	sys
import	readms

File=open(sys.argv[1])
size=float(sys.argv[2])
segsites=readms.get_segsites(File)

print "@NAME:SCAFFOLDS	VERSION:0.4.21	FORMAT:TEXT	CONCATENATED"
print "@NAME	LENGTH"
print "1	"+str(len(segsites)+1)
print "@END_TABLE"
print "@NAME:DATA	VERSION:TYPED	FORMAT:TEXT	CONCATENATED	INDEXED"
print "@SCFNAME	POS	VALUE"

for x in range(1, len(segsites)+1 ):
	print	"1\t"+str(x)+"\t"+str(segsites[x-1]*size*2)
print "@END_TABLE"
