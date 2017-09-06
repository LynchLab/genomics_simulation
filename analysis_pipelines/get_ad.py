import sys
import numpy
import math

File=open(sys.argv[1])
File=File.readlines()

#N=int( (math.sqrt(8*len(File))+1)/2 )
N=int( (math.sqrt(8*len(File))-1)/2 )

#MU=numpy.zeros((N, N))
RS=numpy.zeros((N, N))
#FS=numpy.zeros((N, N))

for z in range(1, N*(N-1) /2+1):
	line=File[z].split(' ')

	x=int(line[0])
	y=int(line[1])

#        MU[x,y]=float(line[5])
#        MU[y,x]=MU[x,y]

        RS[x,y]=float(line[4])
        RS[y,x]=RS[x,y]

#        FS[x,y]=float(line[6])
#        FS[y,x]=FS[x,y]

for x in range(0, N):
	for y in range(0, N):
		print RS[x,y],
	print
