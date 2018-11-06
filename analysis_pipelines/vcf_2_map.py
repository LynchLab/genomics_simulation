import sys

MapFile=open(sys.argv[1])
VcfFile=open(sys.argv[2])

chrmM, posM, cM=MapFile.readline().strip('\n').split('\t')
chrmM, posM, cM=MapFile.readline().strip('\n').split('\t')

cL=0
posL=1

for line in VcfFile:
	if line[0]!='#':
		line=line.strip('\n').split('\t')

		chrm1=line[0]
		pos1=int(line[1])

		if chrm1==chrmM:
			while (pos1>=posM):
				chrmL, posL, cL=chrmM, cM, posM
				posM="NA"
				while ("NA"==posM):
					chrmM, posM, cM=MapFile.readline().strip('\n').split('\t')
				
		
			if pos1<posM:
				#linear interpolation.
				cL=float(cL)
				cM=float(cM)
				posL=int(posL)
				posM=int(posM)
				print chrm1, pos1, cL+(cM-cL)*(pos1-posL)/(posM-posL)
