import sys

loci_file=open(sys.argv[1])
states_file=open(sys.argv[2])
N=int(sys.argv[3])
vcf_file=open(sys.argv[4])

loci=[]
for line in loci_file:
	loci.append(line.split(' ')[1].strip('\n'))

state={}

this_state=[0]*N

for line in states_file:
	states=line.split('\t')[1:]
	name=loci.pop(0)

	for x in range(0, N):
		this_state[x]=int(states[x*2])+int(states[x*2+1])
#	print name, this_state
#	quit()
	state[name]=this_state[:]

QQ={}
for line in vcf_file:
	if line[0]=='#':
		continue
	line=line.split('\t')
	trips=line[9:]
	try:
		#PL=line[8].split(':').index("PL")
		PL=line[8].split(':').index("GP")
	except:
		try:
			PL=line[8].split(':').index("GP")
		except:
			continue
		continue
	name=line[1]
	if (line[4]!='.'):
#		print name, trips[0]
		for x in range(0, N):
			if trips[x][0]!='.':
				t=trips[x].split(':')[PL].split(',')
				try:
					this_state=state[name][x]
				except:
					continue
				for y in range(0, 3):
					Q=int(float(t[y].strip('\n')))
					try:
					#	if( int(Q)>20 and this_state==y):
					#		print name, x, state[name], trips[x]
						QQ[Q][0]+=(this_state==y)
						QQ[Q][1]+=1
					except:
						try:
							QQ[Q]=[0,0]
							QQ[Q][0]+=(this_state==y)
							QQ[Q][1]+=1
						except:
							E=0
for Q in QQ.keys():
	print "A:", Q, QQ[Q][0], QQ[Q][1]
