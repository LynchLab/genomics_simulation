import argparse

parser = argparse.ArgumentParser(description='make a mutation file.')
parser.add_argument('-N', metavar='--number', type=int, default=0,
                        help='column')
parser.add_argument('-s', metavar='--state', type=str, default="states.txt",
                        help='state file')
parser.add_argument('-m', metavar='--mutation', type=str, default="mutation.txt",
                        help='mutation file')

args = parser.parse_args()

state_file=open(args.s)
mutat_file=open(args.m)
N=args.N

states=[]

File_array_0=[]
File_array_1=[]

for n in range (0, N):
	N_str='%03d' % n
	#print N_str
	File_array_0.append(open("seq_"+N_str+".0.poly", "w") )
	File_array_1.append(open("seq_"+N_str+".1.poly", "w") )

for line in state_file:
#	states.append( line.strip('\n').split('\t')[1:]  ) 
	states.append( line[1:]  ) 

for line in mutat_file:
	state=states.pop(0)
	#print "*", state
	for n in range (0, N):
		#print "HI!", n, N
#		print "*", state[n*4], state[n*4]=='1', state[n*4+2], state[n*4+2]=='1'
		if state[n*4]=='1':
			File_array_0[n].write(line)
		if state[n*4+2]=='1':
			File_array_1[n].write(line)
