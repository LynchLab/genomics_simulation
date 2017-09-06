#!/usr/bin/python

import sys
import random

File=open(sys.argv[1])
N=int(sys.argv[2])

f=[]

for line in File:
	f.append(line)

array=list(range(len(f)))

random.shuffle(array)
out=array[1:N]
out=sorted(out)

for x in range(0, N-1):
	print f[out[x]],
