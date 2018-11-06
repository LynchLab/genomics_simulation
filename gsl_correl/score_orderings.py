import sys
import numpy
import random
import itertools
from scipy import misc

#How many diagonal regions are in the image?
N=int(sys.argv[1])

#Get the eucledian distance between two vectors.
def d(X, Y):
	N=min(len(X), len(Y))
	EXY=0
	#Each element in the vector X and Y has an [R]ed, [B]lue, [G]reen, and [A]lpha values, and 
	#we get the squared difference between the red, blue, and green values of X and Y.
	for i in range(0, N):
		EXY+=(X[i][0]/255.-Y[i][0]/255.)**2
		EXY+=(X[i][1]/255.-Y[i][1]/255.)**2
		EXY+=(X[i][2]/255.-Y[i][2]/255.)**2
	return numpy.sqrt(EXY)/(N-1)

#Just declaring the list which will have the pixels along each edge of the regions.
bottom = [[[] for x in range(N)] for y in range(N)] 
left = [[[] for x in range(N)] for y in range(N)] 
top = [[[] for x in range(N)] for y in range(N)] 
right = [[[] for x in range(N)] for y in range(N)] 

#reading and storing the pixels.
for x in range(0, N):
	for y in range(x, N):
		name="image_"+str(x)+"_"+str(y)+".png"
		
		data = misc.imread(name, flatten = False)

		H=data.shape[0]
		W=data.shape[1]

		bottom[x][y]=data[H-1,0:W]
		left[x][y]=data[0:H,0]

		if (x!=y):
			top[x][y] = data[0,0:W]
			right[x][y] = data[0:H,W-1]

#All possible orderings of regions 0 through N-1.
ALL_ORD=itertools.permutations(range(0, N) )

#Get score for each ordering.
for ORD in ALL_ORD:
	sum_dist=0
	for x in range(0, N):
		for y in range(x, N):
			if (ORD[x] < ORD[y]):
				X0, Y0=ORD[x], ORD[y]
			else:
				X0, Y0=ORD[y], ORD[x]

			if(x!=N-1):
				X1=ORD[x+1]
			if(y!=N-1):
				Y1=ORD[y+1]

			if (x!=y):
				if (X1 < Y0):
					sum_dist+=d(right[X0][Y0], left[X1][Y0])
				else:
					sum_dist+=d(right[X0][Y0], left[Y0][X1])
			if(y!=N-1):
				if (X0 < Y1):
					sum_dist+=d(bottom[X0][Y0], top[X0][Y1])
				else:
					sum_dist+=d(bottom[X0][Y0], top[Y1][X0])
	print ORD, sum_dist
