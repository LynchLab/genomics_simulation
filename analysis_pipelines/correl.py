import math
import sys
from scipy import stats
import numpy
import itertools

if len(sys.argv)==2:
	print sys.argv	
	File=open(sys.argv[1])
	File=File.read().split('\n')[:-1]
else:
	File=sys.stdin
	File=File.read().split('\n')[:-1]


print "run freq rho1 rho2 rho3 P OX EX OY EY"
for j in range(0, len(File) ):
	File[j]=map(int, File[j].split('\t')[1:])

N=len(File[0])

for y in range(0, N/2 ):
	for x in range (y+1, N/2 ):
		fd={}
		for j in range(0, len(File) ):
			D=File[j]
			freq=sum(D)-D[x*2]-D[x*2+1]-D[y*2]-D[y*2+1]
			X=D[x*2]+D[x*2+1]
			Y=D[y*2]+D[y*2+1]
			try:
				fd[freq].append([X,Y])
			except:
				fd[freq]=[[X,Y]]

		sum_x={}
		sum_y={}
		freq_array=[]
		x_array=[]
		y_array=[]

		for freq in fd.keys():
			SUM_X=0
			SUM_Y=0
			n=float(len(fd[freq]))
			if n==0:
				continue

			for datum in fd[freq]:
				SUM_X+=datum[0]
				SUM_Y+=datum[1]

			sum_x[freq]=float(SUM_X)/float(n*2)
			sum_y[freq]=float(SUM_Y)/float(n*2)
	
			freq_array.append(float(freq)/float(N-4) )

			x_array.append(float(SUM_X)/float(n*2) )
			y_array.append(float(SUM_Y)/float(n*2) )

#			print (float(SUM_X)/float(N*2) )
#			print (float(SUM_Y)/float(N*2) )

		lm_x=stats.linregress(freq_array[1:-2], x_array[1:-2])
		lm_y=stats.linregress(freq_array[1:-2], y_array[1:-2])

		lm_x2=numpy.polyfit( freq_array[1:-2], x_array[1:-2], 5)
		lm_y2=numpy.polyfit( freq_array[1:-2], y_array[1:-2], 5)

		index=0
		for freq in fd.keys():


			VAR_X=0
			VAR_Y=0
			COV_XY=0
			COV_XY2=0
			COV_XY3=0

			n=float(len(fd[freq]))

			P=float(freq)/float(N-4)
			Q=1-P

			P_X=lm_x[0]*P+lm_x[1]
			P_Y=lm_y[0]*P+lm_y[1]

#			P_X=lm_x2[0]*P**3+lm_x2[1]*P**2+lm_x2[2]*P+lm_x2[3]
#			P_Y=lm_y2[0]*P**3+lm_y2[1]*P**2+lm_y2[2]*P+lm_y2[3]

#			P_X=lm_x2[0]*P**5+lm_x2[1]*P**4+lm_x2[2]*P**3+lm_x2[3]*P**2+lm_x2[4]*P+lm_x2[5]
#			P_Y=lm_y2[0]*P**5+lm_y2[1]*P**4+lm_y2[2]*P**3+lm_y2[3]*P**2+lm_y2[4]*P+lm_y2[5]

			OX=x_array[index]
			OY=y_array[index]
			for datum in fd[freq]:
			#	VAR_X+=( (datum[x*2]+datum[x*2+1] )/2.-SUM_X)**2
			#	VAR_Y+=( (datum[y*2]+datum[y*2+1] )/2.-SUM_Y)**2
			#	COV_XY+=(datum[x]/2.-SUM_X)*(datum[y]/2.-SUM_Y)
				COV_XY+=(float(datum[0])/2.-P_X)*(float(datum[1])/2.-P_Y)
				COV_XY2+=(float(datum[0])/2.-P)*(float(datum[1])/2.-P)
				COV_XY3+=(float(datum[0])/2.-OX)*(float(datum[1])/2.-OY)
			VAR_X=P_X*(1-P_X)
			VAR_Y=P_Y*(1-P_Y)
			VAR_OX=OX*(1-OX)
			VAR_OY=OY*(1-OY)
			if VAR_X==0 or VAR_Y==0 or VAR_OX==0 or VAR_OY==0 or P==0 or Q==0:
				index+=1
		#		print "skipping ", index, P, OX, OY 
				continue
			COV_XY/=n
			COV_XY2/=n
			COV_XY3/=n
			print y*N/2.+x,
			index+=1
			for pret in [freq, COV_XY2/P/Q, COV_XY/math.sqrt(VAR_X*VAR_Y), COV_XY3/math.sqrt(VAR_OX*VAR_OY), P, P_X, x_array[index-1], P_Y, y_array[index-1] ]:
				print "%0.4f" % pret,
			print
