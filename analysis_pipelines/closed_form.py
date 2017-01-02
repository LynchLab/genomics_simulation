import math
import sys
from scipy import stats
import numpy
		
File=open(sys.argv[1])
line=File.readline()
X=map(float, line.split('\t')[2:])
File.close()

print "run freq rho1 rho2 rho3 P OX EX OY EY"

File=open(sys.argv[1])

for y in range(0, len(X) ):
	for x in range (y+1, len(X) ):
		data=[]
		fd={}
		for line in File:
			X=map(float, line.split('\t')[2:])
			freq=sum(X)-X[x]-X[y]
			try:
				fd[freq].append(X)
			except:
				fd[freq]=[X]
			data.append(X)
		File.close()
		sum_x={}
		sum_y={}
		freq_array=[]
		x_array=[]
		y_array=[]

		for freq in fd.keys():
			SUM_X=0
			SUM_Y=0
			N=float(len(fd[freq]))
			if N==0:
				continue

			for datum in fd[freq]:
				SUM_X+=datum[x]
				SUM_Y+=datum[y]

			sum_x[freq]=float(SUM_X)/float(2*N)
			sum_y[freq]=float(SUM_Y)/float(2*N)
	
			freq_array.append(float(freq)/float((len(X)-2)*2) )
			x_array.append(float(SUM_X)/float(2*N) )
			y_array.append(float(SUM_Y)/float(2*N) )

		lm_x=stats.linregress(freq_array[1:-2], x_array[1:-2])
		lm_y=stats.linregress(freq_array[1:-2], y_array[1:-2])

		lm_x2=numpy.polyfit( freq_array[1:-2], x_array[1:-2], 5)
		lm_y2=numpy.polyfit( freq_array[1:-2], y_array[1:-2], 5)

		index=0
		for freq in fd.keys():

			print y*len(X)+x,

			VAR_X=0
			VAR_Y=0
			COV_XY=0
			COV_XY2=0
			COV_XY3=0

			N=float(len(fd[freq]))

			P=float(freq)/float((len(X)-2)*2)
			Q=1-P

			P_X=lm_x[0]*P+lm_x[1]
			P_Y=lm_y[0]*P+lm_y[1]

#			P_X=lm_x2[0]*P**3+lm_x2[1]*P**2+lm_x2[2]*P+lm_x2[3]
#			P_Y=lm_y2[0]*P**3+lm_y2[1]*P**2+lm_y2[2]*P+lm_y2[3]
			P_X=lm_x2[0]*P**5+lm_x2[1]*P**4+lm_x2[2]*P**3+lm_x2[3]*P**2+lm_x2[4]*P+lm_x2[5]
			P_Y=lm_y2[0]*P**5+lm_y2[1]*P**4+lm_y2[2]*P**3+lm_y2[3]*P**2+lm_y2[4]*P+lm_y2[5]

			OX=x_array[index]
			OY=y_array[index]
			for datum in fd[freq]:
			#	VAR_X+=(datum[x]/2.-SUM_X)**2
			#	VAR_Y+=(datum[y]/2.-SUM_Y)**2
			#	COV_XY+=(datum[x]/2.-SUM_X)*(datum[y]/2.-SUM_Y)
				COV_XY+=(float(datum[x])/2.-P_X)*(float(datum[y])/2.-P_Y)
				COV_XY2+=(float(datum[x])/2.-P)*(float(datum[y])/2.-P)
				COV_XY3+=(float(datum[x])/2.-OX)*(float(datum[y])/2.-OY)
			VAR_X=P_X*(1-P_X)
			VAR_Y=P_Y*(1-P_Y)
			VAR_OX=OX*(1-OX)
			VAR_OY=OY*(1-OY)
			if VAR_X==0 or VAR_Y==0:
				continue
			COV_XY/=N
			COV_XY2/=N
			COV_XY3/=N
			try:
				index+=1
				print freq, COV_XY/math.sqrt(VAR_X*VAR_Y), COV_XY2/P/Q, COV_XY3/math.sqrt(VAR_OX*VAR_OY), P, P_X, x_array[index-1], P_Y, y_array[index-1]
			except:
				print "NaN"
			
