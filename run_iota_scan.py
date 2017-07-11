#from subprocess import call
import os
import numpy
os.system("some_command < input_file | another_command > output_file")  


for x in range(0, 10):
	os.system("zcat ~/temp/genomics_simulation/sequences/states.txt.gz | python-2.7.9 states_to_pheno_w_inbred.py 500 0 10000 300 ~/temp/genomics_simulation/analysis_pipelines/inbred -0.9 > temp")


	File=open("temp")
	a=[]
	d=[]
	i=[]
	z=[]

	for line in File:
		line=line.split()
		z.append(float(line[2]) )
		a.append(float(line[3]) )
		d.append(float(line[4]) )
		i.append(float(line[5]) )

	z_var=numpy.var(z)
	a_var=numpy.var(a)
	d_var=numpy.var(d)
	i_var=numpy.var(i)
	print x, z_var, a_var, d_var, i_var
