import sys

class result:
	def __init__ (self, bias, mse):
		self.bias=bias
		self.mse=mse
	def add_times(self, time):
		self.times=' '.join([self.times, time])
	def get_scaling(self):
		time1, time2, time2=self.times.split(' ')
		self.scaling=( (float(time1)/float(time2))+(float(time2)/float(time3)) )/2.

File1=open("../analysis_files/summary.csv")
File2=open("../analysis_files/benchmark.csv")
results={}

for line in File1:
	line=line.replace("\"", "").split(',')
	results[ line[0] ]=result(line[1], line[2])	

for line in File2:
	line=line.split(',')
	results[line[0]].add_time(line[3])

for key in results.keys():
	results[key].get_scaling()
	print key, "&3$\\times$", results[key].bias, "&", results[key].mse, "&", results[key].times, "&", results[key].scaling, "\\\\"
