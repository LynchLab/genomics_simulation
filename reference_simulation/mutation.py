import random
import sys
import argparse

base=['A','C','G','T']
btoi={'A':0,'C':1,'G':2,'T':3,'a':0,'c':1,'g':2,'t':3}

def denovo (string, start, stop, *arg):
	new=[]
	for x in range(start, start+stop):
		new.append(base[random.randint(0,3)])
	return string[:start]+''.join(new)+string[start:]
def denovo_size (size, start, stop, *arg):
	return size+stop

def copy (string, src_start, src_stop, dst_start, forward, *arg):
	if (forward):
		return string[:dst_start]+string[src_start:src_stop]+string[dst_start:]
	else:
		return string[:dst_start]+string[src_stop:src_start:-1]+string[dst_start:]
def copy_size (size, src_start, src_stop, dst_start, forward, *arg):
	return size+min(src_stop, size)-min(src_start, size)

def tdcopy (string, src_start, src_stop, dst_start, forward, n):
	if (forward):
		return string[:dst_start]+string[src_start:src_stop]*n+string[dst_start:]
	else:
		return string[:dst_start]+string[src_stop:src_start:-1]*n+string[dst_start:]
def tdcopy_size (size, src_start, src_stop, dst_start, forward, n):
	return size+(min(src_stop, size)-min(src_start, size))*n

def trans (string, src_start, src_stop, dst_start, forward, *arg):
	if (dst_start<src_start):
		if (forward):
			return string[:dst_start]+string[src_start:src_stop]+string[dst_start:src_start]+string[src_stop:]	
		else:
			return string[:dst_start]+string[src_stop:src_start:-1]+string[dst_start:src_start+1]+string[src_stop:]	
	elif dst_start<src_stop:
		return string
	else:
		if (forward):
			return string[:src_start]+string[src_stop:dst_start]+string[src_start:src_stop]+string[dst_start:]	
		else:
			return string[:src_start]+string[src_stop:dst_start+1]+string[src_stop:src_start:-1]+string[dst_start:]	
def trans_size (size, src_start, src_stop, dst_start, forward, *arg):
	return size

def invert (string, start, stop, *arg):
	return string[:start+1]+string[stop:start:-1]+string[stop:]	
def invert_size (size, start, stop, *arg):
	return size

def delete (string, start, stop, *arg):
	return string[:start]+string[stop:]
def delete_size (size, start, stop, *arg):
	return size-(min(stop, size)-min(start, size) )
	
def static (string, start, stop, *arg):
	return string[:start]+string[start:stop].upper()+string[stop:]
def static_size (size, start, stop, *arg):
	return size
	
def SNP(string, locus, newbase, *arg):
	before=len(string)
	if locus==0:
		string=base[newbase]+string[1:]
	elif(locus<len(string) ):
		oldbase=btoi[string[locus]]
		newbase=(oldbase+newbase)%4
		string=string[:locus]+base[newbase]+string[locus+1:]
	return string
def SNP_size(size, locus, newbase, *arg):
	return size

def SUB(char, newbase, *arg):
	oldbase=btoi[char]
	newbase=(oldbase+newbase)%4
	return base[newbase]
def SUB_size(size, locus, newbase, *arg):
	return size

fns=[denovo, 		copy, 		tdcopy,		trans, 		invert, 	delete, 	static,		SNP, SUB]
fns_size=[denovo_size, 	copy_size, 	tdcopy_size,	trans_size, 	invert_size, 	delete_size, 	static_size, 	SNP_size, SUB_size]
fnsk={"denovo":denovo, 	"copy":copy, 	"tdcopy":tdcopy, "trans":trans, "invert":invert, "delete":delete, "static":static, "SNP":SNP, "SUB":SUB}
