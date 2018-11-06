import sys

t=["filename"]
v=[sys.argv[5]]

File=open(sys.argv[1])

#print (len(sys.argv))
for line in File:
	line=line.strip('\n').split('\t')
	if len(line) > 1:
		v.append(line[1])
		t.append(line[0])

File=open(sys.argv[2])
for line in File:
	if line[0:7]=="V(G)/Vp":
		line=line.strip('\n').split('\t')
		v.append(line[1])
		t.append("V(G)/Vp(gcta)")

File=open(sys.argv[3])
for line in File:
	if line[0:7]=="V(G)/Vp":
		line=line.strip('\n').split('\t')
		v.append(line[1])
		t.append("V(G)/Vp(gcta2)")


File=open(sys.argv[4])
for line in File:
	if line[0:3]=="Var" or line[0:5]=="2*Cov" :
		line=line.strip('\n').split('\t')
		v.append(line[1])
		t.append(line[0])

File=open(sys.argv[5])
for line in File:
	if line[0:6]=="vars 1":
		line=line.split()
		v+=line[-6:]
		t.append("V_MtM")
		t.append("V_MJM")
		t.append("V_HtH")
		t.append("V_HJH")
		t.append("V_MtH")
		t.append("V_E")

if ( len(sys.argv)==7 ):
	for T in t:
		print T,
	print

for V in v:
	print V,
print
