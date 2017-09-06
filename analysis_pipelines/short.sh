A=`cat out.txt`
B=`cat ../analysis_files/test.hsq | grep "V(.)	" | cut -f 2 -d '	'| tr "\n" " "`
echo $A $B
