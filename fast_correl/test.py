import math

def row_index(i, M):
    ii = M*(M+1)/2-1-i
    K = math.floor((math.sqrt(8*ii+1)-1)/2)
    return int(M-1-K)

def column_index(i, M):
    ii = M*(M+1)/2-1-i
    K = math.floor((math.sqrt(8*ii+1)-1)/2)
    return int( (i - M*(M+1)/2 + (K+1)*(K+2)/2 +M-1-K) )
    #return int( (i + (K+1)*(K+2)/2 ) )

def get_i(x,y, M):
	return	(M*(M+1)/2) - (M-x)*((M-x)+1)/2 + y - x 

print get_i(0,0,7)
print get_i(0,6,7)
print get_i(3,3,7)
print get_i(3,6,7)
print get_i(6,6,7)

#print row_index(0,4)
#print row_index(3,4)
#print row_index(4,4)
#print row_index(6,4)
#print row_index(7,4)
#print row_index(8,4)
#print row_index(9,4)

#print column_index(0,4)
#print column_index(3,4)
#print column_index(4,4)
#print column_index(6,4)
#print column_index(7,4)
#print column_index(8,4)
#print column_index(9,4)
