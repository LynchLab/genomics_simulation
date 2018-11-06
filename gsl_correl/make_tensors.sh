echo "#include <cstdint>" > $1
echo "int calc(const double *a_vals, const double *b_vals, double *q_vals, const int32_t &B, const int32_t &N){" >> $1

echo "uint32_t a1_dimension=B;" >> $1
echo "uint32_t a2_dimension=N;" >> $1

echo "uint32_t b1_dimension=B;" >> $1
echo "uint32_t b2_dimension=N;" >> $1

echo "double *ir_vals=new double[B*B*N];" >> $1
echo "double *ic_vals=new double[B*B*N];" >> $1
echo "double  *i_vals=new double[B*B*N];" >> $1

echo "int32_t ir1_dimension=B;" >> $1
echo "int32_t ir2_dimension=B;" >> $1
echo "int32_t ir3_dimension=N;" >> $1

echo "int32_t ic1_dimension=B;" >> $1
echo "int32_t ic2_dimension=B;" >> $1
echo "int32_t ic3_dimension=N;" >> $1

echo "int32_t i1_dimension=B;" >> $1
echo "int32_t i2_dimension=B;" >> $1
echo "int32_t i3_dimension=N;" >> $1

echo "int32_t q1_dimension=N;" >> $1
echo "int32_t q2_dimension=N;" >> $1

taco "ir(k,l,i) = a(k,i)*b(l,i)" #| sed 's/\x1b\[[0-9;]*[a-zA-Z]//g' >> $1
taco "ic(k,l) = ir(k,l,j)" #| sed 's/\x1b\[[0-9;]*[a-zA-Z]//g' >> $1 
taco "i(k,l,i) = ir(k,l,i)-ic(k,l)" #| sed 's/\x1b\[[0-9;]*[a-zA-Z]//g' >> $1 
taco "q(i,j) = i(k,l,i)*i(k,l,j)" #| sed 's/\x1b\[[0-9;]*[a-zA-Z]//g' >> $1

echo "delete [] ir_vals;" >> $1
echo "delete [] ic_vals;" >> $1
echo "delete [] i_vals;" >> $1
echo "}" >> $1
