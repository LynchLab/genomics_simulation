#include <cstdint>
int calc(const double *a_vals, const double *b_vals, double *q_vals, const int32_t &B, const int32_t &N){
uint32_t a1_dimension;
uint32_t a2_dimension;
uint32_t b1_dimension;
uint32_t b2_dimension;
double *ir_vals=new double[B*B*N];
double *ic_vals=new double[B*B*N];
double  *i_vals=new double[B*B*N];
int32_t ir1_dimension=B;
int32_t ir2_dimension=B;
int32_t ir3_dimension=N;
int32_t ic1_dimension=B;
int32_t ic2_dimension=B;
int32_t ic3_dimension=N;
int32_t i1_dimension=B;
int32_t i2_dimension=B;
int32_t i3_dimension=N;
int32_t q1_dimension=N;
int32_t q2_dimension=N;
// Generated by the Tensor Algebra Compiler (tensor-compiler.org)
for (int32_t ka = 0; ka < a1_dimension; ka++) {
  for (int32_t lb = 0; lb < b1_dimension; lb++) {
    int32_t pir2 = ka * ir2_dimension + lb;
    for (int32_t ia = 0; ia < a2_dimension; ia++) {
      int32_t pa2 = ka * a2_dimension + ia;
      int32_t pb2 = lb * b2_dimension + ia;
      int32_t pir3 = pir2 * ir3_dimension + ia;
      ir_vals[pir3] = a_vals[pa2] * b_vals[pb2];
    }
  }
}
// Generated by the Tensor Algebra Compiler (tensor-compiler.org)
for (int32_t kir = 0; kir < ir1_dimension; kir++) {
  for (int32_t lir = 0; lir < ir2_dimension; lir++) {
    int32_t pir2 = kir * ir2_dimension + lir;
    int32_t pic2 = kir * ic2_dimension + lir;
    double tj = 0;
    for (int32_t jir = 0; jir < ir3_dimension; jir++) {
      int32_t pir3 = pir2 * ir3_dimension + jir;
      tj += ir_vals[pir3];
    }
    ic_vals[pic2] = tj;
  }
}
// Generated by the Tensor Algebra Compiler (tensor-compiler.org)
for (int32_t kir = 0; kir < ir1_dimension; kir++) {
  for (int32_t lir = 0; lir < ir2_dimension; lir++) {
    int32_t pir2 = kir * ir2_dimension + lir;
    int32_t pic2 = kir * ic2_dimension + lir;
    int32_t pi2 = kir * i2_dimension + lir;
    double tl = ic_vals[pic2];
    for (int32_t iir = 0; iir < ir3_dimension; iir++) {
      int32_t pir3 = pir2 * ir3_dimension + iir;
      int32_t pi3 = pi2 * i3_dimension + iir;
      i_vals[pi3] = ir_vals[pir3] - tl;
    }
  }
}
// Generated by the Tensor Algebra Compiler (tensor-compiler.org)
for (int32_t pq = 0; pq < q1_dimension * q2_dimension; pq++) {
  q_vals[pq] = 0;
}
for (int32_t ki = 0; ki < i1_dimension; ki++) {
  for (int32_t li = 0; li < i2_dimension; li++) {
    int32_t pi2 = ki * i2_dimension + li;
    int32_t pi20 = ki * i2_dimension + li;
    for (int32_t ii = 0; ii < i3_dimension; ii++) {
      int32_t pi3 = pi2 * i3_dimension + ii;
      double ti = i_vals[pi3];
      for (int32_t ji = 0; ji < i3_dimension; ji++) {
        int32_t pi30 = pi20 * i3_dimension + ji;
        int32_t pq2 = ii * q2_dimension + ji;
        q_vals[pq2] = q_vals[pq2] + ti * i_vals[pi30];
      }
    }
  }
}
delete [] ir_vals;
delete [] ic_vals;
delete [] i_vals;
}
