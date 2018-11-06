#ifndef _TYPES_
#define _TYPES_

#include "Eigen/Core"
#include "Eigen/Dense"

#define t(a)	a.transpose()
#define Av(a,b) VECTOR(a.array() * b.array() )
#define o(a,b) MATRIX( a.array() * b.array() )
//#define d(a,b) MATRIX( a.array() / b.array() )
#define center(a,b) MATRIX( a.colwise() - b )
#define center_scalar(a,b) MATRIX( a.array() - b )
//#define D(a,b) MATRIX( a.array().colwise() / b.array() )

#endif
