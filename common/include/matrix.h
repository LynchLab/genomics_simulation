#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "Eigen/Core"
#include <iostream>

#define ARRAY2 	Eigen::ArrayXXd
#define ARRAY1 	Eigen::ArrayXd
#define MATRIX 	Eigen::MatrixXd
#define VECTOR 	Eigen::VectorXd
#define TYPE 	double

void drop_rows(Eigen::MatrixXf &, const Eigen::Matrix<int, Eigen::Dynamic, 1> &);
void drop_cols(Eigen::MatrixXf &, const Eigen::Matrix<int, Eigen::Dynamic, 1> &);

void drop_rows(Eigen::MatrixXd &, const Eigen::Matrix<int, Eigen::Dynamic, 1> &);
void drop_cols(Eigen::MatrixXd &, const Eigen::Matrix<int, Eigen::Dynamic, 1> &);

void drop_rows(Eigen::Matrix<float, Eigen::Dynamic, 1> &, const Eigen::Matrix<int, Eigen::Dynamic, 1> &);
void drop_cols(Eigen::Matrix<float, Eigen::Dynamic, 1> &, const Eigen::Matrix<int, Eigen::Dynamic, 1> &);

void drop_rows(Eigen::Matrix<double, Eigen::Dynamic, 1> &, const Eigen::Matrix<int, Eigen::Dynamic, 1> &);
void drop_cols(Eigen::Matrix<double, Eigen::Dynamic, 1> &, const Eigen::Matrix<int, Eigen::Dynamic, 1> &);

#endif 
