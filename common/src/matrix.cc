#include "matrix.h"

// These are the 6 template arguments to every Eigen matrix
Eigen::Matrix<int, Eigen::Dynamic, 1> sortMatrix( Eigen::Matrix<int, Eigen::Dynamic, 1> target ) 
{
    // Manually construct a vector of correctly-typed matrix rows
    std::vector <int> matrixRows;
    for (unsigned int i = 0; i < target.rows(); i++) matrixRows.push_back( target(i,0) );
    std::sort(matrixRows.begin(), matrixRows.end() );

    Eigen::Matrix<int, Eigen::Dynamic, 1> sorted(matrixRows.size(),1);
    for (unsigned int i = 0; i < matrixRows.size(); i++)
            sorted(i,0) = matrixRows[i];
    return sorted;
}

void
drop_rows(Eigen::Matrix<float, Eigen::Dynamic, 1> &matrix, const Eigen::Matrix<int, Eigen::Dynamic, 1> &rows_to_remove)
{
	Eigen::Matrix<int, Eigen::Dynamic, 1> sorted_rows_to_remove=sortMatrix(rows_to_remove);

	unsigned int m = matrix.rows();
	unsigned int n = matrix.cols();

	for (size_t z=0; z<sorted_rows_to_remove.rows(); z++)
	{
		int dm=sorted_rows_to_remove(z,0)-z;
		int block_size;
		if (z<sorted_rows_to_remove.rows()-1) 
		{
			block_size=sorted_rows_to_remove(z+1,0)-sorted_rows_to_remove(z,0);
		} else {
			block_size=m-sorted_rows_to_remove(z,0)-1;
		}
		matrix.block(dm,0,block_size,n).noalias() = matrix.block(dm+z+1,0,block_size,n);
	}
	matrix.conservativeResize(m-sorted_rows_to_remove.rows(),n);
}

void 
drop_rows(Eigen::MatrixXf& matrix, const Eigen::Matrix<int, Eigen::Dynamic, 1> &rows_to_remove)
{
	Eigen::Matrix<int, Eigen::Dynamic, 1> sorted_rows_to_remove=sortMatrix(rows_to_remove);

	unsigned int m = matrix.rows();
	unsigned int n = matrix.cols();

	for (size_t z=0; z<sorted_rows_to_remove.rows(); z++)
	{
		int dm=sorted_rows_to_remove(z,0)-z;
		int block_size;
		if (z<sorted_rows_to_remove.rows()-1) 
		{
			block_size=sorted_rows_to_remove(z+1,0)-sorted_rows_to_remove(z,0);
		} else {
			block_size=m-sorted_rows_to_remove(z,0)-1;
		}
		matrix.block(dm,0,block_size,n) = matrix.block(dm+z+1,0,block_size,n);
	}
	matrix.conservativeResize(m-sorted_rows_to_remove.rows(),n);
}

void 
drop_cols(Eigen::MatrixXf& matrix, const Eigen::Matrix<int, Eigen::Dynamic, 1> &cols_to_remove)
{
	Eigen::Matrix<int, Eigen::Dynamic, 1> sorted_cols_to_remove=sortMatrix(cols_to_remove);

	unsigned int n = matrix.rows();
	unsigned int m = matrix.cols();

	for (size_t z=0; z<sorted_cols_to_remove.rows(); z++)
	{
		int dm=cols_to_remove(z,0)-z;
		int block_size;
		if (z<sorted_cols_to_remove.rows()-1) 
		{
			block_size=sorted_cols_to_remove(z+1,0)-sorted_cols_to_remove(z,0);
		} else {
			block_size=m-sorted_cols_to_remove(z,0)-1;
		}
		matrix.block(0,dm,n,block_size) = matrix.block(0,dm+z+1,n,block_size);
	}
	matrix.conservativeResize(n,m-sorted_cols_to_remove.rows() );
} 

void
drop_rows(Eigen::Matrix<double, Eigen::Dynamic, 1> &matrix, const Eigen::Matrix<int, Eigen::Dynamic, 1> &rows_to_remove)
{
	Eigen::Matrix<int, Eigen::Dynamic, 1> sorted_rows_to_remove=sortMatrix(rows_to_remove);

	unsigned int m = matrix.rows();
	unsigned int n = matrix.cols();

	for (size_t z=0; z<sorted_rows_to_remove.rows(); z++)
	{
		int dm=sorted_rows_to_remove(z,0)-z;
		int block_size;
		if (z<sorted_rows_to_remove.rows()-1) 
		{
			block_size=sorted_rows_to_remove(z+1,0)-sorted_rows_to_remove(z,0);
		} else {
			block_size=m-sorted_rows_to_remove(z,0)-1;
		}
		matrix.block(dm,0,block_size,n).noalias() = matrix.block(dm+z+1,0,block_size,n);
	}
	matrix.conservativeResize(m-sorted_rows_to_remove.rows(),n);
}

void 
drop_rows(Eigen::MatrixXd& matrix, const Eigen::Matrix<int, Eigen::Dynamic, 1> &rows_to_remove)
{
	Eigen::Matrix<int, Eigen::Dynamic, 1> sorted_rows_to_remove=sortMatrix(rows_to_remove);

	unsigned int m = matrix.rows();
	unsigned int n = matrix.cols();

	for (size_t z=0; z<sorted_rows_to_remove.rows(); z++)
	{
		int dm=sorted_rows_to_remove(z,0)-z;
		int block_size;
		if (z<sorted_rows_to_remove.rows()-1) 
		{
			block_size=sorted_rows_to_remove(z+1,0)-sorted_rows_to_remove(z,0);
		} else {
			block_size=m-sorted_rows_to_remove(z,0)-1;
		}
		matrix.block(dm,0,block_size,n) = matrix.block(dm+z+1,0,block_size,n);
	}
	matrix.conservativeResize(m-sorted_rows_to_remove.rows(),n);
}

void 
drop_cols(Eigen::MatrixXd& matrix, const Eigen::Matrix<int, Eigen::Dynamic, 1> &cols_to_remove)
{
	Eigen::Matrix<int, Eigen::Dynamic, 1> sorted_cols_to_remove=sortMatrix(cols_to_remove);

	unsigned int n = matrix.rows();
	unsigned int m = matrix.cols();

	for (size_t z=0; z<sorted_cols_to_remove.rows(); z++)
	{
		int dm=sorted_cols_to_remove(z,0)-z;
		int block_size;
		if (z<sorted_cols_to_remove.rows()-1) 
		{
			block_size=sorted_cols_to_remove(z+1,0)-sorted_cols_to_remove(z,0);
		} else {
			block_size=m-sorted_cols_to_remove(z,0)-1;
		}
		matrix.block(0,dm,n,block_size) = matrix.block(0,dm+z+1,n,block_size);
	}
	matrix.conservativeResize(n,m-sorted_cols_to_remove.rows() );
} 
