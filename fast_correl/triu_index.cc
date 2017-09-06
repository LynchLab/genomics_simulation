#include "triu_index.h"

Triangular_index::Triangular_index(size_t size)
{
	n_=size;
	size_= (n_*(n_+1)/2);
	k_=0;
	x_=0;
	y_=0;
	dirty_=false;
}


void 
Triangular_index::clean_(void)
{
	if (k_<=size_){
		int k =n_*(n_+1)/2-1-k_;
  		int K = floor( (sqrt(8*k+1)-1)/2);
		x_ = n_-1-K;
		y_ = k_ -n_*(n_+1)/2 + (K+1)*(K+2)/2+n_-1-K;
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Array index out of bounds.\n"), __FILE__, __LINE__);
	}
	dirty_=false;
}

void 
Triangular_index::get_xy(size_t &x, size_t &y) 
{
	if (dirty_) clean_();
	x=x_;
	y=y_;
}

size_t 
Triangular_index::get_k(const size_t &x, const size_t &y) 
{
	set(x,y);
	return k_;
}

size_t 
Triangular_index::size(void) const
{
	return size_;
}

void 
Triangular_index::set(const size_t &x, const size_t &y)
{
	x_=x;
	y_=y;
	k_ = (n_*(n_+1)/2) - (n_-x)*((n_-x)+1)/2 + y - x;
	dirty_=false;
}

Triangular_index& Triangular_index::operator++()
{
	k_++;
	if (++y_ >= n_)
	{
		y_=(++x_);
	}
	return *this;
}

Triangular_index Triangular_index::operator++(int)
{
	Triangular_index tmp(*this); // copy
	operator++(); // pre-increment
	return tmp;   // return old value
}
