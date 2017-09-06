#ifndef TRIANGULAR_INDEX_H
#define TRIANGULAR_INDEX_H

#include <cstdlib>
#include <math.h>
#include <cstdio>
#include <libintl.h>


class Triangular_index 
{
private:
	size_t k_, n_, x_, y_, size_;
	bool dirty_;
	void clean_(void);
public:
	Triangular_index(size_t size);

	void get_xy(size_t &x, size_t &y);
	
	inline size_t get_k(void) const {return k_;};
	
	size_t get_k(const size_t &x, const size_t &y);
	size_t size(void) const;

	void set(const size_t &x, const size_t &y);


	Triangular_index& operator++();

	Triangular_index operator++(int);

	friend bool operator< (const Triangular_index& lhs, const Triangular_index& rhs);
	friend bool operator> (const Triangular_index& lhs, const Triangular_index& rhs);
	friend bool operator<=(const Triangular_index& lhs, const Triangular_index& rhs);
	friend bool operator>=(const Triangular_index& lhs, const Triangular_index& rhs);
	friend bool operator==(const Triangular_index& lhs, const Triangular_index& rhs);
	friend bool operator!=(const Triangular_index& lhs, const Triangular_index& rhs);

	Triangular_index operator+(const size_t &a) const
	{
		Triangular_index ret=*this;
		ret.k_+=a;
		ret.clean_();
		return ret;	
	}

	int operator-(const Triangular_index &t) const
	{
		return k_-t.k_;
	}


	Triangular_index& operator+=(const int &i){
		this->k_ += i;
		clean_();
		return *this;
	}

	/*Triangular_index operator+(const Triangular_index &t) const
	{
		Triangular_index ret=*this;
		ret.k_+=t.k_;
		ret.clean_();
		return ret;	
	}*/
};

inline bool operator< (const Triangular_index& lhs, const Triangular_index& rhs)
{
	return lhs.k_ < rhs.k_; 
}

inline bool operator> (const Triangular_index& lhs, const Triangular_index& rhs){ return rhs < lhs; }
inline bool operator<=(const Triangular_index& lhs, const Triangular_index& rhs){ return !(lhs > rhs); }
inline bool operator>=(const Triangular_index& lhs, const Triangular_index& rhs){ return !(lhs < rhs); }

inline bool operator==(const Triangular_index& lhs, const Triangular_index& rhs)
{
	return lhs.k_==rhs.k_;
}

inline bool operator!=(const Triangular_index& lhs, const Triangular_index& rhs){ return !(lhs == rhs); }
#endif
