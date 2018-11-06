#ifndef _RELATEDNESS_DATA_H_
#define _RELATEDNESS_DATA_H_

#include <cstring>

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>
#include <sstream>

#include "Eigen/Core"

#include "typedef.h"
#include "data.h"
#include "stream_tools.h"
#include "matrix.h"

#define E_LIM 25

/// Relatedness data.
class Relatedness : public Data{ 
private:
	void write (std::ostream&) const;	//!< use to write Allele. Inherits <<
	void read (std::istream&);		//!< use to read Allele. Inherits >>
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Relatedness(Columns);
	}
	bool sub_set_;
	bool tr_;
	size_t a_, b_;
	size_t n_;
	
	real_t mean_sites_analyzed_;
	real_t sd_sites_analyzed_;


public:
	char delim;	//!< the delimiter used when reading/writing the class in text mode.	

//	std::string X_;	//!< the name of the first (X) sample in the comparison.
//	std::string Y_;	//!< the name of the second (Y) sample in the comparison.

	id0_t _X, _Y;
	size_t _N;

	std::vector <std::string> sample_name;
	id1_t sites;	//!< the number of sites analyzed.

	VECTOR Mt1;		//!< transpose(M) * 1
	VECTOR Ht1;

	VECTOR rMt1;
	VECTOR rHt1;
	VECTOR lMt1;
	VECTOR lHt1;

	MATRIX MtM;		//!< transpose(M)*M
	MATRIX MtH;
	MATRIX HtH;

	MATRIX MJM;		//!< M*J*M
	MATRIX MJH;
	MATRIX HJH;		

	MATRIX tr;		//!< traces of MtM, HtH, etc..

	//Higher order interactions...

	//TODO Vector of descriptors ..?? 
	/*
	VECTOR Xt1;
	MATRIX XtX(N,N);
	*/

	Relatedness();	
	Relatedness(const std::vector <std::string> &); 
	Relatedness(const size_t &); 

	//! The header line of plain text files. 
	std::string header(void) const;
	//! Size in bytes for binary read/write. 
	size_t size(void) const;

	//! zeros statistics and sets names to empty.
	void clear(void); 
	//! zeros statistics, but doesn't set names to empty.
	void zero(void);  

	static const std::string file_name;	//!< The dafualt extention for files.
	static const std::string table_name;	//!< Destination table in Db.

	const std::string get_file_name(void) const;
	const std::string get_table_name(void) const;

	static const bool binary;

	const bool get_binary() const;

	const size_t sample_size(void) const;
	
	inline void set_mean_sites(const real_t &s){mean_sites_analyzed_=s;};
	inline void set_std_sites(const real_t &s){sd_sites_analyzed_=s;};

	inline const real_t& mean_sites(void) const {return mean_sites_analyzed_;};
	inline const real_t& std_sites(void) const {return sd_sites_analyzed_;};

	void make_outer_prod(void);
	void make_trace(void);

	Relatedness block(const int &, const int &, const int &, const int &) const;

	void drop_missing(const Eigen::Matrix<int, Eigen::Dynamic, 1> &);

	Relatedness& operator=(const Relatedness &rhs);
};

#endif
