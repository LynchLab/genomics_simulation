#ifndef _PHENOTYPE_H_
#define _PHENOTYPE_H_

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

/// Phenotype data.
class Phenotype : public Data{ 
private:
	void write (std::ostream&) const;	//!< use to write Allele. Inherits <<
	void read (std::istream&);		//!< use to read Allele. Inherits >>
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Phenotype(Columns);
	}
	size_t _n_samples, _n_traits;
public:
	char delim;	//!< the delimiter used when reading/writing the class in text mode.	

	std::vector <std::string> sample_name;
	std::vector <std::string> trait;

	Eigen::MatrixXd value;

	Phenotype();	
	Phenotype(const std::vector <std::string> &); 
	Phenotype(const size_t &); 

	//! The header line of plain text files. 
	std::string header(void) const;
	//! Size in bytes for binary read/write. 
	size_t size(void) const;

	//! zeros values and sets names to empty.
	void clear(void); 
	//! zeros values, but doesn't set names to empty.
	void zero(void);  

	static const std::string file_name;	//!< The dafualt extention for files.
	static const std::string table_name;	//!< Destination table in Db.

	const std::string get_file_name(void) const;
	const std::string get_table_name(void) const;

	Eigen::Matrix<int, Eigen::Dynamic, 1> get_missing(const int &) const;
	void drop_missing(Eigen::Matrix<int, Eigen::Dynamic, 1> &);

	static const bool binary;

	const bool get_binary() const;

	const size_t sample_size(void) const;

};

#endif
