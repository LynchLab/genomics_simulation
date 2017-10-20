#ifndef _CORREL_DATA_H_
#define _CORREL_DATA_H_

#include <cstring>

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>
#include <sstream>

#include "typedef.h"
#include "data.h"

#define E_LIM 25

/// Correl data.
class Correl : public Data{ 
private:
	void write (std::ostream&) const;	//!< use to write Allele. Inherits <<
	void read (std::istream&);		//!< use to read Allele. Inherits >>
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Correl(Columns);
	}
public:
	char delim;	//!< the delimiter used when reading/writing the class in text mode.	

//	std::string X_;	//!< the name of the first (X) sample in the compairison.
//	std::string Y_;	//!< the name of the second (Y) sample in the compairison.

	id0_t X_, Y_;

	id1_t sites;	//!< the number of sites analyzed.

	float_t b_, d_, bd_;// f_, k_, d1_, d2_, fx_, fy_;
	/*
	 * See ?
	 */

	Correl();	
	//! delegating a neccisary constructor.	
	Correl(const std::vector <std::string> &) : Correl(){}; 
	//! construct with names. 
	Correl(const std::string &, const std::string &);		  

	//! The header line of plain text files. 
	std::string header(void) const;
	//! Size in bytes for binary read/write. 
	size_t size(void) const;

	//! Sets the name of the X sample. 
	void set_X_name(const std::string &);
	void set_X_name(const id0_t &);
	//! Sets the name of the Y sample. 
	void set_Y_name(const std::string &);
	void set_Y_name(const id0_t &);

	static const std::string file_name;	//!< The dafualt extention for files.
	static const std::string table_name;	//!< Destination table in Db.

	const std::string get_file_name(void) const;
	const std::string get_table_name(void) const;

	static const bool binary;

	const bool get_binary() const;

	Correl& operator=(const Correl &rhs);
};

#endif
