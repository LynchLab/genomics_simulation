#ifndef _REAL_DATA_H_
#define _REAL_DATA_H_

#include <cstring>

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>
#include <sstream>

#include "typedef.h"
#include "data.h"

/// Real data.
class Real : public Indexed_data{ 
private:
	void write (std::ostream&) const;	//!< use to write Allele. Inherits <<
	void read (std::istream&);		//!< use to read Allele. Inherits >>
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Real(Columns);
	}
public:
	void write_binary (std::ostream& out) const;
	void read_binary (std::istream& in);

	char delim;	//!< the delimiter used when reading/writing the class in text mode.	

	float_t value;

	Real();	
	//! delegating a neccisary constructor.	
	Real(const std::vector <std::string> &) : Real(){}; 
	//! construct with names. 
	Real(const std::string &, const std::string &);		  

	//! The header line of plain text files. 
	std::string header(void) const;
	//! Size in bytes for binary read/write. 
	size_t size(void) const;

	static const std::string file_name;	//!< The dafualt extention for files.
	static const std::string table_name;	//!< Destination table in Db.

	const std::string get_file_name(void) const;
	const std::string get_table_name(void) const;

	static const bool binary;

	const bool get_binary() const;

	Real& operator=(const Real &rhs);
};

#endif
