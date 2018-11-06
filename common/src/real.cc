#include "real.h"

const std::string Real::file_name=".dat";//!< The destination table in the Db.
const std::string Real::table_name="DATA";//!< The destination table in the Db.
const bool Real::binary=true;

const Registration Real::registered=Registration(Real::table_name, Real::create);

Real::Real (void){}
	
using namespace std;

void Real::read(std::istream& in) {
	std::string line, temp;
	std::getline(in, line);
	std::stringstream line_stream(line);
	line_stream >> value;
}

void Real::write_binary (std::ostream& out) const
{
	out.write( (char *) &value, sizeof(float_t) );
}

void Real::read_binary (std::istream& in)
{
	in.read( (char *) &value, sizeof(float_t) );
}

void Real::write (std::ostream& out) const 
{
	out << value;
}

std::string Real::header(void) const {
	return "@SCFNAME    \tPOS\tVALUE\n"; 
}

size_t Real::size() const{
	return sizeof(float_t); 
}

/*
const std::string Real::sql_header(void) const {
        return "(ABS_POS int, VALUE REAL, PRIMARY KEY (ABS_POS) )";
}

const std::string Real::sql_column_names(void) const {
        return "(ABS_POS, REAL)";
}

const std::string Real::sql_values(void) const {
        char return_buffer[SQL_LINE_SIZE]={0};
#if FLT_EVAL_METHOD == 2
	snprintf(return_buffer, SQL_LINE_SIZE, "(%d, %Lf)",
	abs_pos_
	value);
#elif FLT_EVAL_METHOD == 1
	snprintf(return_buffer, SQL_LINE_SIZE, "(%d, %f)",
	abs_pos_,
	value);
#elif FLT_EVAL_METHOD == 0
	snprintf(return_buffer, SQL_LINE_SIZE, "(%d, %f)",
	abs_pos_,
	value);
#endif
        return std::string(return_buffer);
}
*/

const std::string Real::get_file_name(void) const
{
	return file_name;
}

const std::string Real::get_table_name(void) const
{
	return table_name;
}

const  bool Real::get_binary(void) const
{
	return binary;
}
