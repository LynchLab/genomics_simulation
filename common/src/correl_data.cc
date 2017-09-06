#include "correl_data.h"

const std::string Correl::file_name=".cor";
const std::string Correl::table_name="CORREL";
const bool Correl::binary=false;

const Registration Correl::registered=Registration(Correl::table_name, Correl::create);

Correl::Correl (){
	delim='\t';
}

Correl::Correl (const std::string &X, const std::string &Y){
	delim='\t';
}

void 
Correl::set_X_name(const std::string &X)
{
//	X_=X;
}

void
Correl::set_X_name(const id0_t &X)
{
	X_=X;
}

void
Correl::set_Y_name(const id0_t &Y)
{
	Y_=Y;
}

void 
Correl::set_Y_name(const std::string &Y)
{
//	Y_=Y;
}

void
Correl::read (std::istream& in) 
{
	std::string line;
	std::getline(in, line);
	std::stringstream stream_line(line);	
	stream_line >> X_ ;
	stream_line >> Y_ ;
	stream_line >> b_ ;
	stream_line >> d_ ;
	stream_line >> bd_ ;
}

void
Correl::write (std::ostream& out) const 
{
	out << X_ << delim;
	out << Y_ << delim;
	out << b_ << delim;
	out << d_ << delim;
	out << bd_;
}

std::string 
Correl::header(void) const {
	return std::string("@SAMPLE_X\tSAMPLE_Y\tBETA_XY\tDELTA_XY\tGAMMA_XY\n");
}

size_t 
Correl::size(void) const {
	//18
	return sizeof(id1_t)+sizeof(float_t)*E_LIM*2+sizeof(float_t)*18+sizeof(char)+sizeof(id0_t)*2;
}

Correl& 
Correl::operator=(const Correl &rhs)
{
	//memcpy(e_X_, rhs.e_X_, E_LIM*sizeof(float_t) );
	//memcpy(e_Y_, rhs.e_Y_, E_LIM*sizeof(float_t) );
	b_=rhs.b_;
	d_=rhs.d_;
	bd_=rhs.bd_;
	return *this;
}

const bool 
Correl::get_binary(void) const
{
	return binary;
}

const std::string Correl::get_file_name(void) const
{
	return file_name;
}

const std::string Correl::get_table_name(void) const
{
	return table_name;
}
