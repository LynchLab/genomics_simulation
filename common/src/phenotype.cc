#include "phenotype.h"

const std::string Phenotype::file_name=".phe";
const std::string Phenotype::table_name="PHENOTYPE";
const bool Phenotype::binary=false;

const Registration Phenotype::registered=Registration(Phenotype::table_name, Phenotype::create);

Phenotype::Phenotype ()
{
	delim='\t';
}

Phenotype::Phenotype (const std::vector<std::string> &fields)
{
	delim='\t';
	_n_traits=fields.size()-1;
	_n_samples=0;
	trait=std::vector <std::string> (fields.begin()+1, fields.end() );
	value=Eigen::MatrixXd::Zero(0,_n_traits);
}

Phenotype::Phenotype (const size_t &N)
{
	delim='\t';
	_n_traits=N;
	_n_samples=0;
	value=Eigen::MatrixXd::Zero(0,_n_traits);
}

void
Phenotype::read (std::istream& in) 
{
	std::string line;
	std::getline(in, line);
	std::vector <std::string> fields=split(line, delim);

	sample_name.push_back(fields[0]);
	value.conservativeResize(++_n_samples, _n_traits);

	for (size_t y=0; y<_n_traits; y++)
	{
		value(_n_samples-1,y)=atof(fields[1+y].c_str());
	}
}

void
Phenotype::write (std::ostream& out) const 
{
	for (size_t x=0; x<_n_samples; x++)
	{
		std::cout << sample_name[x] << delim;
		for (size_t y=0; y<_n_traits; y++)
		{
			std::cout << value(x,y) << delim;
		}
		if (x!=_n_samples-1) std::cout << std::endl;
	}
}

std::string 
Phenotype::header(void) const 
{
	std::stringstream ret; 
	ret << "@SAMPLE";
	for (size_t x=0; x<_n_traits; x++)
	{
		 ret << delim << trait[x];
	}
	ret << std::endl;
	return ret.str();
}

size_t 
Phenotype::size(void) const 
{
	//TODO IMPLEMENT
	return 0;
}

void 
Phenotype::clear(void)
{
	zero();
}

void 
Phenotype::zero(void)
{
	//TODO IMPLEMENT
}

const size_t
Phenotype::sample_size(void) const
{
	return _n_samples;
}

const bool 
Phenotype::get_binary(void) const
{
	return binary;
}

const std::string Phenotype::get_file_name(void) const
{
	return file_name;
}

const std::string Phenotype::get_table_name(void) const
{
	return table_name;
}
	

Eigen::Matrix<int, Eigen::Dynamic, 1> 
Phenotype::get_missing(const int &trait_num) const
{
//	std::cerr << "hi!\n";
	std::vector<int> ret_pre;
	Eigen::Matrix<int, Eigen::Dynamic, 1> ret;
	for (int z=0; z<value.rows(); z++)
	{
		if (value(z, trait_num)==0) ret_pre.push_back(z);
	}
	ret=Eigen::Matrix<int, Eigen::Dynamic, 1>::Zero(ret_pre.size(),1);
	for (int z=0; z<ret_pre.size(); z++)
	{
		ret(z,0)=ret_pre[z];
	}
	return ret;
}

void
Phenotype::drop_missing(Eigen::Matrix<int, Eigen::Dynamic, 1> &drops) 
{
	drop_rows(value, drops);
}
