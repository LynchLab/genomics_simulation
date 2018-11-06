#include "relatedness_data.h"

const std::string Relatedness::file_name=".rel";
const std::string Relatedness::table_name="RELATEDNESS";
const bool Relatedness::binary=false;

const Registration Relatedness::registered=Registration(Relatedness::table_name, Relatedness::create);

Relatedness::Relatedness ()
{
	mean_sites_analyzed_=0;
	sd_sites_analyzed_=0;
	delim='\t';
	tr_=false;
}

Relatedness::Relatedness (const std::vector<std::string> &fields)
{
	std::cerr << "Initalizing...\n";
	delim='\t';
	_N=(fields.size()-3)/3;
	std::cerr << _N << std::endl;

	for (size_t x=0; x<_N; x++)
	{
		sample_name.push_back("sample_"+std::to_string(x) );
	}

	std::vector <std::string> three=split(fields[0], '|');
	std::cerr << three[1] << ", " << three[2] << std::endl;
	mean_sites_analyzed_=atof(three[1].c_str() );
	sd_sites_analyzed_=atof(three[2].c_str() );

	Mt1=VECTOR::Zero(_N);
	Ht1=VECTOR::Zero(_N);
	MtM=MATRIX::Zero(_N,_N);
	MtH=MATRIX::Zero(_N,_N);
	HtH=MATRIX::Zero(_N,_N);	

	tr_=false;
}

Relatedness::Relatedness (const size_t &N)
{
	mean_sites_analyzed_=0;
	sd_sites_analyzed_=0;
	delim='\t';
	_N=N;
	for (size_t x=0; x<_N; x++)
	{
		sample_name.push_back("sample_"+std::to_string(x) );
	}

	Mt1=VECTOR::Zero(_N);
	Ht1=VECTOR::Zero(_N);
	MtM=MATRIX::Zero(_N,_N);
	MtH=MATRIX::Zero(_N,_N);
	HtH=MATRIX::Zero(_N,_N);	

	tr_=false;}

void
Relatedness::read (std::istream& in) 
{
	for (size_t x=0; x<_N; x++)
	{
		std::string line;
		std::getline(in, line);
		std::vector <std::string> fields=split(line, delim);

		sample_name[x]=fields[0];

		Mt1(x)=atof(fields[1].c_str() );
		Ht1(x)=atof(fields[2].c_str() );

		for (size_t y=0; y<_N; y++)
		{
			MtH(x,y)=atof(fields[3+y].c_str());
		}
		for (size_t y=0; y<_N; y++)
		{
			MtM(x,y)=atof(fields[3+_N+y].c_str() );
		}
		for (size_t y=0; y<_N-1; y++)
		{
			HtH(x,y)=atof(fields[3+2*_N+y].c_str() );
		}
		HtH(x,_N-1)=atof(fields[3+2*_N+_N-1].c_str() );
	}
}

void
Relatedness::write (std::ostream& out) const 
{
	for (size_t x=0; x<_N; x++)
	{
		out << sample_name[x] << delim;
		out.precision(17);
		out << Mt1(x) << delim;
		out << Ht1(x) << delim;
		for (size_t y=0; y<_N; y++)
		{
			out << MtH(x,y) << delim;
		}
		for (size_t y=0; y<_N; y++)
		{
			out << MtM(x,y) << delim;
		}
		for (size_t y=0; y<_N-1; y++)
		{
			out << HtH(x,y) << delim;
		}
		out << HtH(x,_N-1);
		if (x!=_N-1) out << std::endl;
	}
}

std::string 
Relatedness::header(void) const 
{
	std::stringstream ret; 
	ret << "@SAMPLE|" << mean_sites_analyzed_ << "|" << sd_sites_analyzed_ << delim << "Mt1:" << delim << "Ht1:";
	for (size_t x=0; x<_N; x++) 
	{
		 ret << delim << "MtH_" << x;
	}
	for (size_t x=0; x<_N; x++) 
	{
		 ret << delim << "MtM_" << x;
	}
	for (size_t x=0; x<_N; x++) 
	{
		 ret << delim << "HtH_" << x;
	}
	ret << std::endl;
	return ret.str();
}

size_t 
Relatedness::size(void) const 
{
	//TODO IMPLEMENT
	return 0;
}

void 
Relatedness::clear(void)
{
	zero();
}

void 
Relatedness::zero(void)
{
	//TODO IMPLEMENT
}

Relatedness& 
Relatedness::operator=(const Relatedness &rhs)
{
	sample_name=rhs.sample_name;

	mean_sites_analyzed_=rhs.mean_sites_analyzed_;
	sd_sites_analyzed_=rhs.sd_sites_analyzed_;

	HtH=rhs.HtH;
	MtM=rhs.MtM;
	MtH=rhs.MtH;

	HJH=rhs.HJH;
	MJM=rhs.MJM;
	MJH=rhs.MJH;

	Mt1=rhs.Mt1;
	Ht1=rhs.Ht1;

	_N=rhs._N;

	tr=rhs.tr;
	tr_=rhs.tr_;
	//TODO IMPLEMENT
	return *this;
}

const size_t
Relatedness::sample_size(void) const
{
	return _N;
}

const bool 
Relatedness::get_binary(void) const
{
	return binary;
}

const std::string Relatedness::get_file_name(void) const
{
	return file_name;
}

const std::string Relatedness::get_table_name(void) const
{
	return table_name;
}

void
Relatedness::make_outer_prod(void)
{
	MJM=Mt1*Mt1.transpose();
	MJH=Mt1*Ht1.transpose()+Ht1*Mt1.transpose();
	HJH=Ht1*Ht1.transpose();
}

void
Relatedness::make_trace(void)
{
	if (not (tr_) )
	{
		tr=MATRIX::Zero(7,1);
		tr(0,0)=MtM.trace();
		tr(1,0)=(Mt1*Mt1.transpose() ).trace();
		tr(2,0)=HtH.trace();
		tr(3,0)=(Ht1*Ht1.transpose() ).trace();
		tr(4,0)=MtH.trace();
		tr(5,0)=(Mt1*Ht1.transpose()+Ht1*Mt1.transpose() ).trace();
		tr(6,0)=Mt1.rows();
		tr_=true;
	}
}


	
Relatedness Relatedness::block(const int &a, const int &b, const int &c, const int &d) const
{
	Relatedness ret;
	ret.sample_name=this->sample_name;
	ret._N=this->_N;

	ret.mean_sites_analyzed_=this->mean_sites_analyzed_;
	ret.sd_sites_analyzed_=this->sd_sites_analyzed_;

	ret.MtM=this->MtM.block(a,b,c,d);
	ret.MtH=this->MtH.block(a,b,c,d);
	ret.HtH=this->HtH.block(a,b,c,d);

	ret.MJM=this->MJM.block(a,b,c,d);
	ret.MJH=this->MJH.block(a,b,c,d);
	ret.HJH=this->HJH.block(a,b,c,d);

	ret.rMt1=this->Mt1.block(a,0,c,1);
	ret.lMt1=this->Mt1.block(b,0,d,1);
	ret.rHt1=this->Ht1.block(a,0,c,1);
	ret.lHt1=this->Ht1.block(b,0,d,1);

	return ret;
}
	
void 
Relatedness::drop_missing(const Eigen::Matrix<int, Eigen::Dynamic, 1> &drop)
{
	std::cerr << "Dropping Mt1...\n";
	drop_rows(Mt1, drop);
	std::cerr << "Dropping Ht1...\n";
	drop_rows(Ht1, drop);
	
	std::cerr << "Dropping MtM...\n";
	drop_rows(MtM, drop);
	drop_cols(MtM, drop);

	std::cerr << "Dropping HtH...\n";
	drop_rows(HtH, drop);
	drop_cols(HtH, drop);

	std::cerr << "Dropping MtH...\n";
	drop_rows(MtH, drop);
	drop_cols(MtH, drop);
}
