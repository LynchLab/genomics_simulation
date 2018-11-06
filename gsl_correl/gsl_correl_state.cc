#include "interface.h"
#include "map_file.h"
#include "state.h"
#include "sample_name.h"
#include "correl_data.h"
#include "real.h"

#include <cstring>
#include <sstream>
#include <tuple>
#include <map>
#include <fstream>
#include <ctime>
#include <omp.h>
#include <bitset>

#include "Eigen/Core"

#define UINT uint32_t
#define WORD 32
#define FOUR 2

#define t(a)	a.transpose()
#define Av(a,b) Eigen::VectorXd(a.array() * b.array() )
#define Am(a,b) Eigen::MatrixXd( a.array() * b.array() )
#define center(a,b) Eigen::MatrixXd( a.colwise() - b )

std::vector<size_t> sort_indexes(const UINT *v, const size_t &size) 
{
	std::vector<size_t> idx(size);
	iota(idx.begin(), idx.end(), 0);

	std::sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

	return idx;
}

static UINT Mask[WORD]={0x00000001, 0x00000002, 0x00000004, 0x00000008,
			0x00000010, 0x00000020, 0x00000040, 0x00000080,
			0x00000100, 0x00000200, 0x00000400, 0x00000800,
			0x00001000, 0x00002000, 0x00004000, 0x00008000,
			0x00010000, 0x00020000, 0x00040000, 0x00080000,
			0x00100000, 0x00200000, 0x00400000, 0x00800000,
			0x01000000, 0x02000000, 0x04000000, 0x08000000,
			0x10000000, 0x20000000, 0x40000000, 0x80000000};

inline void check_poly(const UINT *P, const UINT &N, bool *b)
{
	UINT type=P[0];
	UINT M = (1UL << WORD)-1;
	for (size_t x=0; x<N; x++)
	{
		if( (type & M) != (P[x] & M) )
		{
			M &= ~( (type & M) ^ (P[x] & M) ) ;
		}
		if (M == 0 )
			break;
	}
	for (size_t x=0; x<WORD; x++)
		b[x] = ~(M >> x) & 1U;
}

inline size_t c2bin(const double &c)
{
	return int(c*63);
}


inline void get_lambda (const UINT *l1, const UINT *l2, const UINT *m1, const UINT *m2, double *c, const UINT &N, Eigen::MatrixXd &l, Eigen::MatrixXd &lN)
{

	double *d_arr=new double [WORD];
	double **id_arr=new double *[N*2];

	for (size_t x=0; x<N*2; x++)
	{
		id_arr[x]=new double [WORD];
	}

	double **this_id=id_arr;
	double **that_id=this_id+N;
	double *p1 = new double [WORD];
	double *p2 = new double [WORD];

	//std::cerr << "getting double bit sum\n";

	double_bit_sumN(l1, l2, m1, m2, 2*N, p1, p2, d_arr, id_arr, c);

	//std::cerr << "got double bit sum\n";
	//bool test=false;

	for (size_t x=0; x<N; x++)
	{
		for (size_t y=0; y < WORD; y++)
		{
			if ( !isnan(d_arr[y]) && !isnan(p1[y]) && !isnan(p2[y]) )
			{
				//if ( p1[y] > 0.2  && p2[y] > 0.2 && p1[y] < 0.8 && p2[y] < 0.8 ) 
				if ( p1[y] > 0.001  && p2[y] > 0.001 && p1[y] < 0.999 && p2[y] < 0.999) 
				{
					//Variance weighting...
					l(x, c2bin(c[y]) ) += pow( d_arr[y]-p1[y]*p2[y], 1);//*(p1[y]*p2[y])*(1.-p1[y]*p2[y]);//( (*this_id)[y]-p1[y]*p2[y])/(c[y]*d_arr[y]) -1./c[y]+1.;
					//l(x, c2bin(c[y]) ) += pow( (*that_id)[y]-p1[y]*p2[y], 1)*(p1[y]*p2[y])*(1.-p1[y]*p2[y]);//( (*this_id)[y]-p1[y]*p2[y])/(c[y]*d_arr[y]) -1./c[y]+1.;
					//l(x, c2bin(c[y]) ) += ( (*this_id)[y]-d_arr[y]+c[y]*d_arr[y])/(c[y]*d_arr[y] )*( d_arr[y])*(1.-d_arr[y]);//( (*this_id)[y]-p1[y]*p2[y])/(c[y]*d_arr[y]) -1./c[y]+1.;

					lN(x, c2bin(c[y]) ) += 1;//*(p1[y]*p2[y])*(1.-p1[y]*p2[y]);//  (d_arr[y])*(1.-d_arr[y]);
				}
			}
		}
		this_id++;
		that_id++;
	}

	delete [] d_arr;
	for (size_t x=0; x<N*2; x++)
	{
		delete [] id_arr[x];
	}
	delete [] id_arr;

	delete [] p1;
	delete [] p2;
}

inline void get_f(const UINT *l1, const UINT *m1, const UINT &N, Eigen::MatrixXd &f, Eigen::MatrixXd &fN)
{
	//std::cerr << "HI!\n";
	double *p=new double [WORD];

	const UINT *this_place=l1;
	const UINT *that_place=l1+N;

	bit_sumN(l1, m1, 2*N, p);

	for (size_t x=0; x<N; x++)
	{
		for (size_t y=0; y < WORD; y++)
		{
			if ( ( (*m1) & Mask[y] ) != 0 && p[y] > 0.001  && p[y] < 0.999) 
			//if ( p[y] > 0.001  && p[y] < 0.999 ) 
			{
				 f(x,0) += 1.-double( int( (*this_place & Mask[y]) != (*that_place & Mask[y]) ) )/double(2.*p[y]*(1-p[y]) );
				fN(x,0) += 1.;
			}
		}
		this_place++;
		that_place++;
	}

	delete [] p;
}

int main (int argc, char **argv){

	size_t block=10000;

//#define block 20000

	std::cerr << __FILE__ << std::endl;

	std::string names_file="", input_file="", map_file="";
	int indX=-1, indY=-1, REP=8;
	double map_size=36;

	Environment env;

	double *r=new double [WORD];

	env.set_name("call_relatedness");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Calculates Lambda. You don't know what that is. Please direct questions to matthew.s.ackerman@gmail.com");

	env.optional_arg('i',"input",   input_file,      "please provide a number.", "state file.");
	env.optional_arg('c',"chiasma", map_size,      "please provide a number.", "total map size of the genome.");
	env.optional_arg('g',"genmap",  map_file,      "please provide a number.", "list of genetic map positions of SNPs.");
	env.optional_arg('r',"rep",  	REP,      "please provide a number.", "Rep till done.");
	env.positional_arg('n',"names", names_file,      "please provide a number.", "names of individuals.");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);

	// random pair (X,Y) of integers on (0,L), (X,L)
	// sort by X
	//Eigen::initParallel();

	Eigen::setNbThreads(4);
	std::cerr << "using " << Eigen::nbThreads( ) << " threads.\n";

	Indexed_file <Real> genetic_map;

	if (map_file!="")
	{
		genetic_map.open(map_file.c_str(), READ);
	}	

	State Pstates;
 	{
		Flat_file <State> state_file;
		if (input_file=="")
		{
			state_file.open(READ);
		}
		else
		{
			state_file.open(input_file.c_str(), READ);
		}
		Pstates=state_file.read_header();
		state_file.read(Pstates);
		state_file.close();
		std::cerr << "Done reading. Sample: " << Pstates.sample_size() << ", " << " genome: "<<  Pstates.genome_size() << ", masked: "<< Pstates.masked() << std::endl;
	}

	int N(Pstates.sample_size());
	int J=32;
	int LEN(Pstates.genome_size()*32);

	//dist, 11, N, individual
	count[DISTANCE][2][N];

	UINT p_arr[WORD];
	UINT h_arr[WORD];

	std::default_random_engine rng;

	size_t ll=0;
	size_t K1k=0;
	size_t K2k=0;

	UINT *P=new uint32_t [N*2];
	UINT *P2=P+N;
	UINT *P1=P;

	UINT *M=new uint32_t [N*2];
	UINT *M2=M+N;
	UINT *M1=M;

	Eigen::MatrixXd lambda=Eigen::MatrixXd::Zero(N, J);
	Eigen::MatrixXd lambdaN=Eigen::MatrixXd::Zero(N, J);

	Eigen::MatrixXd f=Eigen::MatrixXd::Zero(N,1);
	Eigen::MatrixXd fN=Eigen::MatrixXd::Zero(N,1);

	clock_t t1, t2, t3, s;
	size_t SK1=0, SK2=0, K1=0, K2=0;

	bool b_arr[WORD];

	ll=0;

	Pstates.rewind();
	std::vector <size_t> poly;
	std::vector <double> c;
	double p[32];

	std::cerr << "Scan 1\n";
	Real cm;
	if  (map_file!="")
		genetic_map.read_header();

	while (!Pstates.empty())
	{
		Pstates.uncompress(P1, P2, M1, M2);

		bit_sumN(P1, M1, N*2, p);
		
		for (size_t x=0; x < WORD; x++)
		{
			if  (map_file!="")
				genetic_map.read(cm);
			else
				cm.value=FOUR*ll*double(map_size)/(double(LEN) )*100;
			//if ( ( p[x] > 0.005 and p[x] < 0.995 ) && not(p[x]>0.45 && p[x]<0.55 )  ) 
			if ( ( p[x] > 0.001 and p[x] < 0.999 )   ) 
			{
				poly.push_back(ll);
				c.push_back(cm.value);
				std::cerr << cm.value << std::endl;
			}
			ll++;
		}
	}
	genetic_map.close();

	std::cerr << "POLY SIZE: " << poly.size() << "\n";

	for (size_t i=0; i < REP; i++) 
	{
		Pstates.rewind();
		ll=0;

		while (!Pstates.empty())
		{
			Pstates1.uncompress(P1_1, P1_2, M1_1, M1_2);
			Pstates2.uncompress(P2_1, P2_2, M2_1, M2_2);
			get_correl(counts, P1_1, P2_1, M1_1, M2_1, c);
		}
	}

	std::cout << "NaN" << '\t' << "NaN";
	for ( size_t y=0; y<J; y++)
	{
		std::cout << '\t' << lambdaN(0,y);// << "(" << lambdaN(x,y) << ")";
	}
	std::cout << std::endl;

	for ( size_t x=0; x<N; x++)
	{
		std::cout << x << '\t' << f(x,0)/fN(x,0);
	//	std::cout << x << '\t' << f(x,0) << ", " << fN(x,0);
		for ( size_t y=0; y<J; y++)
		{
			std::cout << '\t' << lambda(x,y)/lambdaN(x,y);// << "(" << lambdaN(x,y) << ")";
		//	std::cout << '\t' << lambdaN(x,y);
		}
		std::cout << std::endl;
	}
}
