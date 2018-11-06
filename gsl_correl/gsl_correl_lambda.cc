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


inline void double_bit_sumN(const UINT *l1, const UINT *l2, const UINT *m1, const UINT *m2, const UINT &N, double *p1, double *p2, double *D, double **Di, double *c)
{
	const UINT *end=l1+N;
	const UINT *this_place=l1;
	const UINT *that_place=l2;
	UINT *Tt=new UINT[WORD];
	double **this_id=Di;
	UINT A, B;

	memset(Tt, 0,  sizeof(UINT) * WORD);
	memset(p1, 0, sizeof(double) * WORD);
	memset(p2, 0, sizeof(double) * WORD);
	memset(D, 0, sizeof(double) * WORD);

//#define X0
#ifdef X0

	UINT x0[WORD];
	UINT x1[WORD];
	UINT x2[WORD];
	UINT x3[WORD];

	memset(x0, 0,  sizeof(UINT) * WORD);
	memset(x1, 0,  sizeof(UINT) * WORD);
	memset(x2, 0,  sizeof(UINT) * WORD);
	memset(x3, 0,  sizeof(UINT) * WORD);
#endif 

	while (this_place!=end) 
	{
		for (size_t x=0; x < WORD; x++)
		{
			//CHECK THING!
			if ( ( ( (*m1) & Mask[x] ) != 0 ) && ( (*m2) & Mask[x] != 0 ) )
			{
				A = ( (*this_place & Mask[x] ) !=0 );
				B = ( (*that_place & Mask[x] ) !=0 );
				(*this_id)[x] = ( ( A & B ) != 0 );
				p1[x] += A;
				p2[x] += B;
				Tt[x] += 1;
				D[x] += ( (A & B ) != 0 );
#ifdef X0
				x0[x] += A * B;
				x1[x] += A * (1-B);
				x2[x] += (1-A) * B;
				x3[x] += (1-A) * (1-B);
#endif 
			}
		}
		this_id++;
		this_place++;
		that_place++;
		m1++;
		m2++;
	}
			
	for (size_t x=0; x < WORD; x++)
	{

#ifdef X0
		if (x0[x]+x1[x] > 0 && x0[x]+x2[x] > 0 )
			std::cerr  << c[x] << ", " << x0[x] <<  ", " << x1[x] << ", " <<  x2[x] << ", " << x3[x] << std::endl;
#endif 

		p1[x] = p1[x]/double(Tt[x]);
		p2[x] = p2[x]/double(Tt[x]);
		D[x]  = D[x]/double(Tt[x]);//-p1[x]*p2[x];
	}

	delete [] Tt;
}

inline void bit_sumN(const UINT *site, const UINT *mask, const UINT &N, double *p)
{
	const UINT *this_place=site;
	const UINT *end=site+N;
	const UINT *m1=mask;
	UINT Tt[32];

	memset(Tt, 0, sizeof(Tt[0])*32);
	memset(p,  0, sizeof( p[0])*32);

	while (this_place!=end) 
	{
		for (size_t x=0; x< WORD; x++)
		{

			if ( ( (*m1) & Mask[x] ) != 0 )
			{
				p[x] += ( (*this_place & Mask[x] ) !=0 );
				Tt[x] += 1;
			}
		}
		this_place++;
		m1++;
	}

	for (size_t x=0; x < WORD; x++)
	{
		p[x]=p[x]/double(Tt[x]);
	}
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

	Eigen::VectorXd s1=Eigen::ArrayXd::Zero(N);
	Eigen::VectorXd s2=Eigen::ArrayXd::Zero(N);
	Eigen::VectorXd s3=Eigen::ArrayXd::Zero(N);

	Eigen::VectorXd fs(N);
	Eigen::VectorXd dm(N);

	Eigen::MatrixXd S4=Eigen::ArrayXXd::Zero(N, N);
	Eigen::MatrixXd S5=Eigen::ArrayXXd::Zero(N, N);
	Eigen::MatrixXd S6=Eigen::ArrayXXd::Zero(N, N);

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

	UINT *L11=new uint32_t [N*2*block];
	UINT *L12=L11+N;
	UINT *L1=L11;
	UINT *L21=new uint32_t [N*2*block];
	UINT *L22=L21+N;
	UINT *L2=L21;

	UINT *mL11=new uint32_t [N*2*block];
	UINT *mL12=mL11+N;
	UINT *mL21=new uint32_t [N*2*block];
	UINT *mL22=mL21+N;

	UINT *L1c=new uint32_t [32*block];
	UINT *L2c=new uint32_t [32*block];

	Eigen::MatrixXd lambda=Eigen::MatrixXd::Zero(N, J);
	Eigen::MatrixXd lambdaN=Eigen::MatrixXd::Zero(N, J);

	Eigen::MatrixXd f=Eigen::MatrixXd::Zero(N,1);
	Eigen::MatrixXd fN=Eigen::MatrixXd::Zero(N,1);

	clock_t t1, t2, t3, s;
	size_t SK1=0, SK2=0, K1=0, K2=0;

	bool b_arr[WORD];

	UINT *polyA = new UINT[block*WORD];
	UINT *polyB = new UINT[block*WORD];
	UINT *S1   =  new UINT[block*WORD];
	UINT *S2   =  new UINT[block*WORD];
	UINT *S1k  =  new UINT[block*WORD];
	UINT *S2k  =  new UINT[block*WORD];

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

	UINT Z=poly.size();

	std::uniform_int_distribution <size_t> unif_L(0, Z-1);
	std::exponential_distribution <double> geom_L( 1./100. );

	//std::uniform_int_distribution <size_t> unif_L2(0, 500);

	for (size_t x=0; x < block*WORD; x++)
	{	
		UINT A=unif_L(rng);
		//poly[unif_L(rng)];
		//UINT B=gsl_ran_geometric (const gsl_rng * r, double p);
		double B=c[A]+geom_L(rng);

		UINT C=A+1;
		while (c[C++] < B && C < Z-1);

		polyA[x]=A;
		polyB[x]=C;

		A=poly[A];
		B=poly[C];

		//if (A > B) std::swap(A,B);
		//if (A == 28046 || B == 28046 )  std::cerr << A << ", " << B << ", " << B-A << std::endl;
		S1[x]=A >> 5;
		S2[x]= int(B) >> 5;
		S1k[x]= A % 0x20;
		S2k[x]= int(B) % 0x20;
	}

	std::vector<size_t> SK1_2_x=sort_indexes(S1, block*WORD);
	std::vector<size_t> SK2_2_x=sort_indexes(S2, block*WORD);

	/*
	for (size_t x=0; x < block; x++)
	{
		std::cerr << SK1_2_x[x] << ", "  << S1[SK1_2_x[x]] << ", " << S1[x] << std::endl;
	/k/}

	for (size_t x=0; x < block; x++)
	{
		std::cerr << SK2_2_x[x] << ", "  << S2[SK2_2_x[x]] << ", " << S2[x] << std::endl;
	}
	*/
	
	for (size_t i=0; i < REP; i++) 
	{
		std::cerr << "i:" << i << std::endl;
		Pstates.rewind();
		ll=0;

		while (!Pstates.empty())
		{
			Pstates.uncompress(P1, P2, M1, M2);
			//TODO I think these can be dropped by moving them into Eigen...
			memset(p_arr, 0, sizeof(UINT)*WORD);

			if(SK1==WORD*block && SK2==WORD*block)
			{
				//std::cerr  << "reallocating\n";
				for (size_t x=0; x < WORD*block; x++)
				{
					UINT A=unif_L(rng);
					//poly[unif_L(rng)];
					//UINT B=gsl_ran_geometric (const gsl_rng * r, double p);
					double tB=c[A]+geom_L(rng);
					//UINT B==tB*Z || B==tB*A;
					//std::cerr << tB-c[A] << std::endl;

					UINT C=A+1;
					while (c[C++] < tB && C < Z-1);

					polyA[x]=A;
					polyB[x]=C;

					A=poly[A];
					double B=poly[C];

					//if( A > B) std::swap(A,B);
					//if (A == 28046 || B == 28046 )  std::cerr << A << ", " <<  B << ", " << B-A << std::endl;

					S1[x]= A >> 5;
					S2[x]= int(B) >> 5;
					S1k[x]= A % 0x20 ;
					S2k[x]= int(B) % 0x20 ;
					
				}
				SK1_2_x=sort_indexes(S1, block*WORD);
				SK2_2_x=sort_indexes(S2, block*WORD);
				SK1=0;
				SK2=0;
			}


			if (SK1 < WORD*block) 
			{
				while ( ll == S1[SK1_2_x[SK1]] )
				{
					size_t k=S1k[SK1_2_x[SK1]];
					//std::cerr << "L1c:" << SK1_2_x[SK1] << " of " << 32*block << ", SK1_2_x:" << SK1 << " of " <<  SK1_2_x.size() << ", c:" << polyA[SK1_2_x[SK1] ] << " of " << c.size() << std::endl;
					{
						{
							{
								L1c[ SK1_2_x[SK1] ]=c[ polyA[SK1_2_x[SK1] ] ];
								UINT z= SK1_2_x[SK1] >> 5;
								UINT zk= SK1_2_x[SK1] % 0x20;
								for ( size_t x = 0; x < N; x++)
								{
									if ( ( (M1[x] & Mask[k]) != 0) && ( (M2[x] & Mask[k]) != 0) )
									{		
										if ( ( P1[x] & Mask[k] ) != 0) 
											L11[ z*N*2+x ] |=   1UL << zk;
										else 
											L11[ z*N*2+x ] &= ~(1UL << zk);
	
										if ( ( P2[x] & Mask[k] ) != 0) 
											L12[ z*N*2+x ] |=   1UL << zk;
										else 
											L12[ z*N*2+x ] &= ~(1UL << zk);
	
										mL11[ z*N*2+x ] |= 1UL << zk;
										mL12[ z*N*2+x ] |= 1UL << zk;

										//mL11[ z*N*2+x ] = 0xFFFFFFFFF; //~(1UL << zk);
										//mL12[ z*N*2+x ] = 0xFFFFFFFFF; //~(1UL << zk);
									}
									
									else
									{
										mL11[x+K1*N] &=  ~(1UL << zk);
										mL12[x+K1*N] &=  ~(1UL << zk);
									}
								}
							}
						}
					}
					SK1++;
					if (SK1==WORD*block) break; 
				}
			}
			if (SK2 < WORD*block) 
			{
				while ( ll == S2[SK2_2_x[SK2]] )
				{
					size_t k=S2k[SK2_2_x[SK2]];
					//std::cerr << "L2c:" << S2[SK2_2_x[SK2]] << " of " << 32*block << ", SK2_2_x:" << SK2 << " of " <<  SK2_2_x.size() << ", c:" << polyB[SK2_2_x[SK2]] << " of " << c.size() << std::endl;
					{
						{
							{
								L2c[ SK2_2_x[SK2] ]=c[ polyB[SK2_2_x[SK2] ] ];
								UINT z= SK2_2_x[SK2] >> 5;
								UINT zk= SK2_2_x[SK2] % 0x20;
								for ( size_t x = 0; x < N; x++)
								{
									if ( (M1[x] & Mask[k]) != 0 && (M2[x] & Mask[k]) != 0 )
									{		
										if ( ( P1[x] & Mask[k] ) != 0) 
											L21[ z*N*2+x ] |=   1UL << zk;
										else 
											L21[ z*N*2+x ] &= ~(1UL << zk);
	
										if ( ( P2[x] & Mask[k] ) != 0) 
											L22[ z*N*2+x ] |=   1UL << zk;
										else 
											L22[ z*N*2+x ] &= ~(1UL << zk);
	
										mL21[ z*N*2+x ] |= 1UL << zk;
										mL22[ z*N*2+x ] |= 1UL << zk;

									}
									else
									{
										mL21[x+K2*N] &=   ~(1UL << K2k);
										mL22[x+K2*N] &=   ~(1UL << K2k);
									}
								}
							}
						}
					}
					SK2++;
					if (SK2==WORD*block) break; 
				}
			}
			ll++;
			if ( SK1 == WORD*block && SK2 == WORD*block ) break;
		}

		std::cerr << "Done\n";

		if (SK1==WORD*block && SK2 == WORD*block )
		{
			/*
			for (size_t x=0; x<block*WORD; x++)
				std::cerr << "SK1:" << SK1_2_x[x] << std::endl;
			for (size_t x=0; x<block*WORD; x++)
				std::cerr << "SK2:" << SK2_2_x[x] << std::endl;
			*/
			for (size_t x=0; x<2*block*N; x+=2*N)
			{
				for (size_t y=0; y<WORD; y++)
				{
					r[y]=double(L2c[ int(floor( double(x)/(2.*N)*WORD ) )+ y ]-L1c[ int(floor( double(x)/(2.*N)*WORD ) )+y ] );
					r[y]=(1.-exp(-( double(2.)*r[y]/100. ) ) )/2.;
					//r[y]=(1.-exp(-( double(2.*FOUR)*r[y]*double(map_size) )/(double(LEN) ) ) )/2.;
				//	std::cerr <<  x / (2*N) *WORD+y << ", " << L1c[x/(2*N)*WORD+y ] << ", " << L2c[x/(2*N)*WORD+y ] << ", " << r[y] << std::endl;
				}
				//std::cerr << "getting lambda " << std::endl;
				get_lambda(L1+x, L2+x, mL11+x, mL21+x, r, N, lambda, lambdaN);
				get_f(L1+x, mL11+x, N, f, fN);
				/*for (size_t y=0; y<WORD; y++)
				{
					std::cerr << r[y] << ", " << c2bin(r[y]) << ", " << lambda(x,y)/lambdaN(x,y) << std::endl;
				}*/
			}
			//return 0;
			K1=0;
			K2=0;
			K2k=0;
			K1k=0;
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
