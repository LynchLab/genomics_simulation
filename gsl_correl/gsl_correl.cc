#include "interface.h"
#include "map_file.h"
#include "state.h"
#include "sample_name.h"
#include "relatedness_data.h"
#include "system.h"

#include <cstring>
#include <sstream>
#include <tuple>
#include <map>
#include <fstream>
#include <ctime>
#include <omp.h>

#include "matrix.h"
#include "types.h"

#include "Eigen/Core"

#define UINT uint32_t
#define WORD 32

//#define BLOCK 4096
#define BLOCK 1024
//#define BLOCK 32

static UINT Mask[WORD]={0x00000001, 0x00000002, 0x00000004, 0x00000008,
			0x00000010, 0x00000020, 0x00000040, 0x00000080,
			0x00000100, 0x00000200, 0x00000400, 0x00000800,
			0x00001000, 0x00002000, 0x00004000, 0x00008000,
			0x00010000, 0x00020000, 0x00040000, 0x00080000,
			0x00100000, 0x00200000, 0x00400000, 0x00800000,
			0x01000000, 0x02000000, 0x04000000, 0x08000000,
			0x10000000, 0x20000000, 0x40000000, 0x80000000};

inline void transposed_bit_sumN(const UINT *place, const UINT &N, UINT *sum)
{
	const UINT *end=place+N;
	const UINT *this_place=place;
	while (this_place!=end) {
		sum[0]+=( (*this_place & Mask[0])!=0 );
		sum[1]+=( (*this_place & Mask[1])!=0 );
		sum[2]+=( (*this_place & Mask[2])!=0 );
		sum[3]+=( (*this_place & Mask[3])!=0 );
		sum[4]+=( (*this_place & Mask[4])!=0 );
		sum[5]+=( (*this_place & Mask[5])!=0 );
		sum[6]+=( (*this_place & Mask[6])!=0 );
		sum[7]+=( (*this_place & Mask[7])!=0 );
		sum[8]+=( (*this_place & Mask[8])!=0 );
		sum[9]+=( (*this_place & Mask[9])!=0 );
		sum[10]+=( (*this_place & Mask[10])!=0 );
		sum[11]+=( (*this_place & Mask[11])!=0 );
		sum[12]+=( (*this_place & Mask[12])!=0 );
		sum[13]+=( (*this_place & Mask[13])!=0 );
		sum[14]+=( (*this_place & Mask[14])!=0 );
		sum[15]+=( (*this_place & Mask[15])!=0 );
		sum[16]+=( (*this_place & Mask[16])!=0 );
		sum[17]+=( (*this_place & Mask[17])!=0 );
		sum[18]+=( (*this_place & Mask[18])!=0 );
		sum[19]+=( (*this_place & Mask[19])!=0 );
		sum[20]+=( (*this_place & Mask[20])!=0 );
		sum[21]+=( (*this_place & Mask[21])!=0 );
		sum[22]+=( (*this_place & Mask[22])!=0 );
		sum[23]+=( (*this_place & Mask[23])!=0 );
		sum[24]+=( (*this_place & Mask[24])!=0 );
		sum[25]+=( (*this_place & Mask[25])!=0 );
		sum[26]+=( (*this_place & Mask[26])!=0 );
		sum[27]+=( (*this_place & Mask[27])!=0 );
		sum[28]+=( (*this_place & Mask[28])!=0 );
		sum[29]+=( (*this_place & Mask[29])!=0 );
		sum[30]+=( (*this_place & Mask[30])!=0 );
		sum[31]+=( (*this_place & Mask[31])!=0 );
		this_place++;
	}
}

inline void transposed_het_sumN(const UINT *place, const UINT &N, UINT *sum)
{
	const UINT *end=place+2*N;
	const UINT *this_place1=place;
	const UINT *this_place2=place+N;
	while (this_place2!=end) {
		for (int x=0; x<WORD; ++x)
		{
			sum[x]+=( (*this_place1 & Mask[x])!=(*(this_place2) & Mask[x]) );
		}
		this_place1++;
		this_place2++;
	}
}

int main (int argc, char **argv){

	std::cerr << __FILE__ << std::endl;

	std::string names_file="", input_file="";
	int indX=-1, indY=-1, model=0;
	std::string namex="", namey="";

	int a=0;
	int b=0;

	bool sub_matrix=false;
	bool binary=false;

	Environment env;
	env.set_name("call_relatedness");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("A POS relatedness caller. Please direct questions to matthew.s.ackerman@gmail.com");

	env.optional_arg('x',"namey",  namex,      "please provide a number.", "number of individuals in the populations.");
	env.optional_arg('y',"namex",  namey,      "please provide a number.", "number of individuals in the populations.");
	env.optional_arg('m',"model",  model,      "please provide a number.", "model of the DoGE to use.");
	env.optional_arg('i',"input",  input_file,      "please provide a number.", "number of individuals in the populations.");
	env.positional_arg('n',"names",  names_file,      "please provide a number.", "names of individuals in the populations.");
	env.optional_arg('A',"parta",  a,      "please provide a number.", "number of individuals in the populations.");
	env.optional_arg('B',"partb",  b,      "please provide a number.", "number of individuals in the populations.");
	env.flag(       'b',"binary",   &binary,        &flag_set,      "an error occurred while displaying the version message", "binary output");             //DONE


	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);

	//Eigen::initParallel();
	Eigen::setNbThreads(4);
	std::cerr << "using " << Eigen::nbThreads( ) << " threads.\n";

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
		std::cerr << "Sample: " << Pstates.sample_size() << ", " << " genome: "<<  Pstates.genome_size() << std::endl;
	}

	std::cerr << "allocating matrices..." << std::endl;

	int N=Pstates.sample_size();
	int n=N;
	if(sub_matrix)
	{
		n=n/2;
	}
	int LEN=Pstates.genome_size()*32;

			  //		      MtM, MtH, HtH`
	size_t mem_needed =5*n*sizeof(double)+3*n*n*sizeof(double)+2*n*2*sizeof(uint32_t) + 3*BLOCK*n*sizeof(double)+3*BLOCK*sizeof(double)+Pstates.buffer_size();
	std::cerr << 5*n*sizeof(double) << std::endl;  
	std::cerr << 3*n*n*sizeof(double) << std::endl; //big
	std::cerr << 2*n*2*sizeof(uint32_t) << std::endl; 
	std::cerr << 3*BLOCK*n*sizeof(double) << std::endl; //big
	std::cerr << 3*BLOCK*sizeof(double) << std::endl;
	std::cerr << Pstates.buffer_size() << std::endl; //big
	size_t mem_avail = system_memory();
	std::cerr << "Total memory needed:" << mem_needed << "(" << float(mem_needed)/float(mem_avail) << ") of " << mem_avail << std::endl;
	if (mem_needed > mem_avail) 
	{
		std::cerr << "not enough memory" << std::endl;
	}


	VECTOR s1=ARRAY1::Zero(n);
	VECTOR s2=ARRAY1::Zero(n);
	VECTOR ll=ARRAY1::Zero(n);

	VECTOR fs(N);
	VECTOR dm(N);

	MATRIX S4=ARRAY2::Zero(n, n);
	MATRIX S5=ARRAY2::Zero(n, n);
	MATRIX S6=ARRAY2::Zero(n, n);

	UINT p_arr[WORD];

	UINT *P1=new uint32_t [N*2];
	UINT *P2=P1+N;

	UINT *N1=new uint32_t [N*2];
	UINT *N2=N1+N;

	MATRIX M(BLOCK,N), H(BLOCK,N), Q=MATRIX::Constant(BLOCK, N, 1);
	VECTOR EM(BLOCK), EH(BLOCK), EQ(BLOCK), F(BLOCK), PQ(BLOCK);


	double f=0, pq=0;

	clock_t t1, t2, t3, s;
	size_t K=0;
	s=clock();

	while (!Pstates.empty())
	{

		Pstates.uncompress(P1, P2, N1, N2);

		memset(p_arr, 0, sizeof(UINT)*WORD);

		transposed_bit_sumN(P1, 2*N, p_arr);
	
		for (size_t k=0; k<WORD; k++)
		{
			double p=double(p_arr[k])/double(2*N);
			if (p>0 and p<1)
			{
				{
				for (size_t x=0; x<N; x++)
				{
					Q(K,x)=( (N1[x] & (1 << k) !=0 ) & (N2[x] & (1<<k) != 0 ) );
					if ( Q(K,x) )
					{
						M(K,x)=( (P1[x] & (1<<k) ) !=0)+( (P2[x] & (1 << k) ) !=0);
						H(K,x)=( (P1[x] & (1 << k) )!=(P2[x] & (1<<k) ) );
					} else {
						M(K,x)=0;
						H(K,x)=0;
					}
				}
				if  ( 2*Q.block(K,0,1,N).sum()!=M.block(K,0,1,N).sum() ) K++;
				}
				if(K==BLOCK)
				{
					EQ = Q.rowwise().sum();

					EH = H.rowwise().sum().array() / EQ.array();
					EM = M.rowwise().sum().array() / EQ.array();

					PQ = (EM.array()/2.)*(1.-EM.array()/2. );
					F = ( 1.-EH.array()/(2.*(PQ.array() ) ) );

					/*
					for (int x=0; x<50; x++)
					{
					std::cerr << 2*Q.block(x,0,1,N).sum() << "/" << M.block(x,0,1,N).sum() << ":";
					std::cerr << EH(x,0) << ", ";
					std::cerr << EM(x,0) << ", ";
					std::cerr << PQ(x,0) << ", ";
					std::cerr << F(x,0) << std::endl ;
					}
					*/

					f += F.sum();
					pq += BLOCK; 

					//PQ.sum();
					//EH = (H.rowwise().sum() );
					//EM = (M.rowwise().sum() );

					M = center(M, EM).cwiseProduct(Q);
					H = center(H, EH).cwiseProduct(Q);

					if (sub_matrix){
						MATRIX M1=M.block(0, a*n, BLOCK, n);
						MATRIX M2=M.block(0, b*n, BLOCK, n);
						MATRIX H1=H.block(0, a*n, BLOCK, n);
						MATRIX H2=H.block(0, b*n, BLOCK, n);

						S4 += ( t(H1)*M2+t(M1)*H2 );
						S5 += 2.*( t(M1)*M2 );
						S6 += 2.*( t(H1)*H2 );
					} else {
						S4 += ( t(H)*M+t(M)*H );
						S5 += 2.*( t(M)*M );
						S6 += 2.*( t(H)*H );
					}


					K=0;
	
					if (sub_matrix)
					{
						s1+=M.block(0,a*n,BLOCK,n).colwise().sum();
						s2+=H.block(0,a*n,BLOCK,n).colwise().sum();
						ll+=Q.block(0,a*n,BLOCK,n).colwise().sum();
					}  else {
						s1+=M.colwise().sum();
						s2+=H.colwise().sum();
						ll+=Q.colwise().sum();
					}
					t2=clock();
					std::cerr << BLOCK << ", t1:" << float(BLOCK)/float(t2-t3) << ", t2:" << float(t2-t3) << ", F=" << f/pq << std::endl;
					t3=clock();
				}
			}
		}
		//if(Pstates.genome_size() < 1000) break;
	}

	if(K!=0)
	{
		std::cerr << K << std::endl;
		MATRIX tM=M.block(0,0,K,N);
		MATRIX tH=H.block(0,0,K,N);
		MATRIX tQ=Q.block(0,0,K,N);

		M=tM;
		H=tH;
		Q=tQ;

		EQ = Q.rowwise().sum();
		EH = EQ.asDiagonal().inverse() * (H.rowwise().sum() );
		EM = EQ.asDiagonal().inverse() * (M.rowwise().sum() );

		M=center(M, EM).cwiseProduct(Q);
		H=center(H, EH).cwiseProduct(Q);

		if (sub_matrix){
			MATRIX M1=M.block(0, a*n, K, n);
			MATRIX M2=M.block(0, b*n, K, n);
			MATRIX H1=H.block(0, a*n, K, n);
			MATRIX H2=H.block(0, b*n, K, n);

			S4 += ( t(H1)*M2+t(M1)*H2 );
			S5 += ( 2.*t(M1)*M2 );
			S6 += ( 2.*t(H1)*H2 );
		} else {
			S4 += ( t(H)*M+t(M)*H );
			S5 += ( 2.*t(M)*M );
			S6 += ( 2.*t(H)*H );
		}

		K=0;

		if (sub_matrix)
		{
			s1+=M.block(0,a*n,BLOCK,n).colwise().sum();
			s2+=H.block(0,a*n,BLOCK,n).colwise().sum();
			ll+=Q.block(0,a*n,BLOCK,n).colwise().sum();
		} else {
			s1+=M.colwise().sum();
			s2+=H.colwise().sum();
			ll+=Q.colwise().sum();
		}
	}
	
	M.resize(0,0);
	H.resize(0,0);
	Q.resize(0,0);

	//double S5T=(S5.trace()-(s1*t(s1) ).trace() )/double(N);
	//double S6T=(S6.trace()-(s2*t(s2) ).trace() )/double(N);

	Relatedness rel(n);
	Flat_file <Relatedness> rel_out;

	if (binary)
		rel_out.open(WRITE | BINARY );
	else
		rel_out.open(WRITE);

	VECTOR l1=ll;
	MATRIX denom = ll*t(ll);

	for (int x=0; x<N; x++)	
	{
		for (int y=x; y<N; y++)	
		{
			float t=std::min(ll(x), ll(y) );
			std::cerr << x << "t:" << t << std::endl;
			denom(x,y) = 1./(2.*(N-1) ); //t / (2.*(t-1.) );
			denom(y,x) = 1./(2.*(N-1) ); //t / (2.*(t-1.) );
		}
	}

	for (int x=0; x<N; x++)	
	{
		std::cerr << "Meep Meep:" << x << ":" << l1(x) << std::endl;
		l1(x) =  1./ sqrt( 2*(N-1)*l1(x) ); // sqrt(N);// l1(x) ;
	}

	std::cerr << "ll=" << ll.mean() << ", l1=" << l1.mean() << ", F=" << f/pq << std::endl;

	rel.MtH = S4.cwiseProduct(denom);
	rel.MtM = S5.cwiseProduct(denom);
	rel.HtH = S6.cwiseProduct(denom);

	rel.Mt1 = s1.cwiseProduct(l1);
	rel.Ht1 = s2.cwiseProduct(l1);

	rel.set_mean_sites(ll.mean() );
	
	rel_out.write_header(rel);
	rel_out.write(rel);
	rel_out.close();
}
