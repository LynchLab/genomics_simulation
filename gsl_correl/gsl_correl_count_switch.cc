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

	uint32_t *count;

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

	//std::cerr << "allocating matrices..." << std::endl;

	//return 0;

	int N=Pstates.sample_size();
	int n=N;

	count=new uint32_t [N];

	for (int x=0; x<N; x++)
		count[x]=0;

	if(sub_matrix)
	{
		n=n/2;
	}
	int LEN=Pstates.genome_size()*32;

			  //		      MtM, MtH, HtH`
	//size_t mem_needed =5*n*sizeof(double)+3*n*n*sizeof(double)+2*n*2*sizeof(uint32_t) + 3*BLOCK*n*sizeof(double)+3*BLOCK*sizeof(double)+Pstates.buffer_size();

	//std::cerr << Pstates.buffer_size() << std::endl; //big

	//size_t mem_avail = system_memory();

	//std::cerr << "Total memory needed:" << mem_needed << "(" << float(mem_needed)/float(mem_avail) << ") of " << mem_avail << std::endl;

	//if (mem_needed > mem_avail) 
	//{
	//	std::cerr << "not enough memory" << std::endl;
	//}

	UINT p_arr[WORD];

	UINT *P1=new uint32_t [N*2];
	UINT *P2=P1+N;

	UINT *N1=new uint32_t [N*2];
	UINT *N2=N1+N;

	while (!Pstates.empty())
	{

		Pstates.uncompress(P1, P2, N1, N2);

		for (size_t k=0; k < WORD-1; k++)
		{
			for (int x=0; x<N; x++)
			{
				//if (x==0)
				//	std::cerr << x << ", " << k << ", " << ( ((P1[x] & Mask[k])!=0) != ((P1[x] & Mask[k+1])!=0) ) << ", " << ( ( (P2[x] & Mask[k])!=0) != ((P2[x] & Mask[k+1])!=0) ) << std::endl;
				count[x] += ( ((P1[x] & Mask[k])!=0) != ((P1[x] & Mask[k+1])!=0) ) + ( ( (P2[x] & Mask[k])!=0) != ((P2[x] & Mask[k+1])!=0) );
				//std::cerr << x << "\t" <<  int( (P1[x] & Mask[k] == 0) != (P1[x] & Mask[k+1] == 0) ) << std::endl;
			}
		}
	}


	for (int x=0; x<N; x++)	
	{
		std::cout << x << '\t' << count[x] << std::endl;
	}
}
