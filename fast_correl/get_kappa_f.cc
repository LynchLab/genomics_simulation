#include "circular_list.h"
#include "interface.h"
#include "map_file.h"
#include "state.h"
#include "sample_name.h"
#include "relatedness_data.h"

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cstring>
#include <sstream>
#include <tuple>
#include <map>
#include <fstream>

static uint32_t Mask[32]={0x00000001, 0x00000002, 0x00000004, 0x00000008,
			0x00000010, 0x00000020, 0x00000040, 0x00000080,
			0x00000100, 0x00000200, 0x00000400, 0x00000800,
			0x00001000, 0x00002000, 0x00004000, 0x00008000,
			0x00010000, 0x00020000, 0x00040000, 0x00080000,
			0x00100000, 0x00200000, 0x00400000, 0x00800000,
			0x01000000, 0x02000000, 0x04000000, 0x08000000,
			0x10000000, 0x20000000, 0x40000000, 0x80000000};

inline int get(const uint32_t *place, const uint32_t &x, const uint32_t &bit)
{
	if ((*(place+x) & Mask[bit])!=0 ){
		return 1;
	} else  {
		return 0;
	}
}

inline void fast_bit_sumN(const uint32_t *place, const uint32_t &N, double *sum)
{
	const uint32_t *end=place+N;
	const uint32_t *this_place=place;
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
	sum[0]/=double(N);
	sum[1]/=double(N);
	sum[2]/=double(N);
	sum[3]/=double(N);
	sum[4]/=double(N);
	sum[5]/=double(N);
	sum[6]/=double(N);
	sum[7]/=double(N);
	sum[8]/=double(N);
	sum[9]/=double(N);
	sum[10]/=double(N);
	sum[11]/=double(N);
	sum[12]/=double(N);
	sum[13]/=double(N);
	sum[14]/=double(N);
	sum[15]/=double(N);
	sum[16]/=double(N);
	sum[17]/=double(N);
	sum[18]/=double(N);
	sum[19]/=double(N);
	sum[20]/=double(N);
	sum[21]/=double(N);
	sum[22]/=double(N);
	sum[23]/=double(N);
	sum[24]/=double(N);
	sum[25]/=double(N);
	sum[26]/=double(N);
	sum[27]/=double(N);
	sum[28]/=double(N);
	sum[29]/=double(N);
	sum[30]/=double(N);
	sum[31]/=double(N);
}

inline void fast_het_sumN(const uint32_t *place, const uint32_t &N, double *sum)
{
	const uint32_t *end=place+2*N;
	const uint32_t *this_place=place;
	const uint32_t *that_place=place+N;

	while (that_place!=end) {
		for (size_t x=0; x<32; ++x) sum[x]+=( (*this_place & Mask[x])!=(*that_place & Mask[x]) );
		++this_place;
		++that_place;
	}
	for (size_t x=0; x<32; ++x) sum[x]/=double(N);
}

int main (int argc, char **argv){


	std::string names_file="";
	int indX=-1, indY=-1;
	std::string namex="", namey="";

	Environment env;
	env.set_name("call_relatedness");
        env.set_version(VERSION);
        env.set_author("Matthew Ackerman");
        env.set_description("A POS relatedness caller. Please direct questions to matthew.s.ackerman@gmail.com");
        env.positional_arg('n',"names",  names_file,      "please provide a number.", "number of individuals in the populations.");

	State Pstates;
 	{
		Flat_file <State> state_file;
		state_file.open(READ);
		Pstates=state_file.read_header();
		state_file.read(Pstates);
		std::cerr << "done.\n";
		state_file.close();
		std::cerr << Pstates.sample_size() << ", " << Pstates.genome_size() << std::endl;
	}

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);

	Flat_file <Sample_name> in_names;
	Sample_name names;
	in_names.open(names_file.c_str(), std::ios::in);
	names=in_names.read_header();
	while (in_names.read(names).table_is_open() );
	in_names.close();

	uint32_t N=names.sample_names.size();

	uint32_t BLOCK_SIZE=2*N;

	uint32_t *P1=new uint32_t [N*2];
	uint32_t *P2=P1+N;

	Pstates.rewind();

	size_t size=Pstates.genome_size();
			
	double kappa=0, f=0;
	double polyK=0, polyF=0;
	for (size_t l=0;l<size;++l)
	{

		Pstates.uncompress(P1, P2);

		size_t x, y;
		double p[32]={0};
		double h[32]={0};
		fast_bit_sumN(P1, 2*N, p);
		fast_het_sumN(P1, N, h);
		for (size_t k=0; k<32; k++)
		{
			if (p[k]<1. && p[k]>0.)
			{
				kappa+=(1./p[k]+1./(1.-p[k])-3.)*p[k]*(1.-p[k]);
				f+=(1.-h[k]/(2.*p[k]*(1.-p[k]) ) );
				polyF+=1;
				polyK+=p[k]*(1.-p[x]);
			}
		}
	}
	std::cout << "kappa, f, count" << std::endl;
	std::cout << kappa/double(polyK) << ", " << f/double(polyF) << ", " << polyF << std::endl;
	//close();
	return 0;
}
