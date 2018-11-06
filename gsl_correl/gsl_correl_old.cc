#include "interface.h"
#include "map_file.h"
#include "state.h"
#include "sample_name.h"
#include "correl_data.h"

#include <cstring>
#include <sstream>
#include <tuple>
#include <map>
#include <fstream>

#include "Eigen/Core"

#define UINT uint32_t
#define WORD 32

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

	std::string names_file="", input_file="";
	int indX=-1, indY=-1, model=0;
	std::string namex="", namey="";

	Environment env;
	env.set_name("call_relatedness");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("A POS relatedness caller. Please direct questions to matthew.s.ackerman@gmail.com");

	env.optional_arg('x',"namey",  namex,      "please provide a number.", "number of individuals in the populations.");
	env.optional_arg('y',"namex",  namey,      "please provide a number.", "number of individuals in the populations.");
	env.optional_arg('m',"model",  model,      "please provide a number.", "model of the DoGE to use.");
	env.optional_arg('i',"input",  input_file,      "please provide a number.", "number of individuals in the populations.");
	env.positional_arg('n',"names",  names_file,      "please provide a number.", "number of individuals in the populations.");

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);

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
		std::cerr << Pstates.sample_size() << ", " << Pstates.genome_size() << std::endl;
	}

	int N(Pstates.sample_size());
	int LEN(Pstates.genome_size()*32);
	
	Eigen::MatrixXf G(LEN, N*2);

	UINT *dm=new UINT[N];
	UINT p_arr[WORD];
	UINT h_arr[WORD];

	float a=1, d=-1;

	switch (model)
	{
		case 0:
			a=1;
			d=0;
		break;
		case 1:
			a=0;
			d=1;
		break;
		case 2:
			a=1;
			d=1;
		break;
		case 3:
			a=-1;
			d=0;
		break;
		case 4:
			a=-1;
			d=1;
		break;
		case 5:
			a=-1;
			d=-1;
		break;
	}	

	size_t l=0;

	UINT *P=new uint32_t [N*2];
	UINT *P2=P+N;

	while (!Pstates.empty())
	{	
		Pstates.uncompress(P, P2);
		memset(p_arr, 0, sizeof(UINT)*WORD);
		memset(h_arr, 0, sizeof(UINT)*WORD);
		transposed_bit_sumN(P, 2*N, p_arr);
		transposed_het_sumN(P, N, h_arr);
		for (size_t k=0; k<WORD; k++)
		{
			float p=float(p_arr[k])/float(2*N);
			//std::cerr << p_arr[k] << std::endl;
			if (p>0 and p<1)
			{
				float q=1.-p;
	
				for (size_t x=0; x<N; x++)
					dm[x]=( (P[x] & (1 << k) ) !=0)+( (P2[x] & (1<<k) ) !=0);

				float mean_a=a*(2.*p-1.);
				float mean_d=d*float(h_arr[k])/float(N);
	
				float this_beta[3]={-a-mean_a, -mean_a, a-mean_a};
				float this_delta[3]={-mean_d, d-mean_d, -mean_d };
 		 
				for (size_t x=0; x<N; x++)
				{
					G(l,x)=this_beta[dm[x]];
					G(l,N+x)=this_delta[dm[x]];
				}
				l++;
			}
		}
		//std::cerr << "read " << std::endl;
	}
	delete [] dm;

	G.conservativeResize(l,N*2);
	
	Eigen::VectorXf gsum=G.colwise().sum();
	Eigen::VectorXf mua=gsum.head(N);
	Eigen::VectorXf mud=gsum.tail(N);

	//Eigen::MatrixXf Gp = 
	G.rowwise()-=(gsum.transpose()/l);

	Eigen::MatrixXf cov=G.transpose() * G;
	Eigen::MatrixXf A=cov.block(0, 0, N, N);
	Eigen::MatrixXf D=cov.block(N,N,N,N);
	Eigen::MatrixXf AD=cov.block(0,N,N,N)+cov.block(N,0,N,N);

	Eigen::MatrixXf S1=mua * mua.transpose();//N*A/A.trace();
	Eigen::MatrixXf S2=mud * mud.transpose();//N*A/A.trace();

	Eigen::MatrixXf S3=AD;//N*AD/AD.trace();
	Eigen::MatrixXf S4=A;//N*A/A.trace();
	Eigen::MatrixXf S5=D;//N*D/D.trace();

         char del='\t';
	 double l2=1;//l*l;
         for (size_t x=0; x<N; x++)
         {
		std::cout << mua(x) << del;
		std::cout << mud(x) << del;
                 for (size_t y=0; y<N; y++)
                 {
                         std::cout << S3(x,y)/l2 << del;
                 }
                 for (size_t y=0; y<N; y++)
                 {
                         std::cout << S4(x,y)/l2 << del;
                 }
                 for (size_t y=0; y<N-1; y++)
                 {
                         std::cout << S5(x,y)/l2 << del;
                 }
                 std::cout << S5(x,N-1)/l2 << std::endl;
         }

/*	
	COV=tG*G;
	matrix A=COV.block();
	D=COV.block();
	AD=COV.block();
*/
	//close();
}
