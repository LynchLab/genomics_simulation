#include "interface.h"
#include "map_file.h"
#include "state.h"
#include "sample_name.h"
#include "correl_data.h"
#include "relatedness_data.h"

#include <cstring>
#include <sstream>
#include <tuple>
#include <map>
#include <fstream>

#include "Eigen/Core"
#include "types.h"
#include "matrix.h"

#define UINT uint32_t
#define WORD 32

MATRIX makev(const size_t &n)
{
	MATRIX I=MATRIX::Identity(n,n); 
	MATRIX J=MATRIX::Constant(n,n, -1./(float(n)-1.) );
	return float(n)/(float(n)-1)*I+J;
}

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
	TYPE Z=1.;
	int indX=-1, indY=-1, model=0;
	std::string namex="", namey="";

	TYPE ma=0, md=0, va=1, vd=1, sad=0;

	Environment env;
	env.set_name("call_relatedness");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("A POS relatedness caller. Please direct questions to matthew.s.ackerman@gmail.com");

	env.optional_arg('A',"ma",  ma,      "please provide a number.", "mean of a.");
	env.optional_arg('D',"md",  md,      "please provide a number.", "mean of d.");
	env.optional_arg('a',"va",  va,      "please provide a number.", "var of d.");
	env.optional_arg('d',"vd",  vd,      "please provide a number.", "var of a.");
	env.optional_arg('e',"sad",  sad,      "please provide a number.", "a d covariance.");
	env.optional_arg('z',"vz",  Z,      "please provide a number.", "Phenotypic variance.");


	env.optional_arg('x',"namey",  namex,      "please provide a number.", "number of individuals in the populations.");
	env.optional_arg('y',"namex",  namey,      "please provide a number.", "number of individuals in the populations.");
	env.optional_arg('m',"model",  model,      "please provide a number.", "model of the DoGE to use.");
	env.optional_arg('i',"input",  input_file,      "please provide a number.", "number of individuals in the populations.");
	env.positional_arg('n',"names",  names_file,      "please provide a number.", "number of individuals in the populations.");

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);

	TYPE mud2=powf(md, 2);
	TYPE mua2=powf(ma, 2);

	TYPE sa=va, sd=vd;

	/*if (va > 0)
		sa=sqrt(va);
	else 
		sa=0;

	if (vd > 0)
		sd=sqrt(vd);
	else 
		sd=0;*/

	std::cerr << "Input: " << ma << ", " << md << ", " << sa << ", " << sd << ", " << sad << std::endl;

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
	
	VECTOR s1=ARRAY1::Zero(N);
	VECTOR s2=ARRAY1::Zero(N);

	VECTOR mua=ARRAY1::Zero(N);
	VECTOR mud=ARRAY1::Zero(N);

	VECTOR fs(N);
	VECTOR dm(N);

	MATRIX S4=ARRAY2::Zero(N, N);
	MATRIX S5=ARRAY2::Zero(N, N);
	MATRIX S6=ARRAY2::Zero(N, N);

	MATRIX SA=ARRAY2::Zero(N, N);
	MATRIX SD=ARRAY2::Zero(N, N);
	MATRIX SAD=ARRAY2::Zero(N, N);

	UINT p_arr[WORD];
	//UINT h_arr[WORD];

	size_t ll=0;

	UINT *P=new uint32_t [N*2];
	UINT *P2=P+N;

	VECTOR h1(N), h2(N);// h3(N), h4(N);

	while (!Pstates.empty())
	{	
		Pstates.uncompress(P, P2);
		memset(p_arr, 0, sizeof(UINT)*WORD);
		//memset(h_arr, 0, sizeof(UINT)*WORD);
		transposed_bit_sumN(P, 2*N, p_arr);
		//transposed_het_sumN(P, N, h_arr);
		for (size_t k=0; k<WORD; k++)
		{
			TYPE p=TYPE(p_arr[k])/TYPE(2*N);
			if (p>0 and p<1)
			{
				TYPE q=1.-p;
	
				for (size_t x=0; x<N; x++) 
				{
					h1[x] = ( ( P[x] & (1 << k ) ) !=0 );
					h2[x] = ( ( P2[x] & (1 << k ) ) !=0 );
					dm[x] = ( ( P[x] & (1 << k ) ) != ( P2[x] & (1 << k ) ) );
				}

				TYPE dmbar = TYPE( dm.sum() )/TYPE( N );

				h1 = center_scalar( h1, p );
				h2 = center_scalar( h2, p );

				dm = center_scalar(dm, dmbar);

				h1 = h1+h2;

				s1+=VECTOR(h1)*ma;
				s2+=VECTOR(dm)*md;

				//ma*h1+md*dm
				mua+=VECTOR(h1*(ma+md*(q-p) ) );
				mud+=VECTOR(md*dm-h1*md*(q-p) );

				ll++;
			}
		}
	}

	MATRIX T(N,N);
	MATRIX G(N,N);
	MATRIX D(N,N);
//	MATRIX Fs=( (t(s2.replicate(1,N) ) ).array()+s2.replicate(1,N).array() )/TYPE(ll);


	TYPE l1=sqrt( ll ) / sqrt(2.*ll*(ll-1.) );

	s1=s1*l1;
	s2=s2*l1;

	mua=mua*l1;
	mud=mud*l1;

	Pstates.rewind();

	TYPE A1=0, B1=0, B2=0, B3=0, C1=0, C2=0;


	int count=0;	
	while (!Pstates.empty())
	{
		Pstates.uncompress(P, P2);
		memset(p_arr, 0, sizeof(UINT)*WORD);
		//memset(h_arr, 0, sizeof(UINT)*WORD);
		transposed_bit_sumN(P, 2*N, p_arr);
		//transposed_het_sumN(P, N, h_arr);
		for (size_t k=0; k<WORD; k++)
		{
			TYPE p=TYPE(p_arr[k])/TYPE(2*N);
			if (p>0 and p<1)
			{
				TYPE q=1.-p;
	
				for (size_t x=0; x<N; x++)
				{
					h1[x]=( (P[x] & (1<<k) ) !=0);
					h2[x]=( (P2[x] & (1 << k) ) !=0);
					dm[x]=( (P[x] & (1 << k) )!=(P2[x] & (1<<k) ) );
				}

				TYPE dmbar=TYPE(dm.sum())/TYPE(N);

				h1=center_scalar(h1, p);
				h2=center_scalar(h2, p);

				dm=center_scalar(dm, dmbar);

				h1=h1+h2;

				T=h1*t(h1);
				D=dm*t(dm);
				G=h1*t(dm)+dm*t(h1);

				S4.noalias() += ( T*( sa+sd*powf(q-p, 2)+2*sad*(q-p) ) );
				S5.noalias() += ( sd*D+T*sd*powf(q-p, 2)-G*(q-p)*sd );
				S6.noalias() += ( G*( (q-p)*sd+sad )+T*(-2*powf(q-p, 2)*sd -2*sad*(q-p) ) );

				SA.noalias() += T*sa;
				SD.noalias() += D*sd;
				SAD.noalias() += G*sad;

				//std::cerr << p << ", " << T.trace() << ", " << D.trace() << ", " << G.trace() << std::endl;

				TYPE ET = 2.*N*1.*(p-powf(p, 2) );
				TYPE ED = 2.*N*2.*(1./16.-powf(0.5-p, 4) );
				TYPE EG = 2.*N*4.*(powf(p, 3)-1.5*powf(p, 2)+0.5*p );

				A1 += ET * ( (sa)+(sd)*powf(q-p, 2)+2*(sad)*(q-p) );
				B1 += ED * (sd);
				B2 += ET * (sd)*powf(q-p, 2);
                                B3 += EG * (q-p)*(sd);
				C1 += EG * ( (q-p)*(sd)+(sad) );
				C2 += ET * (-2*powf(q-p, 2)*(sd) -2*(sad)*(q-p) );
			
				count+=1.;	
			}
		}
	}
	
	std::cerr << A1 << ", " << B1 << ", " << B2 << ", " << B3 << ", " << C1 << ", " << C2 << std::endl;

	char del='\t';

	MATRIX V = makev(N);

	TYPE denom = TYPE(ll)/( (TYPE(ll)-1.) );

	std::cerr << "denom : " << denom << std::endl;

	//Wait, shouldn't this be Beta and Delta or something?

	TYPE S4T = S4.trace()*denom+t(mua)*V*mua;
	TYPE S5T = S5.trace()*denom+t(mud)*V*mud;
	TYPE S6T = S6.trace()*denom+(t(mua)*V*mud+t(mud)*V*mua)(0,0);

	TYPE SAT = SA.trace()*denom+t(s1)*V*s1;
	TYPE SDT = SD.trace()*denom+t(s2)*V*s2;
	TYPE SADT = SAD.trace()*denom+(t(s1)*V*s2+t(s2)*V*s1)(0,0);

	std::cerr << A1/N/Z << std::endl;
	std::cerr << (B1+B2-B3)/N/Z << std::endl;
	std::cerr << (C1+C2)/N/Z << std::endl;

	Relatedness rel(N);
	Flat_file <Relatedness> rel_out;
	rel_out.open(WRITE);

	rel.Mt1 = mua;
	rel.Ht1 = mud;

	rel.MtM = S4/S4T;
	rel.HtH = S5/S5T;
	rel.MtH = S6/S6T;

	rel_out.write_header(rel);
	rel_out.write(rel);
	rel_out.close();

	//std::cerr << "E(A)\t" << s1 << std::endl;
	//std::cerr << "E(d)\t" << s2 << std::endl;

	//std::cerr << "E(beta)\t" << mua << std::endl;
	//std::cerr << "E(delta)\t" << mud << std::endl;

	std::cerr << "Var(A)/Var(Z)\t" << SAT/N/Z << std::endl;
	std::cerr << "Var(D)/Var(Z)\t" << SDT/N/Z << std::endl;
	std::cerr << "2*Cov(A,D)/Var(Z)\t" << SADT/N/Z << std::endl;

	std::cerr << "Var(Beta)/Var(Z)\t" << S4T/N/Z << std::endl; 
	std::cerr << "Var(Delta)/Var(Z)\t" << S5T/N/Z << std::endl;
	std::cerr << "2*Cov(Beta,Delta)/Var(Z)\t" << S6T/N/Z << std::endl;

	std::cerr << "Var(A2)/Var(Z)\t" << SA.trace()/N/Z << std::endl; 
	std::cerr << "Var(D2)/Var(Z)\t" << SD.trace()/N/Z << std::endl;
}
