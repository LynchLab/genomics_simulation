#include "circular_list.h"
#include "interface.h"
#include "map_file.h"
#include "state.h"
#include "sample_name.h"
#include "correl_data.h"
#include "triu_index.h"

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cstring>
#include <sstream>
#include <tuple>
#include <map>
#include <fstream>

#define BUFFER_SIZE2	1000
#define WORD	32
#define UINT	uint32_t
#define STEP	N9 << 1	
#define STEP_F	N << 1	
#define BLOCK	256

uint64_t LIM;
uint64_t N;
uint64_t N9;
uint64_t SIZE;
uint64_t XMIN;

inline size_t hash(const Triangular_index &T, const Triangular_index &t, const uint16_t &d, const uint16_t &i, const uint16_t &j)
{

	return ( (T.get_k()-t.get_k() ) << 1)*N9+9*d+3*i+j;
}

inline size_t hash2(const uint16_t &d, const uint16_t &i, const uint16_t &j)
{

	return 9*d+3*i+j;
}

inline size_t hash_F(const Triangular_index &T, const Triangular_index &t, const uint16_t &d)
{

	return ( (T.get_k()-t.get_k() ) << 1)*N+d;
}

inline size_t hash_F2(const uint16_t &d)
{

	return d;
}

//STATE FILE
/*			0x0000000000000001, 0x0000000000000002, 0x0000000000000004, 0x0000000000000008, 
			0x0000000000000010, 
			0x0000000000000100, 
			0x0000000000001000, 
			0x0000000000010000, 
			0x0000000000100000, 
			0x0000000001000000, 
			0x0000000010000000, 
			0x0000000100000000, 
			0x0000001000000000, 
			0x0000010000000000, 
			0x0000100000000000, 
			0x0001000000000000, 
			0x0010000000000000, 
			0x0100000000000000, 
			0x1000000000000000, 
}*/

static UINT Mask[WORD]={0x00000001, 0x00000002, 0x00000004, 0x00000008,
			0x00000010, 0x00000020, 0x00000040, 0x00000080,
			0x00000100, 0x00000200, 0x00000400, 0x00000800,
			0x00001000, 0x00002000, 0x00004000, 0x00008000,
			0x00010000, 0x00020000, 0x00040000, 0x00080000,
			0x00100000, 0x00200000, 0x00400000, 0x00800000,
			0x01000000, 0x02000000, 0x04000000, 0x08000000,
			0x10000000, 0x20000000, 0x40000000, 0x80000000};

inline int get(const UINT *place, const UINT &x, const UINT &bit)
{
	if ((*(place+x) & Mask[bit])!=0 ){
		return 1;
	} else  {
		return 0;
	}
}

inline void fast_bit_sumN(const UINT *place, const UINT &N, UINT *sum)
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
		/*for (int x=0; x<WORD; ++x) 
		{
			std::cerr << "C " << sum[x] << std::endl;
		}*/
		this_place++;
	}
}

inline void fast_het_sumN(const UINT *place, const UINT &N, UINT *sum)
{
	const UINT *end=place+2*N;
	const UINT *this_place1=place;
	const UINT *this_place2=place+N;
	while (this_place2!=end) {
		for (int x=0; x<WORD; ++x)
		{
			sum[x]+=( (*this_place1 & Mask[x])!=(*(this_place2) & Mask[x]) );
			//std::cerr << "H " << sum[x] << std::endl;
		}
		this_place1++;
		this_place2++;
	}
}

int main (int argc, char **argv){

	std::string names_file="", input_file="";
	int indX=-1, indY=-1;
	std::string namex="", namey="";

	Environment env;
	env.set_name("call_relatedness");
        env.set_version(VERSION);
        env.set_author("Matthew Ackerman");
        env.set_description("A POS relatedness caller. Please direct questions to matthew.s.ackerman@gmail.com");

        env.optional_arg('x',"namey",  namex,      "please provide a number.", "number of individuals in the populations.");
        env.optional_arg('y',"namex",  namey,      "please provide a number.", "number of individuals in the populations.");
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


	Flat_file <Sample_name> in_names;
	Sample_name names;
	in_names.open(names_file.c_str(), std::ios::in);
	names=in_names.read_header();
	while (in_names.read(names).table_is_open() );

	N=names.sample_names.size();
	N9=N*9;

	if (namex!="") {
		indX=find(names.sample_names.begin(), names.sample_names.end(), namex) - names.sample_names.begin();
	}
	if (namey!="") {
		indY=find(names.sample_names.begin(), names.sample_names.end(), namey) - names.sample_names.begin();
	}

	//Maximum number of uint16_t to allocate
	//4 Gb =   2000000000
	SIZE=9000000;
	UINT BLOCK_SIZE=2*N;

	UINT *P1=new UINT [N*2*BLOCK];
	UINT *P2=P1+N;
	UINT P_step=N*2;

	//Need to be private
	Triangular_index T_MIN(N), T_MAX(N);
	Triangular_index T(N);
	uint16_t *D=new uint16_t [SIZE*BLOCK];
	uint16_t *F=new uint16_t [SIZE*BLOCK/9];

	Correl buffer_rel[BUFFER_SIZE2];        

	Flat_file <Correl> rel_file;
	rel_file.open(WRITE);
	rel_file.write_header(buffer_rel[0]);        

	size_t J=std::min(BUFFER_SIZE2, int(T.size()-T.get_k()) );
	Pstates.rewind();

	while (J>0)
	{
	std::cerr << "start:" << J << std::endl;

	for (size_t z=0; z<J; ++z)
	{
		gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (pow(N*2-1, 1), 2);
		gsl_matrix *X = gsl_matrix_alloc (pow(N*2-1, 1), 2);
		gsl_matrix *XINV = gsl_matrix_alloc (pow(N*2-1, 1), 2);
		gsl_vector *Y = gsl_vector_alloc (pow(N*2-1, 1) );
		gsl_vector *w = gsl_vector_alloc (pow(N*2-1, 1) );
		gsl_vector *c = gsl_vector_alloc (2);
		gsl_matrix *cov = gsl_matrix_alloc (2, 2);

		if (indX!=-1 or indY!=-1)
		{
			T.set(indX, indY);
		}
	
		if (T>=T_MAX ) 
		{
			State_stream state_ptr;//=
			Pstates.set_stream(state_ptr);
			memset(D, 0, sizeof(uint16_t)*SIZE);
			memset(F, 0, sizeof(uint16_t)*SIZE/9);
			T_MIN=T;
			T_MAX=T+SIZE/(2*9*N);
			if (T_MAX.get_k() > T.size() ) T_MAX.set(N-1, N);
			//CHANGE TO LOOP
			size_t size=Pstates.genome_size();
			
			//UINT (T_MAX.get_k()-T_MIN.get_k() )*32		
	
			//size_t inc_size=SIZE/(2*9*N)*32;
			//size_t *inc=new size_t [inc_size]{0};
			clock_t time1, time2;

			//#pragma omp parallel for private(D)
			for (size_t l=0;l<size;l+=BLOCK)
			{
#ifdef DEBUG
				time1=clock();
#endif
				UINT this_block=BLOCK;
				if(size-l<BLOCK)
					this_block=size-l;
				
				{	
				UINT *P1_ptr=P1;
				UINT *P2_ptr=P2;

				for (size_t ll=0;ll<this_block;++ll)
				{
					Pstates.uncompress(P1_ptr, P2_ptr, state_ptr);
					P1_ptr+=P_step;
					P2_ptr+=P_step;
				}
				}

				UINT d1[WORD*BLOCK];
				memset(d1, 0, WORD*BLOCK*sizeof(UINT) );
				UINT h1[WORD*BLOCK];
				memset(h1, 0, WORD*BLOCK*sizeof(UINT) );

				UINT d_step=WORD;
				UINT h_step=WORD;

				{
				UINT *P1_ptr=P1;
				UINT *d1_ptr=d1;
				UINT *h1_ptr=h1;
				for (size_t ll=0;ll<this_block;++ll)
				{
					fast_bit_sumN(P1_ptr, 2*N, d1_ptr);
					fast_het_sumN(P1_ptr, N, h1_ptr);
					d1_ptr+=d_step;
					h1_ptr+=h_step;
					P1_ptr+=P_step;
				}
				}
				
				#pragma omp parallel for 
				for (Triangular_index t=T_MIN; t<T_MAX; ++t)
				{ 
					size_t x, y;
					t.get_xy(x,y);
					uint16_t *D_ptr=D+(t.get_k()-T_MIN.get_k())*(STEP);
					uint16_t *F_ptr=F+(t.get_k()-T_MIN.get_k())*(STEP_F);

					UINT *h1_ptr=h1;
					UINT *d1_ptr=d1;
					UINT *P1_ptr=P1;
					UINT *P2_ptr=P2;

					for (size_t ll=0;ll<this_block;++ll)
					{
						for (size_t k=0; k<WORD; ++k)
						{
							D_ptr[hash2(d1_ptr[k], get(P1_ptr, x, k)+get(P2_ptr, x, k), get(P1_ptr, y, k)+get(P2_ptr, y, k))]++;
							F_ptr[hash_F2(d1_ptr[k])]+=h1_ptr[k];
							//D[hash(t, T_MIN, d1_ptr[k], get(P1_ptr, x, k)+get(P2_ptr, x, k), get(P1_ptr, y, k)+get(P2_ptr, y, k))]++;
						}
						h1_ptr+=h_step;
						d1_ptr+=d_step;
						P1_ptr+=P_step;
						P2_ptr+=P_step;
					}
				}
#ifdef DEBUG
				time2=clock();
				std::cerr << l << " : of " << size << " : " << (SIZE*BLOCK)/double(time2-time1) << std::endl;
#endif
			}
			//delete [] inc;
		}

		double k_w, k, f_p, f_w, f_X=0, f_X_w=0, f_Y=0, f_Y_w=0, Theta=0, Theta_w=0, gamma_XY=0, gamma_XY_w=0, gamma_YX=0, gamma_YX_w=0, Z=0, Z_w=0;
		double mu, k1, k2;

		double chisq;

		Theta_w=0;
		k=0;
		f_p=0;
		f_w=0;
		k_w=0;
		Theta=0;
		f_X=0;
		f_X_w=0;
		f_Y=0;
		f_Y_w=0;
		gamma_XY=0;
		gamma_XY_w=0;
		gamma_YX=0;
		gamma_YX_w=0;
		size_t c2_min = 0;

		for (size_t c2=(1+c2_min); c2<((N-1)*2-c2_min); c2++)
		{
				//x,y
				double WD, CD, d2=double(c2)/double(N*2-4);

				double mmmm=double(D[hash(T, T_MIN, c2, 0, 0)]);
				double Mmmm=double(D[hash(T, T_MIN, c2, 1, 0)]);
				double MMmm=double(D[hash(T, T_MIN, c2, 2, 0)]);
				double mmMm=double(D[hash(T, T_MIN, c2, 0, 1)]);
				double MmMm=double(D[hash(T, T_MIN, c2, 1, 1)]);
				double MMMm=double(D[hash(T, T_MIN, c2, 2, 1)]);
				double mmMM=double(D[hash(T, T_MIN, c2, 0, 2)]);
				double MmMM=double(D[hash(T, T_MIN, c2, 1, 2)]);
				double MMMM=double(D[hash(T, T_MIN, c2, 2, 2)]);

				double Den=mmmm+Mmmm+MMmm+mmMm+MmMm+MMMm+mmMM+MmMM+MMMM;

				double f=1.-double(F[hash_F(T, T_MIN, c2)]/double(Den*N) )/(2.*d2*(1.-d2) );


				double OX=((Mmmm+MmMm+MmMM)/2.+MMmm+MMMm+MMMM+1)/Den;
				double OY=((mmMm+MmMm+MMMm)/2.+mmMM+MmMM+MMMM+1)/Den;

				OX=d2;
				OY=d2;

				double M=(OX+OY)/2.;


				//std::cerr << d2 << ", " << F[hash_F(T, T_MIN, c2)]/double(Den*N) << ", " << f << std::endl;

				//std::cerr <<  mmmm << ", " << Mmmm << ", " << MMmm << ", " << mmMm << ", " << MmMm << ", " << MMMm << ", " << mmMM << ", " <<  MmMM << ", " << MMMM << ", " << d2 << std::endl;

				if (OX > 0 && OY > 0 && OX < 1 && OY <1 && M!=0.5 && Den > 0)
				{
					Theta+=(1.+f)*(1+f)*d2*(1.-d2)*(mmmm*(0-OX)*(0-OY)+Mmmm*(1.0-OX)*(0-OY)/2.+MMmm*(1.-OX)*(0-OY)
								  +Mmmm*(0.0-OX)*(0-OY)/2.
						   +mmMm*(0-OX)*(1.0-OY)/2.+MmMm*(0.0-OX)*(0.0-OY)/4.+MMMm*(1.-OX)*(0.0-OY)/2.
						   +mmMm*(0-OX)*(0.0-OY)/2.+MmMm*(0.0-OX)*(1.0-OY)/4.+MMMm*(1.-OX)*(1.0-OY)/2.
						   			   +MmMm*(1.0-OX)*(0.0-OY)/4.
						   			   +MmMm*(1.0-OX)*(1.0-OY)/4.
						   +mmMM*(0-OX)*(1.-OY)+MmMM*(0.0-OX)*(1.-OY)/2.+MMMM*(1.-OX)*(1.-OY) 
						   		       +MmMM*(1.0-OX)*(1.-OY)/2. )/sqrt(OX*(1-OX)*(1-OY)*OY )/2.;
					Theta_w+=Den*(1.+f)*d2*(1.-d2);

					CD=1;//pow(OX*(1-OX), 1);
					WD=4*d2*(1.-d2)*(d2*(1.-d2)+f*(4*pow(d2,2)-4*d2+1)-pow(f,2)*d2*(1.-d2) );
					k+=CD*WD*( (1./(1.-d2)+1./d2-3) )*Den;
					f_p+=f*CD*WD*Den;
					f_w+=Den*WD;
					k_w+=Den*WD;

					f_X+=CD*WD*(mmmm*(0-OX)*(0-OX)+Mmmm*(0.-OX)*(1.-OX)+MMmm*(1.-OX)*(1.-OX)
						   +mmMm*(0-OX)*(0-OX)+MmMm*(0.-OX)*(1.-OX)+MMMm*(1.-OX)*(1.-OX)
						   +mmMM*(0-OX)*(0-OX)+MmMM*(0.-OX)*(1.-OX)+MMMM*(1.-OX)*(1.-OX) )/(OX*(1-OX) );
					f_X_w+=Den*WD;
					

					f_Y+=CD*WD*(mmmm*(0-OY)*(0-OY)+Mmmm*(0-OY)*(0-OY)+MMmm*(0-OY)*(0-OY)
						   +mmMm*(0.-OY)*(1.-OY)+MmMm*(0.-OY)*(1.-OY)+MMMm*(0.-OY)*(1.-OY)
						   +mmMM*(1.-OY)*(1.-OY)+MmMM*(1.-OY)*(1.-OY)+MMMM*(1.-OY)*(1.-OY) )/(OY*(1-OY ) );
					f_Y_w+=Den*WD;


					gamma_XY+=(pow(d2,3)*(1-d2)-pow(1.-d2,3)*d2 )*(mmmm*(0-OX)*(0-OX)*(0-OY)+Mmmm*(0.-OX)*(1.-OX)*(0-OY)+MMmm*(1.-OX)*(1.-OX)*(0-OY)
							+(mmMm*(0-OX)*(0-OX)*(1.-OY)+MmMm*(1.-OX)*(0.-OX)*(1.-OY)+MMMm*(1.-OX)*(1.-OX)*(1.-OY) )/2.
							+(mmMm*(0-OX)*(0-OX)*(0.-OY)+MmMm*(1.-OX)*(0.-OX)*(0.-OY)+MMMm*(1.-OX)*(1.-OX)*(0.-OY) )/2.
							+mmMM*(0-OX)*(0-OX)*(1.-OY)+MmMM*(1.-OX)*(0.-OX)*(1.-OY)+MMMM*(1.-OX)*(1.-OX)*(1.-OY) )/(M*(1-M)*(1-2*M)); 
					gamma_XY_w+=Den*(pow(d2,3)*(1-d2)-pow(1.-d2,3)*d2 )*f;

					gamma_YX+=(pow(d2,3)*(1-d2)-pow(1.-d2,3)*d2 )*(mmmm*(0-OY)*(0-OY)*(0-OX)+Mmmm*(0.-OY)*(0.-OY)*(0.0-OX)/2.+MMmm*(0.-OY)*(0.-OY)*(1.-OX)
									    +Mmmm*(0.-OY)*(0.-OY)*(1.0-OX)/2.
							+mmMm*(1.-OY)*(0.-OY)*(0.-OX)+MmMm*(1.-OY)*(0.-OY)*(0.0-OX)/2.+MMMm*(1.-OY)*(0.-OY)*(1.-OX)
									   	     +MmMm*(1.-OY)*(0.-OY)*(1.0-OX)/2.
							+mmMM*(1.-OY)*(1.-OY)*(0.-OX)+MmMM*(1.-OY)*(1.-OY)*(0.0-OX)/2.+MMMM*(1.-OY)*(1.-OY)*(1.-OX) 
									    	     +MmMM*(1.-OY)*(1.-OY)*(1.0-OX)/2.)/(M*(1-M)*(1-2*M) );
					gamma_YX_w+=Den*(pow(d2,3)*(1-d2)-pow(1.-d2,3)*d2 )*f;

					mu=(mmmm*(0-OY)*(0-OY)*(0-OX)*(0-OX)+Mmmm*(0.-OY)*(0.-OY)*(1.-OX)*(0.0-OX)+MMmm*(0.-OY)*(0.-OY)*(1.-OX)*(1.-OX)
							+mmMm*(1.-OY)*(0.-OY)*(0.-OX)*(0.-OX)+MmMm*(0.-OY)*(1.-OY)*(0.-OX)*(1.-OX)+MMMm*(0.-OY)*(1.-OY)*(1.-OX)*(1.-OX)
							+mmMM*(1.-OY)*(1.-OY)*(0.-OX)*(0.-OX)+MmMM*(1.-OY)*(1.-OY)*(1.-OX)*(0.-OX)+MMMM*(1.-OY)*(1.-OY)*(1.-OX)*(1.-OX) )/Den;

					double k1x=pow(OX*(1-OX), 2);
					double k2x=OX*(1-OX)*(3*OX*OX-3*OX+1);
					double k1y=pow(OY*(1-OY), 2);
					double k2y=OY*(1-OY)*(3*OY*OY-3*OY+1);

					k1=sqrt(k1x*k1y);
					k2=sqrt(k2x*k2y);
					gsl_matrix_set (X, c2, 0, k1);
					gsl_matrix_set (X, c2, 1, k2);
      
					gsl_vector_set (Y, c2, mu);
					gsl_vector_set (w, c2, Den);
				}
				if (indX!=-1 or indY!=-1){

					Theta_w=0;
				        Theta=0;
				        f_X=0;
				        f_X_w=0;
				        f_Y=0;
				        f_Y_w=0;
				        gamma_XY=0;
				        gamma_XY_w=0;
				        gamma_YX=0;
			        	gamma_YX_w=0;
				}
		}
		gsl_multifit_wlinear (X, w, Y, c, cov, &chisq, work);
			
#define beta(i) (gsl_vector_get(c,(i)))
		if(indX!=-1 or indY!=-1) return 0;

		buffer_rel[z].b_=2*Theta/Theta_w;
		k=k/k_w;
		f_p=f_p/f_w;
		//std::cerr << (f_p*(f_p-f_X/f_X_w-f_Y/f_Y_w) << beta(1) << beta(0) << ", " << k << std::endl;
		buffer_rel[z].d_=(f_p*(f_p-f_X/f_X_w-f_Y/f_Y_w)+k*beta(1)+beta(0) )/(1.-f_p-f_p*f_p+f_p*k);
		buffer_rel[z].bd_=(gamma_XY+gamma_YX)/(gamma_XY_w+gamma_YX_w);
		size_t x, y;
		T.get_xy(x,y);
		buffer_rel[z].set_X_name(x);
		buffer_rel[z].set_Y_name(y);

		if (T.get_k()!=T.size() ) ++T;

		gsl_vector_free(Y);
		gsl_vector_free(w);
		gsl_vector_free(c);

		gsl_matrix_free(X);
		gsl_matrix_free(XINV);
		gsl_matrix_free(cov);
		
		gsl_multifit_linear_free(work);
	}
	for (size_t z=0; z<J; ++z)
		rel_file.write(buffer_rel[z]);
		J=std::min(BUFFER_SIZE2, int(T.size()-T.get_k()) );
	}
	delete [] P1;
	delete [] D;
	rel_file.close();
	//close();
}
