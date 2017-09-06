#include "circular_list.h"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cstring>
#include <tuple>
#include <map>

uint64_t hash(const uint16_t &d, const uint16_t &x, const uint16_t &i)
{
	uint64_t d64=d, x64=x, i64=i;
	return (d64<<32)+(x64<<16)+i64;
}

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

inline void fast_bit_sumN(const uint32_t *place, const uint32_t &N, uint32_t *sum)
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
}

inline void fast_het_sumN(const uint32_t *place, const uint32_t &N, const uint32_t &bit, uint32_t *sum)
{
	const uint32_t *end=place+N;
	const uint32_t *this_place=place;
	uint32_t mask=Mask[bit];
	while (this_place!=end) {
		*sum+=( (*this_place & mask)!=(*(++this_place) & mask) );
		this_place++;
	}
}

int main (int argc, char **argv){

uint64_t N=atoi(argv[1]);

uint64_t SIZE=3*2*N*N; 
uint32_t BLOCK_SIZE=2*N;

std::cerr << ": " << N << " " << float(SIZE*16)/8000000000. << " GB" << std::endl;

//HUGE!!!
//uint16_t *D=new uint16_t [SIZE];
std::map <uint64_t , uint16_t> D;

//SMALL
uint32_t * bit=new uint32_t [N*2];


std::istream *in=&std::cin;
	
uint32_t readed=0;
uint64_t byteread=0;

while(in->read((char *)(bit), BLOCK_SIZE*sizeof(uint32_t) ) )
{
	//if (byteread%(BLOCK_SIZE*sizeof(uint32_t)*1000 )==0) 
	std::cerr << byteread << std::endl;
	uint32_t d1[32]={0};
	fast_bit_sumN(bit, N*2, d1);
#pragma 
	for (size_t x=0; x<N; x++)
	{
		for (size_t k=0; k<32; k++)
		{
			int i=get(bit, x, k)+get(bit, N+x, k);
			D[hash(d1[k]-i, x, i) ]++;
		}
	}
	byteread+=BLOCK_SIZE*sizeof(uint32_t);
}

std::cerr << "done reading." << std::endl;

std::cout << "@NAME:RELATEDNESS	VERSION:0.4.27-2016-12-20	FORMAT:TEXT	CONCATENATED \n@SAMPLE_X	SAMPLE_Y	f_X	f_X_ll	f_Y	f_Y_ll	θ_XY	θ_ll	γ_XY	γ_XY_ll	γ_YX	γ_YX_ll	δ	δ_ll	Δ	Δ_ll	null_ll	fit\n";

double f_X=0, f_X_w=0, f_Y=0, f_Y_w=0, Theta=0, Theta_w=0, gamma_XY=0, gamma_XY_w=0, gamma_YX=0, gamma_YX_w=0, Z=0, Z_w=0;
double mu, k1, k2;

gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (pow(N*2-1, 1), 2);
gsl_matrix *X = gsl_matrix_alloc (pow(N*2-1, 1), 2);
gsl_matrix *XINV = gsl_matrix_alloc (pow(N*2-1, 1), 2);
gsl_vector *y = gsl_vector_alloc (pow(N*2-1, 1) );
gsl_vector *w = gsl_vector_alloc (pow(N*2-1, 1) );
gsl_vector *c = gsl_vector_alloc (2);
gsl_matrix *cov = gsl_matrix_alloc (2, 2);
double chisq;


for (size_t x=0; x<N; x++)
{
	for (size_t Y=x+1; Y<N; Y++)
	{
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

		//std::cerr <<  x << ", " << Y << std::endl;
		for (size_t c2=1; c2<(N-1)*2; c2++)
		{
				double d2=double(c2)/double(N*2-4);

				double mmmm=double(D[hash(c2, x, 0)]+D[hash(c2, Y, 0)]);
				double Mmmm=double(D[hash(c2, x, 1)]+D[hash(c2, Y, 0)]);
				double MMmm=double(D[hash(c2, x, 2)]+D[hash(c2, Y, 0)]);
				double mmMm=double(D[hash(c2, x, 0)]+D[hash(c2, Y, 1)]);
				double MmMm=double(D[hash(c2, x, 1)]+D[hash(c2, Y, 1)]);
				double MMMm=double(D[hash(c2, x, 2)]+D[hash(c2, Y, 1)]);
				double mmMM=double(D[hash(c2, x, 0)]+D[hash(c2, Y, 2)]);
				double MmMM=double(D[hash(c2, x, 1)]+D[hash(c2, Y, 2)]);
				double MMMM=double(D[hash(c2, x, 2)]+D[hash(c2, Y, 2)]);

				double Den=mmmm+Mmmm+MMmm+mmMm+MmMm+MMMm+mmMM+MmMM+MMMM;

				double OX=((Mmmm+MmMm+MmMM)/2.+MMmm+MMMm+MMMM)/Den;
				double OY=((mmMm+MmMm+MMMm)/2.+mmMM+MmMM+MMMM)/Den;
				double M=(OX+OY)/2.;

				//OX=d2;
				//OY=d2;

				//std::cerr <<  mmmm << ", " << Mmmm << ", " << MMmm << ", " << mmMm << ", " << MmMm << ", " << MMMm << ", " << mmMM << ", " <<  MmMM << ", " << MMMM << ", " << d2 << std::endl;

				if (OX > 0 && OY > 0 && OX < 1 && OY <1 && M!=0.5)
				{
/*					std::cerr << d2 << " "<< OX << " " << OY << " " << (mmmm*(0-OX)*(0-OY)+Mmmm*(1.0-OX)*(0-OY)/2.+MMmm*(1.-OX)*(0-OY)
								  +Mmmm*(0.0-OX)*(0-OY)/2.
						   +mmMm*(0-OX)*(1.0-OY)/2.+MmMm*(0.0-OX)*(0.0-OY)/4.+MMMm*(1.-OX)*(0.0-OY)/2.
						   +mmMm*(0-OX)*(0.0-OY)/2.+MmMm*(0.0-OX)*(1.0-OY)/4.+MMMm*(1.-OX)*(1.0-OY)/2.
						   			   +MmMm*(1.0-OX)*(0.0-OY)/4.
						   			   +MmMm*(1.0-OX)*(1.0-OY)/4.
						   +mmMM*(0-OX)*(1.-OY)+MmMM*(0.0-OX)*(1.-OY)/2.+MMMM*(1.-OX)*(1.-OY) 
						   		       +MmMM*(1.0-OX)*(1.-OY)/2. )/sqrt(OX*(1-OX)*(1-OY)*OY )/Den << " " << Den << std::endl;
*/
					Theta+=(mmmm*(0-OX)*(0-OY)+Mmmm*(1.0-OX)*(0-OY)/2.+MMmm*(1.-OX)*(0-OY)
								  +Mmmm*(0.0-OX)*(0-OY)/2.
						   +mmMm*(0-OX)*(1.0-OY)/2.+MmMm*(0.0-OX)*(0.0-OY)/4.+MMMm*(1.-OX)*(0.0-OY)/2.
						   +mmMm*(0-OX)*(0.0-OY)/2.+MmMm*(0.0-OX)*(1.0-OY)/4.+MMMm*(1.-OX)*(1.0-OY)/2.
						   			   +MmMm*(1.0-OX)*(0.0-OY)/4.
						   			   +MmMm*(1.0-OX)*(1.0-OY)/4.
						   +mmMM*(0-OX)*(1.-OY)+MmMM*(0.0-OX)*(1.-OY)/2.+MMMM*(1.-OX)*(1.-OY) 
						   		       +MmMM*(1.0-OX)*(1.-OY)/2. )/sqrt(OX*(1-OX)*(1-OY)*OY );
					Theta_w+=Den;

					f_X+=(mmmm*(0-OX)*(0-OX)+Mmmm*(0.-OX)*(1.-OX)+MMmm*(1.-OX)*(1.-OX)
						   +mmMm*(0-OX)*(0-OX)+MmMm*(0.-OX)*(1.-OX)+MMMm*(1.-OX)*(1.-OX)
						   +mmMM*(0-OX)*(0-OX)+MmMM*(0.-OX)*(1.-OX)+MMMM*(1.-OX)*(1.-OX) )/(OX*(1-OX) );
					f_X_w+=Den;
					

					f_Y+=(mmmm*(0-OY)*(0-OY)+Mmmm*(0-OY)*(0-OY)+MMmm*(0-OY)*(0-OY)
						   +mmMm*(0.-OY)*(1.-OY)+MmMm*(0.-OY)*(1.-OY)+MMMm*(0.-OY)*(1.-OY)
						   +mmMM*(1.-OY)*(1.-OY)+MmMM*(1.-OY)*(1.-OY)+MMMM*(1.-OY)*(1.-OY) )/(OY*(1-OY ) );
					f_Y_w+=Den;
		/*		}
		}
		for (size_t c2=1; c2<(N-1)*2; c2++)
		{*/
					gamma_XY+=(mmmm*(0-OX)*(0-OX)*(0-OY)+Mmmm*(0.-OX)*(1.-OX)*(0-OY)+MMmm*(1.-OX)*(1.-OX)*(0-OY)
							+(mmMm*(0-OX)*(0-OX)*(1.-OY)+MmMm*(1.-OX)*(0.-OX)*(1.-OY)+MMMm*(1.-OX)*(1.-OX)*(1.-OY) )/2.
							+(mmMm*(0-OX)*(0-OX)*(0.-OY)+MmMm*(1.-OX)*(0.-OX)*(0.-OY)+MMMm*(1.-OX)*(1.-OX)*(0.-OY) )/2.
							+mmMM*(0-OX)*(0-OX)*(1.-OY)+MmMM*(1.-OX)*(0.-OX)*(1.-OY)+MMMM*(1.-OX)*(1.-OX)*(1.-OY) )/(M*(1-M)*(1-2*M)); //p(1-p)(1-2p)
					gamma_XY_w+=Den;

					gamma_YX+=(mmmm*(0-OY)*(0-OY)*(0-OX)+Mmmm*(0.-OY)*(0.-OY)*(0.0-OX)/2.+MMmm*(0.-OY)*(0.-OY)*(1.-OX)
									    +Mmmm*(0.-OY)*(0.-OY)*(1.0-OX)/2.
							+mmMm*(1.-OY)*(0.-OY)*(0.-OX)+MmMm*(1.-OY)*(0.-OY)*(0.0-OX)/2.+MMMm*(1.-OY)*(0.-OY)*(1.-OX)
									   	     +MmMm*(1.-OY)*(0.-OY)*(1.0-OX)/2.
							+mmMM*(1.-OY)*(1.-OY)*(0.-OX)+MmMM*(1.-OY)*(1.-OY)*(0.0-OX)/2.+MMMM*(1.-OY)*(1.-OY)*(1.-OX) 
									    	     +MmMM*(1.-OY)*(1.-OY)*(1.0-OX)/2.)/(M*(1-M)*(1-2*M) );
					gamma_YX_w+=Den;
		/*}
		for (size_t c2=1; c2<(N-1)*2; c2++)
		{*/

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
      
					gsl_vector_set (y, c2, mu);
					gsl_vector_set (w, c2, Den);
			}
		}

/*
		std::cerr << "HI!\n";
		gsl_permutation *p = gsl_permutation_alloc(pow(N-1,2));
		std::cerr << "This is bullshit1!\n";
		int s;
		gsl_linalg_LU_decomp(X, p, &s);
		std::cerr << "This is bullshit2!\n";
		gsl_linalg_LU_solve(X, p, y, c);
		std::cerr << "This is bullshit3!\n";
*/
		/*gsl_blas_dgemm (CblasTrans, CblasNoTrans, CblasNonUnit, *X, *X, CblasNonUnit, *XINV);
		gsl_linalg_cholesky_decomp1(*XINV);
		gsl_linalg_cholesky_invert (*XINV);
		gsl_blas_dgemm (CblasNoTrans, CblasTrans, CblasNonUnit, *XINV, *X, CblasNonUnit, *XR);
		gsl_blas_dgemv (CblasNoTrans, CblasNoTrans, CblasNonUnit, *XR, *y, CblasNonUnit, *c);*
*///		gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
		gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
//		gsl_multifit_linear_free (work);
			
#define beta(i) (gsl_vector_get(c,(i)))

		std::cout << x << "\t" << Y << "\t" << f_X/f_X_w << "\t" << f_X_w << "\t" << f_Y/f_Y_w << "\t" << f_Y_w << "\t" << Theta/Theta_w << "\t" <<Theta_w <<"\t" << gamma_XY/gamma_XY_w << "\t" << gamma_XY_w << "\t" << gamma_YX/gamma_YX_w << "\t" << gamma_YX_w << "\t" << beta(1) << "\t" << "0" << "\t" << beta(0) << "\t" << std::endl;
//			return 0;
	}
}	
	std::cout << "@END_TABE\n";
}	
