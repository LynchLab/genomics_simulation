#include "circular_list.h"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cstring>

class Slow_index
{
	private:
	size_t _1, _2, _3, _4, _5;
	public:
	Slow_index(const size_t &_N)
	{	
		_1=9*_N*_N;
		_2=9*_N; 
		_3=9;
		_4=3;
	}

	//d1(N), x(N), y(N), i(3), j(3) 

	const size_t 	
	operator () (const uint32_t &c, const uint32_t &x, const uint32_t &y, const uint32_t &i, const uint32_t &j)	{return c*_1+x*_2+y*_3+i*_4+j;};
};

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
	//nstd::cerr << "(" << x << ")" << std::endl;
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
		//if((*this_place & mask)==0) std::cerr << 0; 
		//else std::cerr << 1;
		//if((*(this_place+1) & mask)==0) std::cerr << 0; 
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
//	for (size_t x=0; x<N-2; ++x) (*sum)+=get(place, x, bit);
}

inline void fast_het_sumN(const uint32_t *place, const uint32_t &N, const uint32_t &bit, uint32_t *sum)
{
	const uint32_t *end=place+N;
	const uint32_t *this_place=place;
	uint32_t mask=Mask[bit];
	while (this_place!=end) {
		//if((*this_place & mask)==0) std::cerr << 0; 
		//else std::cerr << 1;
		//if((*(this_place+1) & mask)==0) std::cerr << 0; 
		//else std::cerr << 1;
		*sum+=( (*this_place & mask)!=(*(++this_place) & mask) );
		this_place++;
	}
}

/*const size_t & fast_index(const uint &d, const uint32_t &i, const uint32_t &j)
{
	return _index[d][i][j];
}*/

int main (int argc, char **argv){

uint64_t N=atoi(argv[1]);

uint64_t SIZE=N; 
uint32_t BLOCK_SIZE=2*N;
//HUGE!!!
uint16_t *D=new uint16_t [SIZE];
double *F=new double [SIZE];
//SMALL
uint32_t * bit=new uint32_t [N*2];


std::istream *in=&std::cin;
	
uint32_t readed=0;

while(in->read((char *)(bit), BLOCK_SIZE*sizeof(uint32_t) ) )
{
	uint32_t d1[32]={0};
	fast_bit_sumN(bit, N*2, d1);

	for (size_t k=0; k<32; k++)
	{
		D[d1[k]]++;
		F[d1[k]]+=H[k];
	}
}

std::cerr << "done reading." << std::endl;

/* rho_xy|z=rho_xy-rho_xz*rho_yz/( sqrt(1-rho_xz^2)*rho(1-rho_yz^2) */

for (size_t x=0; x<N; x++)
{
	for (size_t Y=x+1; Y<N; Y++)

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
	}
}	
	std::cout << "@END_TABE\n";
}	
