#include "circular_list.h"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cstring>

class Slow_index
{
	private:
	size_t _1, _2;
	public:
	Slow_index(const size_t &_N)
	{	
		_1=3*_N;
		_2=3; 
	}

	//d1(N), x(N), y(N), i(3), j(3) 

	const size_t 	
	operator () (const uint32_t &c, const uint32_t &x, const uint32_t &i,)	{return c*_1+x*_2+i;};
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

uint64_t SIZE=9*2*N*N*N; 
uint32_t BLOCK_SIZE=2*N;
//HUGE!!!
uint16_t *D=new uint16_t [SIZE];
//SMALL
uint32_t * bit=new uint32_t [N*2];


Slow_index slow_index(N);


std::istream *in=&std::cin;
	
uint32_t readed=0;

while(in->read((char *)(bit), BLOCK_SIZE*sizeof(uint32_t) ) )
{
	uint32_t d1[32]={0};
	fast_bit_sumN(bit, N*2, d1);
#pragma 
	for (size_t x=0; x<N; x++)
	{
		for (size_t k=0; k<32; k++)
		{

			int i=get(bit, x, k)+get(bit, N+x, k);
			(*(D+slow_index(d1[k], x, i) ) )++;
		}
	}
}

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


for (size_t n=1; n<N; n++)
{
	size_t max=nChose(N,n);
	for (size_t k=0; k<N-n; k++) {
	x=?;
	for (size_t c2=1; c2<(N-n)*2; c2++)
	{
			double d2=double(c2)/(double(N-n)*2);

			double Den=?;

			double O?=()/Den;
			double M=(O?)/?.;

			if (OX > 0 && OY > 0 && OX < 1 && OY <1 && M!=0.5)
			{
				

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
	}
			
	std::cout << "@END_TABE\n";
}	
