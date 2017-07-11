#include "circular_list.h"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cstring>
#include <sstream>
#include "map-file.h"
//#include "
#include "interface.h"

class Slow_index
{
	private:
	size_t _1, _2;
	public:
	Slow_index(const size_t &_N)
	{	
		_1=9;
		_2=3; 
	}

	//d1(N), x(N), y(N), i(3), j(3) 
	const size_t 	
	operator () (const uint32_t &c, const uint32_t &i, const uint32_t &j)	{return c*_1+i*_2+j;};
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

//TODO make this fast. Right now it is SLOW SLOW SLOW!!!
inline void fast_bit_sumN(const uint32_t *place, const uint32_t &N, uint32_t *sum)
{
	const uint32_t *end=place+N;
	const uint32_t *this_place=place;

/*	uint32_t SIZE=N/8+N%8!=0;
	char *bits=new char SIZE;

	for (int k=0; k<32; k++){
		bits=
	}*/
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

if (argc!=2) std::cerr << "usage: " << argv[0] << " <NUMBER OF SAMPLES>\n";

uint64_t N=atoi(argv[1]);

uint64_t SIZE=9*(N+1); 
uint32_t BLOCK_SIZE=2*N;
uint16_t *D=new uint16_t [SIZE];
uint32_t * bit=new uint32_t [N*2];

Slow_index slow_index(N);

//std::cerr << ": " << N << "--" << SIZE << " " << slow_index(N, 3-1, 3-1) << std::endl;

std::istream *in=&std::cin;
	
uint32_t readed=0;

std::cout << "@NAME:RELATEDNESS	VERSION:0.4.27-2016-12-20	FORMAT:TEXT	CONCATENATED \n@SAMPLE_X	SAMPLE_Y	f_X	f_X_ll	f_Y	f_Y_ll	θ_XY	θ_ll	γ_XY	γ_XY_ll	γ_YX	γ_YX_ll	δ	δ_ll	Δ	Δ_ll	null_ll	fit\n";

double f_X=0, f_X_w=0, f_Y=0, f_Y_w=0, Theta=0, Theta_w=0, gamma_XY=0, gamma_XY_w=0, gamma_YX=0, gamma_YX_w=0, Z=0, Z_w=0;
double mu, k1, k2;

gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (pow(N*2-1, 1), 2);
gsl_matrix *Xm = gsl_matrix_alloc (pow(N*2-1, 1), 2);
gsl_matrix *XmINV = gsl_matrix_alloc (pow(N*2-1, 1), 2);
gsl_vector *yv = gsl_vector_alloc (pow(N*2-1, 1) );
gsl_vector *wv = gsl_vector_alloc (pow(N*2-1, 1) );
gsl_vector *cv = gsl_vector_alloc (2);
gsl_matrix *cov = gsl_matrix_alloc (2, 2);
double chisq;

//while(in->read((char *)(bit), BLOCK_SIZE*sizeof(uint32_t) ) )
//TODO may want to move fast_bit_sumN here

std::stringstream buffer;

while(in->read((char *)(bit), BLOCK_SIZE*sizeof(uint32_t) ) )
{
	buffer.write( (char *)(bit), BLOCK_SIZE*sizeof(uint32_t) );
}
uint32_t BLOCKS=buffer.str().size()/(BLOCK_SIZE);
std::stringstream copy1( buffer.str() );
in=&copy1;

std::cerr << BLOCKS << "\n";

bool *flip_array=new bool [BLOCKS*32];
uint32_t *d1_array=new uint32_t [BLOCKS*32];

bool *flip_array_ptr=flip_array;
uint32_t *d1_array_ptr=d1_array;

//TODO Move this?
memset((char *)(flip_array), 0, BLOCKS*sizeof(bool)*32 );
memset((char *)(d1_array), 0, BLOCKS*sizeof(uint32_t)*32 );
while(in->read((char *)(bit), BLOCK_SIZE*sizeof(uint32_t) ) )
{
	//std::cerr << "settting " << int(d1_array_ptr-d1_array) << std::endl;
	uint32_t *d1=d1_array_ptr+=32;
	bool *flip=flip_array_ptr+=32;
{
	fast_bit_sumN(bit, N*2, d1);
	for (size_t k=0; k<32; k++)
	{
		if (d1[k]>N) 
		{
			d1[k]=2*N-d1[k];
			flip[k]=true;
		} 
	}
}


}
//std::cerr << "DONE READING IN\n";

for (size_t x=0; x<N*2; x+=2)
{
	for (size_t y=x+2; y<N*2; y+=2)
	{
//obj=buffer.string();
std::stringstream copy2( buffer.str() );
in=&copy2;
memset((char *)(D), 0, SIZE*sizeof(uint16_t) );
d1_array_ptr=d1_array;
flip_array_ptr=flip_array;
while(in->read((char *)(bit), BLOCK_SIZE*sizeof(uint32_t) ) )
{
	//std::cerr << int(d1_array_ptr-d1_array) << std::endl;
	uint32_t *d1=d1_array_ptr+=32;
	bool *flip=flip_array_ptr+=32;
	for (size_t k=0; k<32; k++)
	{
		uint32_t mask=Mask[k];
		int i=get(bit, x, k)+get(bit, x+1, k);
		int j=get(bit, y, k)+get(bit, y+1, k);
		if(flip[k])
		{
			i=2-i;
			j=2-j;
		}
		//std::cerr << d1[k] << ", " << i << ", " << j << ", " << std::endl;
		(*(D+slow_index(d1[k]-(i+j), i, j) ) )++;
	}
}

//std::cerr << "DONE READING IN2\n";

	size_t X=x/2;
	size_t Y=y/2;
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
		for (size_t c2=1; c2<(N-1); c2++)
		{
				double d2=double(c2)/double(N*2-4);

				double mmmm=double(*(D+slow_index(c2, 0, 0)));
				double Mmmm=double(*(D+slow_index(c2, 1, 0)));
				double MMmm=double(*(D+slow_index(c2, 2, 0)));
				double mmMm=double(*(D+slow_index(c2, 0, 1)));
				double MmMm=double(*(D+slow_index(c2, 1, 1)));
				double MMMm=double(*(D+slow_index(c2, 2, 1)));
				double mmMM=double(*(D+slow_index(c2, 0, 2)));
				double MmMM=double(*(D+slow_index(c2, 1, 2)));
				double MMMM=double(*(D+slow_index(c2, 2, 2)));

				double Den=mmmm+Mmmm+MMmm+mmMm+MmMm+MMMm+mmMM+MmMM+MMMM;

				double OX=((Mmmm+MmMm+MmMM)/2.+MMmm+MMMm+MMMM)/Den;
				double OY=((mmMm+MmMm+MMMm)/2.+mmMM+MmMM+MMMM)/Den;

				if (Den>1 && OX > 0 && OY > 0 && OX < 1 && OY <1 && OX!=0.5 && OY != 0.5)
				{i
					double M=(OX+OY)/2.;

					Theta+=(mmmm*(0-OX)*(0-OY)+Mmmm*(1.0-OX)*(0-OY)/2.+MMmm*(1.-OX)*(0-OY)
								  +Mmmm*(0.0-OX)*(0-OY)/2.
						   +mmMm*(0-OX)*(1.0-OY)/2.+MmMm*(0.0-OX)*(0.0-OY)/4.+MMMm*(1.-OX)*(0.0-OY)/2.
						   +mmMm*(0-OX)*(0.0-OY)/2.+MmMm*(0.0-OX)*(1.0-OY)/4.+MMMm*(1.-OX)*(1.0-OY)/2.
						   			   +MmMm*(1.0-OX)*(0.0-OY)/4.
						   			   +MmMm*(1.0-OX)*(1.0-OY)/4.
						   +mmMM*(0-OX)*(1.-OY)+MmMM*(0.0-OX)*(1.-OY)/2.+MMMM*(1.-OX)*(1.-OY) 
						   		       +MmMM*(1.0-OX)*(1.-OY)/2. )/sqrt(OX*(1-OX)*(1-OY)*OY );
					Theta_w+=Den;

					std::cerr << M << ", " << (mmmm*(0-OX)*(0-OY)+Mmmm*(1.0-OX)*(0-OY)/2.+MMmm*(1.-OX)*(0-OY)
                                                                   +Mmmm*(0.0-OX)*(0-OY)/2.
                                                    +mmMm*(0-OX)*(1.0-OY)/2.+MmMm*(0.0-OX)*(0.0-OY)/4.+MMMm*(1.-OX)*(0.0-OY)/2.
                                                    +mmMm*(0-OX)*(0.0-OY)/2.+MmMm*(0.0-OX)*(1.0-OY)/4.+MMMm*(1.-OX)*(1.0-OY)/2.
                                                                            +MmMm*(1.0-OX)*(0.0-OY)/4.
                                                                            +MmMm*(1.0-OX)*(1.0-OY)/4.
                                                    +mmMM*(0-OX)*(1.-OY)+MmMM*(0.0-OX)*(1.-OY)/2.+MMMM*(1.-OX)*(1.-OY)
                                                                        +MmMM*(1.0-OX)*(1.-OY)/2. )/sqrt(OX*(1-OX)*(1-OY)*OY )/Den;


					f_X+=(mmmm*(0-OX)*(0-OX)+Mmmm*(0.-OX)*(1.-OX)+MMmm*(1.-OX)*(1.-OX)
						   +mmMm*(0-OX)*(0-OX)+MmMm*(0.-OX)*(1.-OX)+MMMm*(1.-OX)*(1.-OX)
						   +mmMM*(0-OX)*(0-OX)+MmMM*(0.-OX)*(1.-OX)+MMMM*(1.-OX)*(1.-OX) )/(OX*(1-OX) );
					f_X_w+=Den;
					

					f_Y+=(mmmm*(0-OY)*(0-OY)+Mmmm*(0-OY)*(0-OY)+MMmm*(0-OY)*(0-OY)
						   +mmMm*(0.-OY)*(1.-OY)+MmMm*(0.-OY)*(1.-OY)+MMMm*(0.-OY)*(1.-OY)
						   +mmMM*(1.-OY)*(1.-OY)+MmMM*(1.-OY)*(1.-OY)+MMMM*(1.-OY)*(1.-OY) )/(OY*(1-OY ) );
					f_Y_w+=Den;
					if (M<0.45){
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
										    	     +MmMM*(1.-OY)*(1.-OY)*(1.0-OX)/2.
														)/(M*(1-M)*(1-2*M) );
						gamma_YX_w+=Den;
					}
					mu=(mmmm*(0-OY)*(0-OY)*(0-OX)*(0-OX)+Mmmm*(0.-OY)*(0.-OY)*(1.-OX)*(0.0-OX)+MMmm*(0.-OY)*(0.-OY)*(1.-OX)*(1.-OX)
							+mmMm*(1.-OY)*(0.-OY)*(0.-OX)*(0.-OX)+MmMm*(0.-OY)*(1.-OY)*(0.-OX)*(1.-OX)+MMMm*(0.-OY)*(1.-OY)*(1.-OX)*(1.-OX)
							+mmMM*(1.-OY)*(1.-OY)*(0.-OX)*(0.-OX)+MmMM*(1.-OY)*(1.-OY)*(1.-OX)*(0.-OX)+MMMM*(1.-OY)*(1.-OY)*(1.-OX)*(1.-OX) )/Den;
					double k1x=pow(OX*(1-OX), 2);
					double k2x=OX*(1-OX)*(3*OX*OX-3*OX+1);
					double k1y=pow(OY*(1-OY), 2);
					double k2y=OY*(1-OY)*(3*OY*OY-3*OY+1);
					k1=sqrt(k1x*k1y);
					k2=sqrt(k2x*k2y);
					gsl_matrix_set (Xm, c2, 0, k1);
					gsl_matrix_set (Xm, c2, 1, k2);
      
					gsl_vector_set (yv, c2, mu);
					gsl_vector_set (wv, c2, Den);
			}
		}
		#define beta(i) (gsl_vector_get(cv,(i)))
		gsl_multifit_wlinear (Xm, wv, yv, cv, cov, &chisq, work);
		std::cout << X << "\t" << Y << "\t" << f_X/f_X_w << "\t" << f_X_w << "\t" << f_Y/f_Y_w << "\t" << f_Y_w << "\t" << Theta/Theta_w << "\t" <<Theta_w <<"\t" << gamma_XY/gamma_XY_w << "\t" << gamma_XY_w << "\t" << gamma_YX/gamma_YX_w << "\t" << gamma_YX_w << "\t" << beta(1) << "\t" << "0" << "\t" << beta(0) << "\t" << std::endl;

	}
}	
}

/* rho_xy|z=rho_xy-rho_xz*rho_yz/( sqrt(1-rho_xz^2)*rho(1-rho_yz^2) */
	std::cout << "@END_TABE\n";
}	
