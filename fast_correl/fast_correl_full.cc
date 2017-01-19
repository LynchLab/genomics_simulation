#include "circular_list.h"

class Slow_index
{
	private:
	size_t _1, _2, _3, _4;
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
	operator () (const uint32_t &d, const uint32_t &x, const uint32_t &y, const uint32_t &i, const uint32_t &j)	{return d*_1+x*_2+y*_3+i*_4+j;};
};

static uint32_t Mask[32]={0x00000001, 0x00000002, 0x00000004, 0x00000008,
			0x00000010, 0x00000020, 0x00000040, 0x00000080,
			0x00000100, 0x00000200, 0x00000400, 0x00000800,
			0x00001000, 0x00002000, 0x00004000, 0x00008000,
			0x00010000, 0x00020000, 0x00040000, 0x00080000,
			0x00100000, 0x00200000, 0x00400000, 0x00800000,
			0x01000000, 0x02000000, 0x04000000, 0x08000000,
			0x10000000, 0x20000000, 0x40000000, 0x80000000};

inline void fast_bit_sumN(const uint32_t *place, const uint32_t &N, const uint32_t &bit, uint32_t *sum)
{
	const uint32_t *end=place+N;
	const uint32_t *this_place=place;
	uint32_t mask=Mask[bit];
	while (this_place!=end) {
		*sum+=*this_place & mask;
		this_place++;
	}
}

inline int get(const uint32_t *place, const uint32_t &x, const uint32_t &bit)
{
	//nstd::cerr << "(" << x << ")" << std::endl;
	if ((*(place+x) & Mask[bit])!=0 ){
		return 1;
	} else  {
		return 0;
	}
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

uint32_t N=atoi(argv[1]);

uint32_t SIZE=XMAX*N*N*2*2*DIST*2*2;
uint32_t BLOCK_SIZE=N;
uint32_t *D=new uint32_t [SIZE];

uint32_t * bit=new uint32_t [N];

for (int x=0; x<DIST; x++)
{
	bit=buffer.insert(bit, new uint32_t [N] );
}

Slow_index slow_index(N, DIST, XMAX);

bit=buffer.begin();
std::istream *in=&std::cin;
	
uint32_t readed=0;

while(in->read((char *)(*bit), BLOCK_SIZE*sizeof(uint32_t) ) )
{
	for (size_t x=0; x<N; x+=2)
	{
		for (size_t y=x+1; y<N; y+=2)
			int i=(get(*bit, x, k)+get(*bit, x+1, k) );
			int j=(get(*bit, y, k)+get(*bit, y+1, k) );
			int d=0;
			fast_bit_sumN(*bit, N, k, &d);
			d-=(j+i);
			(*(D+slow_index(d, x, y, i, j) ) )++;
	}
}

std::cout << "x y f_X f_Y Theta gamma_XY gamma_YX Delta delta\n";

for (size_t x=0; x<XMAX; x++)
{
	for (size_t y=1; y<DIST; y++)
	{
		for (size_t c1=1; c1<N-1; c1++)
		{
			for (size_t c2=1; c2<N-1; c2++)
			{
				double d1=double(c1)/double(N/2-2);
				double d2=double(c2)/double(N/2-2);

				double mmmm=double(*(D+slow_index(c1, c2, x, y, 0, 0)));
				double Mmmm=double(*(D+slow_index(c1, c2, x, y, 1, 0)));
				double MMmm=double(*(D+slow_index(c1, c2, x, y, 2, 0)));
				double mmMm=double(*(D+slow_index(c1, c2, x, y, 0, 1)));
				double MmMm=double(*(D+slow_index(c1, c2, x, y, 1, 1)));
				double MMMm=double(*(D+slow_index(c1, c2, x, y, 2, 1)));
				double mmMM=double(*(D+slow_index(c1, c2, x, y, 0, 2)));
				double MmMM=double(*(D+slow_index(c1, c2, x, y, 1, 2)));
				double MMMM=double(*(D+slow_index(c1, c2, x, y, 2, 2)));

				if (Den>10)
				{
					double OX=;
					double OY=;
					double Theta=mmmm*double(0-OX)*double(0-OY)+Mmmm*double(0.5-OX)*double(0-OY)+MMmm*double(1.-OX)*double(0-OY)
						   +mmMm*double(0-OX)*double(0.5-OY)+MmMm*double(0.5-OX)*double(0.5-OY)+MMMm*double(1.-OX)*double(0.5-OY)
						   +mmMM*double(0-OX)*double(1.-OY)+MmMM*double(0.5-OX)*double(1.-OY)+MMMM*double(1.-OX)*double(1.-OY);

					double f_X=mmmm*double(0-OX)*double(0-OX)+Mmmm*double(0.5-OX)*double(0.5-OX)+MMmm*double(1.-OX)*double(1.-OX)
						   +mmMm*double(0-OX)*double(0-OX)+MmMm*double(0.5-OX)*double(0.5-OX)+MMMm*double(1.-OX)*double(1.-OX)
						   +mmMM*double(0-OX)*double(0-OX)+MmMM*double(0.5-OX)*double(0.5-OX)+MMMM*double(1.-OX)*double(1.-OX);

					double f_Y=mmmm*double(0-OY)*double(0-OY)+Mmmm*double(0-OY)*double(0-OY)+MMmm*double(0-OY)*double(0-OY)
						   +mmMm*double(0.5-OY)*double(0.5-OY)+MmMm*double(0.5-OY)*double(0.5-OY)+MMMm*double(0.5-OY)*double(0.5-OY)
						   +mmMM*double(1.-OY)*double(1.-OY)+MmMM*double(1.-OY)*double(1.-OY)+MMMM*double(1.-OY)*double(1.-OY);

					double f_Y=mmmm*double(0-OY)*double(0-OY)+Mmmm*double(0-OY)*double(0-OY)+MMmm*double(0-OY)*double(0-OY)
						   +mmMm*double(0.5-OY)*double(0.5-OY)+MmMm*double(0.5-OY)*double(0.5-OY)+MMMm*double(0.5-OY)*double(0.5-OY)
						   +mmMM*double(1.-OY)*double(1.-OY)+MmMM*double(1.-OY)*double(1.-OY)+MMMM*double(1.-OY)*double(1.-OY);

					if (OX>0 && OX<1 && OY>0 && OY<1) 
						std::cout << x << " " << y << " " << c1 << " " << c2 << " " << Den << " " << Num/Den/sqrt((1-d1)*d1*d2*(1-d2) ) << " " << OX << " " << OY <<" " << Num2/Den/sqrt((1-OX)*OX*OY*(1-OY) ) << std::endl;
					else
						std::cout << x << " " << y << " " << c1 << " " << c2 << " " << Den << " " << Num/Den/sqrt((1-d1)*d1*d2*(1-d2) ) << " " << OX << " " << OY <<" NaN\n";
				}
			}
		}
	}
}	
}	
