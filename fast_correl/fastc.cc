#include "circular_list.h"

class Slow_index
{
	private:
	size_t _1, _2, _3, _4, _5;
	public:
	Slow_index(const size_t &_N, const size_t &_dist, const size_t &_xmax)
	{	
		_1=4*_dist*_N*_xmax;
		_2=4*_dist*_xmax; 
		_3=4*_dist;
		_4=4; 
		_5=2;
	}

	//d1(N), d2(N), x(N), dist(dist), i(2), j(2) 
	//slow_index(c1, c2, x, y, 1, 1));
	//size_t&, size_t&, size_t&, size_t&, int, int

	const size_t 	
	operator () (const uint32_t &d1, const uint32_t &d2, const uint32_t &x, const uint32_t &dist, const uint32_t &i, const uint32_t &j)	{return d1*_1+d2*_2+x*_3+dist*_4+i*_5+j;};
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

int DIST=atoi(argv[1]);
uint32_t N=atoi(argv[2]);
uint32_t XMAX=atoi(argv[3]);

uint32_t SIZE=XMAX*N*N*2*2*DIST*2*2;
uint32_t BLOCK_SIZE=N;
uint32_t *D=new uint32_t [SIZE];
//Sbuffer=[];

//TO DO : Circular buffer
//

Circular_list <uint32_t *> buffer(new uint32_t [N]);
Circular_list <uint32_t *>::iterator bit=buffer.begin(), tit, end;

for (int x=0; x<DIST; x++)
{
	bit=buffer.insert(bit, new uint32_t [N] );
}

Slow_index slow_index(N, DIST, XMAX);

bit=buffer.begin();
std::istream *in=&std::cin;
	
//std::cerr << "AT" << *bit << std::endl;
uint32_t readed=0;

while(in->read((char *)(*bit), BLOCK_SIZE*sizeof(uint32_t) ) )
{
	//std::cerr << "IN" << *bit << std::endl;
	for (size_t x=0; x<XMAX; x+=2)
	{
		int y=0;
		uint32_t d1, d2;
		uint32_t D_I, D_O, DIST_O=(DIST>>5);
		tit=bit;
		end=bit-DIST_O-1;
		//std::cerr << "OD: " << DIST_O << "<<" << DIST << std::endl;
		if (readed>=DIST_O){
			while (tit!=end)
			{
				if (y==0) {
				for (int d=1; d< std::min(32, DIST); d++){
					for (int k=0; k<32-d; k++){
						int i=(get(*bit, x, k)!=get(*bit, x+1, k) );
						int j=(get(*tit, x, k+d)!=get(*tit, x+1, k+d) );
						d1=0;
						/*
						for (int y=0; y<N; y+=2){
							std::cerr << get(*bit, y, k); 
							std::cerr << get(*bit, y+1, k); 
							std::cerr << ","; 
						}
						std::cerr << "|";
						for (int y=0; y<N; y+=2){
							std::cerr << get(*tit, y, k+d); 
							std::cerr << get(*tit, y+1, k+d); 
							std::cerr << ","; 
						}
						std::cerr << "::";
						*/
						fast_het_sumN(*bit, N, k, &d1);
						d1-=i;
						d2=0;
						fast_het_sumN(*tit, N, k+d, &d2);
						d2-=j;
						//std::cerr << d1 << ", " << d2 << ", " << x << ", (dist)" << d << ", " << i << ", " << j << std::endl;
						(*(D+slow_index(d1, d2, x, d, i, j) ) )++;
					}
				}
				}
				--tit;
				++y;
				for (int d=0; d<std::min(32, DIST-32); d++){
					for (int k=0; k<32; k++){
						int i=(get(*bit, x, k)!=get(*bit, x+1, k) );
						int j=(get(*tit, x, d)!=get(*tit, x+1, d) );
						d1=0;
						/*
						for (int y=0; y<N; y+=2){
							std::cerr << get(*bit, y, k); 
							std::cerr << get(*bit, y+1, k); 
							std::cerr << ","; 
						}
						std::cerr << "|";
						for (int y=0; y<N; y+=2){
							std::cerr << get(*tit, y, k+d); 
							std::cerr << get(*tit, y+1, k+d); 
							std::cerr << ","; 
						}
						std::cerr << "::";
						*/
						fast_het_sumN(*bit, N, k, &d1);
						d1-=i;
						d2=0;
						fast_het_sumN(*tit, N, d, &d2);
						d2-=j;
						//std::cerr << d1 << ", " << d2 << ", " << x << ", (dist)" << y*32+d << ", " << i << ", " << j << std::endl;
						if(y*32-d+k<DIST)
							(*(D+slow_index(d1, d2, x, y*32-d+k, i, j) ) )++;
					}
				}
			}
		}
	}
	++bit;
	++readed;
}

std::cout << "x dist c1 c2 Den rho1 OX OY rho2\n";

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
				double HH=double(*(D+slow_index(c1, c2, x, y, 1, 1)));
				double Hh=double(*(D+slow_index(c1, c2, x, y, 1, 0)));
				double hH=double(*(D+slow_index(c1, c2, x, y, 0, 1)));
				double hh=double(*(D+slow_index(c1, c2, x, y, 0, 0)));
				double Den=double(HH+Hh+hH+hh);
				if (Den>2)
				{
					double OX=double( HH+Hh )/double(Den);
					double OY=double( HH+hH )/double(Den);
					double Num=hh*double(0-d1)*double(0-d2)+Hh*double(1.0-d1)*double(0-d2)+hH*double(0-d1)*double(1.-d2)+HH*double(1.-d1)*double(1.-d2);
					double Num2=hh*double(0-OX)*double(0-OY)+Hh*double(1.0-OX)*double(0-OY)+hH*double(0-OX)*double(1.-OY)+HH*double(1.-OX)*double(1.-OY);
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
