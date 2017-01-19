#include "circular_list.h"

class Slow_index
{
	private:
	size_t _1, _2, _3, _4, _5;
	public:
	Slow_index(const size_t &_N, const size_t &_dist)
	{	
		_1=4*_dist*_N*_N;
		_2=4*_dist*_N; 
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
		std::cerr << 1;
		return 1;
	} else  {
		std::cerr << 0;
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
		*sum+=( (*this_place & mask)==(*(++this_place) & mask) );
		this_place++;
	}
}

/*const size_t & fast_index(const uint &d, const uint32_t &i, const uint32_t &j)
{
	return _index[d][i][j];
}*/

int main (int argc, char **argv){

Circular_list <uint32_t> buffer(0);
Circular_list <uint32_t>::iterator bit=buffer.begin(), tit, end;

for (int x=1; x<20; x++)
{
	bit=buffer.insert(bit, x);
}

bit=buffer.begin();
end=bit;
bit++;
	

while (bit!=end){
	tit=bit.next();
	while (tit!=bit)
	{
		std::cout << *bit << ", " << *tit << std::endl;
		tit++;
	}
	std::cout << std::endl;
	bit++;
}
}	
