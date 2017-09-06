#include "circular_list.h"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cstring>
#include <sstream>
#include <tuple>
#include <map>

inline int get(const uint32_t *place, const uint32_t &x, const uint32_t &bit)
{
	if ((*(place+x) & Mask[bit])!=0 ){
		return 1;
	} else  {
		return 0;

	}
}

uint16_t up32(uint16_t value)
{
	if (in%2==0) return value;

	unsigned pos = 0;
	while (!(value & 0x8000))
	{
		value <<= 1;
		++pos;
	}
	return 0x8000 >> (pos-1);
}

permute()
{
	
}

int main (int argc, char **argv){

	uint16_t N=atoi(argv[1]);
	uint32_t BLOCK_SIZE=2*N;
	uint32_t OUT_SIZE=2*up32(N);
	
	uint32_t *buffer=new uint32_t [OUT_SIZE];

	while(in->read((char *)(bit), BLOCK_SIZE*sizeof(uint32_t) ) )
	{
		permute(bit, buffer, BLOCK_SIZE);
		mem.write((char *)(buffer), BLOCK_SIZE*sizeof(uint32_t) )
	}
	str_mem=mem.str();

}	
