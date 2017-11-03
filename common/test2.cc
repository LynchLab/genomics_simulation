#include <stdio.h>
#include <cstdint>
#include <iostream>       // std::cout
#include <string>         // std::string
#include <bitset>         // std::bitset



#include <assert.h>
#include <emmintrin.h>
#include <stdint.h>

//#define II  ( i ^ 7 )
#define II  i 

void
sse_trans(uint8_t const *inp, uint8_t *out, int nrows, int ncols)
{
#   define INP(x,y) inp[(x)*ncols/8 + (y)/8]
#   define OUT(x,y) out[(y)*nrows/8 + (x)/8]
    int rr, cc, i, h;
    union { __m128i x; uint8_t b[16]; } tmp;
    assert(nrows % 8 == 0 && ncols % 8 == 0);

    // Do the main body in 16x8 blocks:
    for (rr = 0; rr <= nrows - 16; rr += 16) {
        for (cc = 0; cc < ncols; cc += 8) {
            for (i = 0; i < 16; ++i)
                tmp.b[i] = INP(rr + II, cc);
            for (i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1))
                *(uint16_t*)&OUT(rr,cc+II)= _mm_movemask_epi8(tmp.x);
        }
    }
    if (rr == nrows) return;

    // The remainder is a block of 8x(16n+8) bits (n may be 0).
    //  Do a PAIR of 8x8 blocks in each step:
    for (cc = 0; cc <= ncols - 16; cc += 16) {
        for (i = 0; i < 8; ++i) {
            tmp.b[i] = h = *(uint16_t const*)&INP(rr + II, cc);
            tmp.b[i + 8] = h >> 8;
        }
        for (i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1)) {
            OUT(rr, cc + II) = h = _mm_movemask_epi8(tmp.x);
            OUT(rr, cc + II + 8) = h >> 8;
        }
    }
    if (cc == ncols) return;

    //  Do the remaining 8x8 block:
    for (i = 0; i < 8; ++i)
        tmp.b[i] = INP(rr + II, cc);
    for (i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1))
        OUT(rr, cc + II) = _mm_movemask_epi8(tmp.x);
}

//Copied from Hacker's Delight, by 

void transpose32a(unsigned long a[32]) 
{
	int j, k;
	unsigned long m, t;
	for (j = 16, m = 0x0000FFFF; j; j >>= 1, m ^= m << j) 
	{
		for (k = 0; k < 32; k = ((k | j) + 1) & ~j) 
		{
			t = (a[k] ^ (a[k|j] >> j)) & m;
			a[k] ^= t;
			a[k|j] ^= (t << j);
		}
	}
}

void transpose8(uint32_t A[], const int &n)
{
	
} 

void slow_trans(uint8_t const *inp, uint8_t *out, int nrows, int ncols)
{
	uint8_t bmask[8]={1,2,4,8,16,32,64,128};
	for (int x=0; x<nrows/8; x++) 
	{
		for (int y=0; y<ncols; y++) 
		{
			out[y]=out[y] | inp[x] & bmask[y];
		}
	}
} 



int main (int argc, char *argv[] )
{
	int N=16;
	uint8_t A[16];
	uint8_t B[16]= {0};
	A[0]=0xf0;
	A[1]=0x00;
	A[2]=0x00;
	A[3]=0x00;
	A[4]=0x00;
	A[5]=0x00;
	A[6]=0x00;
	A[7]=0x00;
	A[8]=0x0f;
	A[9]=0x0f;
	A[10]=0x0f;
	A[11]=0x0f;
	A[12]=0xf0;
	A[13]=0xf0;
	A[14]=0xf0;
	A[15]=0xf0;

	std::bitset<8> bits;     // mybits: 0000
	for (int x=0; x<N; x++)
	{
		bits=A[x];
		std::cout << bits.to_string() << std::endl;
	}

	std::cout << std::endl;

	std::bitset<16> bitsT;     // mybits: 0000
	int step=2;
	void *v;
	sse_trans(A, B, 300000, 8);
	//transpose32a();
	for (int x=0; x<N; x+=step)
	{
		v=B+x;
		bitsT=*(uint16_t *)(v);
		std::cout << bitsT.to_string() << std::endl;
	}
	
	return 0;
};
