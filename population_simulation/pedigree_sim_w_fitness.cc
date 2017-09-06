#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <bitset>
#include <time.h>
#include <cstring>
#include <omp.h>
#include <healpix_map.h>
#include <iomanip>
#include <map>

#include "interface.h"
#include "state.h"
#include "map_file.h"

//#define REC		//Define if you want linkage between sites, leave undefined for free recombination. 
#define MUT		//Define if you want to introduce mutation each generation.

#define LOCI	32	//Number of loci in a byte. Leave it at 32 unless you are me.
#define BINS	20

#define DEAD	-1

#define NUM_ENV		8
#define NUM_ENV_H	4

//#define TRAITS	1000

void make_map_data (arr <int> &data, const int &N) {
	if (N%12==0) {
		double d_sqrt = sqrt( N/12 );
		int i_sqrt = d_sqrt;
		if ( d_sqrt == i_sqrt ) {
			data=arr <int> (N);
			for (int x=0; x<N; x++) data[x]=x;
			return;
		} 
	}
	std::cerr << "Invalid number of Individuals for geometric structure. Try 12*i^2 for any i." << std::endl;
	exit(0);
}

class Individual
{
private:
	static uint32_t ID;
public:
	uint32_t P1, P2;
	double z;
	uint32_t e;

	Individual() {
		P1=-1;
		P2=-1;
		z=0;
		e=0;
	};

	void
	update(std::mt19937 &mt, uint32_t *Pstate[], uint32_t *Ostate[], uint32_t x) {
		uint_fast32_t r1 = mt();
		uint_fast32_t r2 = mt();
		*(Ostate[0]+x) = *(Pstate[0]+P1) & r1 | *(Pstate[1]+P1) & ~r1; 
		*(Ostate[1]+x) = *(Pstate[0]+P2) & r2 | *(Pstate[1]+P2) & ~r2;
	};

	void
	update2(const uint32_t &r1, const uint32_t &r2, uint32_t *Pstate[], uint32_t *Ostate[], const size_t &x)//, std::uniform_int_distribution<int> &r32, std::uniform_int_distribution<int>&r10, std::mt19937 &mt ){
	{
		*(Ostate[0]+x) = *(Pstate[0]+P1) & r1 | *(Pstate[1]+P1) & ~r1; 
		*(Ostate[1]+x) = *(Pstate[0]+P2) & r2 | *(Pstate[1]+P2) & ~r2;
	};

};

inline void
half_update1(const Individual &ind, const uint32_t &r1, uint32_t *Pstate[], uint32_t &Ostate)
{
	Ostate = *(Pstate[0]+ind.P1) & r1 | *(Pstate[1]+ind.P1) & ~r1; 
}

inline void
half_update2(const Individual &ind, const uint32_t &r1, uint32_t *Pstate[], uint32_t &Ostate)
{
	Ostate = *(Pstate[0]+ind.P2) & r1 | *(Pstate[1]+ind.P2) & ~r1; 
}

size_t zfunc(const double *zbreak, const size_t &disc_size, std::mt19937 &mt)
{
	std::uniform_real_distribution<double> dis(0.0, zbreak[disc_size-1]);
	double z = dis(mt);
	return size_t(std::upper_bound(zbreak, zbreak+disc_size, z)-zbreak );
}

void disc_mating(Individual *p, Individual *c, const size_t &size, const Healpix_Map <int> &map, std::mt19937 &mt)
{
	std::vector <int> disc;
	int disc_size;

	//A five hundred mile dispersal disc
	//double rad=0.050478006;

	double rad=std::max(atan(2*sqrt(8.0/size) ), 0.050478006);
	double *zbreak=new double [ int(pow(tan(rad)/2, 2)*size*2) ];
	std::cerr << "size :" << int(pow(tan(rad)/2, 2)*size*2) << std::endl;  
	for (size_t x=0; x<size; ++x)
        {
		map.query_disc( map.pix2ang(x), rad, disc);
		disc_size=disc.size();
		disc_size=disc.size();
		//TODO These calls are a little expensive, cut them?
		zbreak[0]=exp(p[map[disc[0]]].z);
		for (size_t y=1; y<size_t (disc_size); ++y) {
			//TODO These calls are a little expensive, cut them?
			zbreak[y]=zbreak[y-1]+exp(p[map[disc[y]]].z);
		}	
                c[x].P1=map[disc[zfunc(zbreak, disc_size, mt)]];
                c[x].P2=map[disc[zfunc(zbreak, disc_size, mt)]];
        }
	delete [] zbreak;
}

void
rand_mating(Individual *p, Individual *c, const size_t &size, std::mt19937 &mt)
{
	double *zbreak=new double [size];
	zbreak[0]=exp(p[0].z);
	for (size_t y=1; y<size; ++y) {
		zbreak[y]=zbreak[y-1]+exp(p[y].z);
	}	

        for (size_t x=0; x<size; ++x)
        {
                c[x].P1=zfunc(zbreak, size, mt);
                c[x].P2=zfunc(zbreak, size, mt);
        }
	delete [] zbreak;
}

void print_head(std::ostream &out, const int &length)
{
	out << "@HD	VN:0.4	SO:coordinate\n@SQ	SN:scaffold_1	LN:" << length*32 << "\n";
}

void print_names(uint32_t N, uint32_t *mask, uint32_t t, std::ostream &out)
{
	
	out << "@NAME:FILES	VERSION:TYPED	FORMAT:TEXT\n";
	out << "@FILE_NAME	SAMPLE_NAME	...\n";
	out << "../analysis_files/mpileup.txt.gz";
	for (size_t x=0; x<N; ++x){
		if ( mask[x >> 5] & (1 << (x & 0x1F) ) ) 
			out << '\t' << t*N+x;
	}
	out << "\n@END_TABLE\n";
}

void print_cord(uint32_t N, uint32_t *mask, uint32_t t, const Healpix_Map <int> &map, std::ostream &out)
{
	
	out << "@NAME:TRAITS	VERSION:TYPED	FORMAT:TEXT\n";
	out << "@SAMPLE_NAME	LONG	LAT\n";
	for (size_t x=0; x<N; ++x){
		pointing p=map.pix2ang(x);
		if ( mask[x >> 5] & (1 << (x & 0x1F) ) ) 
			out << t*N+x << '\t' <<  p.theta << '\t' << p.phi << std::endl;
	}
	out << "@END_TABLE\n";
}

void print_trait(uint32_t N, uint32_t *mask, uint32_t t, Individual *ind, std::ostream &out)
{
	out << "@NAME:TRAITS	VERSION:TYPED	FORMAT:TEXT\n";
	out << "@SAMPLE_NAME	TRAIT1\n";
	for (size_t x=0; x<N; ++x){
		if ( mask[x >> 5] & (1 << (x & 0x1F) ) ) 
			out << t*N+x << '\t' <<  ind[x].z << std::endl;
	}
	out << "@END_TABLE\n";
}



class Gcounter {
	private:
	std::geometric_distribution<int> *geom;
	std::mt19937 *mt;
	int N;
	public:
	Gcounter(std::mt19937 &_mt, double rate)
	{
		mt=&_mt;
		if (rate>1) rate=1;
		geom=new std::geometric_distribution<int> (rate);
		N=(*geom)(*mt);
	}
	~Gcounter(void)
	{
		delete geom;
	}	
	bool draw(void){
		//--N==0 ? return false : {N=(*geom)(*mt); return true;};
		if (N--!=0) return false;
		else {
			N=(*geom)(*mt);
			return true;
		}
//		return false;
	}
};

void mutate (uint32_t *ind, std::uniform_int_distribution<int> &r32, std::uniform_int_distribution<int> &randN, std::poisson_distribution<int> &poisson, std::mt19937 &mt)
{
	const uint32_t MASK[LOCI]={0x00000001,	0x00000002,	0x00000004,	0x00000008,
				0x00000010,	0x00000020,	0x00000040,	0x00000080,
				0x00000100,	0x00000200,	0x00000400, 	0x00000800,	
				0x00001000, 	0x00002000, 	0x00004000, 	0x00008000,
				0x00010000,	0x00020000,	0x00040000,	0x00080000,
				0x00100000,	0x00200000,	0x00400000,	0x00800000,
				0x01000000,	0x02000000,	0x04000000, 	0x08000000,	
				0x10000000, 	0x20000000, 	0x40000000, 	0x80000000};

	int mut_num=poisson(mt);
	for (int x=0; x<mut_num; ++x){
		ind[randN(mt)]^=(MASK[r32(mt)]);
	}
}

int N_REC;

void
recombine2 (uint32_t *recs, uint32_t *r_end, std::uniform_int_distribution<int> &r32, Gcounter &gcounter, std::mt19937 &mt)
{
	const uint32_t MASK[32]={0x00000001,	0x00000003,	0x00000007,	0x0000000f,
				0x0000001f,	0x0000003f,	0x0000007f,	0x000000ff,
				0x000001ff,	0x000003ff,	0x000007ff, 	0x00000fff,	
				0x00001fff, 	0x00003fff, 	0x00007fff, 	0x0000ffff,
				0x0001ffff,	0x0003ffff,	0x0007ffff,	0x000fffff,
				0x001fffff,	0x003fffff,	0x007fffff,	0x00ffffff,
				0x01ffffff,	0x03ffffff,	0x07ffffff, 	0x0fffffff,	
				0x1fffffff, 	0x3fffffff, 	0x7fffffff, 	0xffffffff};

	/*const uint32_t MASK_BAD[32]={0x7ffffff,	0x7fffffff,	0x7ffffffd,	0x7ffffff9,
				0x7ffffff1,	0x7fffffe1,	0x7fffffc1,	0x7fffff81,
				0x7fffff01,	0x7ffffe01,	0x7ffffc01, 	0x7ffff801,	
				0x7ffff001, 	0x7fffe001, 	0x7fffc001, 	0x7fff8001,
				0x7fff0001,	0x7ffe0001,	0x7ffc0001,	0x7ff80001,
				0x7ff00001,	0x7fe00001,	0x7fc00001,	0x7f800001,
				0x7f000001,	0x7e000001,	0x7c000001, 	0x78000001,	
				0x70000001, 	0x60000001, 	0x40000001, 	0x00000001};
	*/
	uint32_t *r=recs;
	while (r!=r_end){
		if( gcounter.draw() ) {
			if (*r & 0x80000000) {
				*r=(MASK[r32(mt)]);
			}
			else {
				*r=~MASK[r32(mt)];
			}
		} else {
			if (*r & 0x80000000) *r=0xFFFFFFFF;
			else *r=0;
		}
		++r;
	}
}

void
recombine (uint32_t &r, std::uniform_int_distribution<int> &ri, std::mt19937 &mt)
{
	//Switch based on the value of the lowest bit
	const uint32_t MASK[32]={0x00000001,	0x00000003,	0x00000007,	0x0000000f,
					0x0000001f,	0x0000003f,	0x0000007f,	0x000000ff,
					0x000001ff,	0x000003ff,	0x000007ff, 	0x00000fff,	
					0x00001fff, 	0x00003fff, 	0x00007fff, 	0x0000ffff,
					0x0001ffff,	0x0003ffff,	0x0007ffff,	0x000fffff,
					0x001fffff,	0x003fffff,	0x007fffff,	0x00ffffff,
					0x01ffffff,	0x03ffffff,	0x07ffffff, 	0x0fffffff,	
					0x1fffffff, 	0x3fffffff, 	0x7fffffff, 	0xffffffff};
	if (r & 0x80000000) r=~(MASK[ri(mt)]);
	else r=MASK[ri(mt)];
}

void
rand_rec(uint32_t &r, std::uniform_int_distribution<int> &r32, std::uniform_int_distribution<int> &r10, std::mt19937 &mt)
{
	if (!r10(mt)){
		recombine(r, r32, mt);
	} else {
		if (r & 0x80000000) r=0xFFFFFFFF;
		else r=0;	
	}
}

class Trait
{
private:
public:
	uint32_t *k, *x2k, *x2z;
	double **W1, **W2;
	uint32_t ntraits;
	std::map <uint32_t, uint32_t> k2x;
	Trait (const size_t &NLOCI, const size_t &NTRAITS, std::mt19937 &mt, const double &var_a, const int &mode) 
	{

		ntraits=NTRAITS;

		k=new uint32_t[NLOCI];
		memset(k, 0, sizeof(uint32_t)*NLOCI);

		x2k=new uint32_t[NTRAITS];
		x2z=new uint32_t[NTRAITS];

		W1=new double *[NTRAITS];
		W2=new double *[NTRAITS];

		std::vector<uint32_t> v(NLOCI*32);
		for (size_t x=0; x<(NLOCI-1)*32; x++) v[x]=x;
		std::random_shuffle(v.begin(), v.end());

		std::vector<uint32_t> ks(v.begin(), v.begin()+NTRAITS);	

		std::sort (ks.begin(), ks.end() );	

		double mean_a=-var_a, mean_d=var_a;
		
		switch (mode)
		{
			case 4:
			case 5:
				mean_a=0;
				break;
			case 2:
			case 3:
				mean_d=0;
				break;
			case 0:
			case 1:
			default:
				break;
			
		}

		std::normal_distribution<double> a_norm(mean_a, var_a);
		std::normal_distribution<double> d_norm(mean_d, 0.75*var_a);

		for (size_t x=0; x<NTRAITS; ++x)
		{
			uint32_t K=ks[x] >> 5;
			uint32_t z=ks[x] & 0x0000001f;

			//std::cerr << "Trait " << ks[x] << " byte: " << K << ", bit" << z << std::endl;

			if ( k2x.find(K) == k2x.end() ) {	
				k2x[K]=x;
			}

			k[K]+=(1 << z);

			x2k[x]=K;
			x2z[x]=z;

			W1[x]=new double[3];
			W1[x][0]=-a_norm(mt);
			switch (mode)
			{
				case 0:
				case 2:
				case 4:
					W1[x][1]=-W1[x][0]; 
					break;
				case 1:
				case 3:
				case 5:
					W1[x][1]=d_norm(mt);
					break;
				default:
					W1[x][1]=0;
					break;
			}
			W1[x][2]=-W1[x][0];
			
			W2[x]=new double[3];
			W2[x][0]=-W1[x][0];
			W2[x][1]=-W1[x][1];
			W2[x][2]=-W1[x][2];
		}	
	};
	~Trait(void)
	{
		delete [] k;
		delete [] x2k;
		delete [] x2z;
		
		for (size_t x=0; x<ntraits; x++) 
		{
			delete [] W1[x];
			delete [] W2[x];
		}

		delete [] W1;
		delete [] W2;
		k2x.clear();
	};
	
};


void inc(Individual *ind, Trait &trait, uint32_t *P[], const size_t &N, const size_t &k) 
{
	const uint32_t *P1=P[0];
	const uint32_t *P2=P[1];
	static uint32_t MASK[32]={0x00000001,	0x00000002,	0x00000004,	0x00000008,
					0x00000010,	0x00000020,	0x00000040,	0x00000080,
					0x00000100,	0x00000200,	0x00000400, 	0x00000800,	
					0x00001000, 	0x00002000, 	0x00004000, 	0x00008000,
					0x00010000,	0x00020000,	0x00040000,	0x00080000,
					0x00100000,	0x00200000,	0x00400000,	0x00800000,
					0x01000000,	0x02000000,	0x04000000, 	0x08000000,	
					0x10000000, 	0x20000000, 	0x40000000, 	0x80000000};
	if(trait.k[k]!=0) 
	{
		uint32_t Z2=trait.k[k];
		uint32_t x=trait.k2x[k];
		for (uint8_t Z=0; Z<32; Z++)
		{
			uint32_t mask=MASK[Z];
			if ( (mask & Z2) != 0)
			{
				double W1[3], W2[3];
				std::memcpy(W1, trait.W1[x], sizeof(double)*3);
				std::memcpy(W2, trait.W2[x], sizeof(double)*3);
			
				for (size_t i=0; i<N; i++)
				{
					//if (ind[i].e >= NUM_ENV_H) 
					ind[i].z+=W1[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
					//else ind[i].z+=W2[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
				}
				x++;
			}
		}	
	}
}

uint32_t Individual::ID=1;

void 
initR (uint32_t **Rand, const size_t & N, std::mt19937 &mt)
{
	uint32_t *r1_ptr, *r2_ptr, *r1_end;
	
	std::uniform_int_distribution<int> r32(0,32);
	std::uniform_int_distribution<int> r10(0,20);

	r1_ptr=Rand[0];
	r2_ptr=Rand[1];
	r1_end=Rand[0]+N;
	while (r1_ptr!=r1_end)
	{
		*r1_ptr=mt();
		*r2_ptr=mt();
		rand_rec(*r1_ptr, r32, r10, mt);
		rand_rec(*r2_ptr, r32, r10, mt);
		++r1_ptr;
		++r2_ptr;
	}	
}

double get_h(uint32_t **Results, 
	const size_t &generation_size,
	const uint32_t &Z)
{
	double H=0;
	const uint32_t MASK[LOCI]={0x00000001,	0x00000002,	0x00000004,	0x00000008,
					0x00000010,	0x00000020,	0x00000040,	0x00000080,
					0x00000100,	0x00000200,	0x00000400, 	0x00000800,	
					0x00001000, 	0x00002000, 	0x00004000, 	0x00008000,
					0x00010000,	0x00020000,	0x00040000,	0x00080000,
					0x00100000,	0x00200000,	0x00400000,	0x00800000,
					0x01000000,	0x02000000,	0x04000000, 	0x08000000,	
					0x10000000, 	0x20000000, 	0x40000000, 	0x80000000};
	{
		for (size_t y=0; y<generation_size; ++y){
			H+=( (Results[0][y] & MASK[Z])!=(Results[1][y] & MASK[Z]) );
		}
	}
	return H/double(generation_size);
}

uint32_t 
count_e(const Individual *ind, const size_t &N, const size_t &e)
{
	uint32_t sum=0;
	for (size_t y=0; y<N; ++y) sum+=ind[y].e==e;
	return sum;
}
double get_freq(uint32_t **Results, 
	const size_t &generation_size,
	const uint32_t &Z)
{
	double p=0;
	const uint32_t MASK[LOCI]={0x00000001,	0x00000002,	0x00000004,	0x00000008,
					0x00000010,	0x00000020,	0x00000040,	0x00000080,
					0x00000100,	0x00000200,	0x00000400, 	0x00000800,	
					0x00001000, 	0x00002000, 	0x00004000, 	0x00008000,
					0x00010000,	0x00020000,	0x00040000,	0x00080000,
					0x00100000,	0x00200000,	0x00400000,	0x00800000,
					0x01000000,	0x02000000,	0x04000000, 	0x08000000,	
					0x10000000, 	0x20000000, 	0x40000000, 	0x80000000};
	//std::cerr << " bit " << Z << "	";
	for (size_t y=0; y<generation_size; ++y){
	//	std::cerr << int((Results[0][y] & MASK[Z])!=0) << "	"<< int( (Results[1][y] & MASK[Z])!=0) << "	";
		p+=( ( (Results[0][y] & MASK[Z])!=0)+( (Results[1][y] & MASK[Z])!=0) );
	}
	//std::cerr << std::endl;
	return p/double(generation_size*2.);
}

double get_freq_if(uint32_t **Results, 
	const size_t &generation_size,
	const Individual *ind,
	const uint32_t e,
	const uint32_t &Z)
{
	double p=0;
	double N=0;
	const uint32_t MASK[LOCI]={0x00000001,	0x00000002,	0x00000004,	0x00000008,
					0x00000010,	0x00000020,	0x00000040,	0x00000080,
					0x00000100,	0x00000200,	0x00000400, 	0x00000800,	
					0x00001000, 	0x00002000, 	0x00004000, 	0x00008000,
					0x00010000,	0x00020000,	0x00040000,	0x00080000,
					0x00100000,	0x00200000,	0x00400000,	0x00800000,
					0x01000000,	0x02000000,	0x04000000, 	0x08000000,	
					0x10000000, 	0x20000000, 	0x40000000, 	0x80000000};
	for (size_t y=0; y<generation_size; ++y){
		if (ind[y].e==e)
		{
			p+=( (Results[0][y] & MASK[Z])!=0)+( (Results[1][y] & MASK[Z])!=0);
			N++;
		}
	}
	return p/double(N*2.);
}

double get_freq_if_2(uint32_t **Results, 
	const size_t &generation_size,
	const Individual *ind,
	const uint32_t e1,
	const uint32_t e2,
	const uint32_t &Z)
{
	double p=0;
	double N=0;
	const uint32_t MASK[LOCI]={0x00000001,	0x00000002,	0x00000004,	0x00000008,
					0x00000010,	0x00000020,	0x00000040,	0x00000080,
					0x00000100,	0x00000200,	0x00000400, 	0x00000800,	
					0x00001000, 	0x00002000, 	0x00004000, 	0x00008000,
					0x00010000,	0x00020000,	0x00040000,	0x00080000,
					0x00100000,	0x00200000,	0x00400000,	0x00800000,
					0x01000000,	0x02000000,	0x04000000, 	0x08000000,	
					0x10000000, 	0x20000000, 	0x40000000, 	0x80000000};
	for (size_t y=0; y<generation_size; ++y){
		if (ind[y].e==e1 || ind[y].e==e2)
		{
			p+=( (Results[0][y] & MASK[Z])!=0)+( (Results[1][y] & MASK[Z])!=0);
			N++;
		}
	}
	return p/double(N*2.);
}

void 
iter (Individual *parents, Individual *children,
	const size_t &generation_size, const char &type,
	const Healpix_Map <int> & map, std::mt19937 &mt)
{
	if (type=='g')
	disc_mating(parents, children, generation_size, map, mt);
	else 
	rand_mating(parents, children, generation_size, mt);
}

void 
open_ped(std::ostream &out)
{
	out << "Time" << '\t' << "Individual ID" << '\t' << "Paternal ID" << '\t' << "Maternal ID" << '\t' << "Sex" << '\t' <<  "z" << std::endl;
}

void 
push_ped(Individual *ind,
	const size_t &generation_size, 
	const size_t &time,
	std::ostream &out)
{
	for (size_t x=0; x<generation_size; ++x)
	{
		out << time << '\t' << x+(time+1)*generation_size << '\t' << (ind->P1)+time*generation_size << '\t' << (ind->P2)+time*generation_size << '\t' << 3 << '\t' << ind->z;
		out << std::endl;
		++ind;
	}
}

void
print_results4(std::ostream &out, const uint32_t &N, uint32_t **Results)
{ 
	for (size_t b=0; b<32; ++b) {
		for (size_t y=0; y<N; ++y){
			std::bitset <32> bits0(Results[0][y]);
			std::bitset <32> bits1(Results[1][y]);
			if (!bits0[b]) {
				if (!bits1[b]) {
					out << "	0	0";
				} else {
					out << "	0	1";
				}
			} else {
				if (!bits1[b]) {
					out << "	1	0";
				} else {
					out << "	1	1";
				}
			}
		}
		out << std::endl;
	}
}

void
print_results_bin(std::ostream &out, const uint32_t &N, uint32_t **Results)
{ 
	for (size_t y=0; y<N; ++y){
		out.write((char *)(&(Results[0][y])), sizeof(uint32_t) );
		out.write((char *)(&(Results[1][y])), sizeof(uint32_t) );
	}
}

void
print_results3(std::ostream &out, const uint32_t &k, const uint32_t &N, uint32_t **Results)
{ 
for (size_t K=0; K<k; ++K){
	for (size_t b=0; b<32; ++b) {
		for (size_t y=0; y<N; ++y){
			std::bitset <32> bits0(Results[0][N*K+y]);
			std::bitset <32> bits1(Results[1][N*K+y]);
			if (!bits0[b]) {
				if (!bits1[b]) {
					out << "	0	0";
				} else {
					out << "	0	1";
				}
			} else {
				if (!bits1[b]) {
					out << "	1	0";
				} else {
					out << "	1	1";
				}
			}
		}
		out << std::endl;
	}
}
}

void
print_results2(std::ostream &out, const uint32_t &k, const uint32_t &N, uint32_t **Results)
{ 
for (size_t K=0; K<k; ++K){
	for (size_t b=0; b<32; ++b) {
		out << "scaffold_1	" << b+K*32+1;
		for (size_t y=0; y<N; ++y){
			std::bitset <32> bits0(Results[0][N*K+y]);
			std::bitset <32> bits1(Results[1][N*K+y]);
			if (!bits0[b]) {
				if (!bits1[b]) {
					out << "	0";
				} else {
					out << "	1";
				}
			} else {
				if (!bits1[b]) {
					out << "	1";
				} else {
					out << "	2";
				}
			}
		}
		out << std::endl;
	}
}
}

void 
print_results(std::ostream &out, const uint32_t &k, const uint32_t &N, uint32_t **Results)
{ 
	for (size_t K=0; K<k; ++K){
		for (size_t b=0; b<32; ++b) {
			out << "scaffold_1	" << b+K*32+1 << "	G";
			for (size_t y=0; y<N; ++y){
				std::bitset <32> bits0(Results[0][N*K+y]);
				std::bitset <32> bits1(Results[1][N*K+y]);
				if (!bits0[b]) {
					if (!bits1[b]) {
						out << "	12	TtTtTtTtTtTt	~~~~~~~~~~~";
					} else {
						out << "	12	TtTtTtAaAaAa	~~~~~~~~~~~";
					}
				} else {
					if (!bits1[b]) {
						out << "	12	TtTtTtAaAaAa	~~~~~~~~~~~";
					} else {
						out << "	12	AaAaAaAaAaAa	~~~~~~~~~~~";
					}
				}
			}
			out << std::endl;
		}
	}
}

/*print_cov is bad?*/
void 
print_cov(const Trait &traits, State &sites, uint32_t *mask2, std::ostream &out) 
{
	const uint32_t MASK[LOCI]={0x00000001,	0x00000002,	0x00000004,	0x00000008,
					0x00000010,	0x00000020,	0x00000040,	0x00000080,
					0x00000100,	0x00000200,	0x00000400, 	0x00000800,	
					0x00001000, 	0x00002000, 	0x00004000, 	0x00008000,
					0x00010000,	0x00020000,	0x00040000,	0x00080000,
					0x00100000,	0x00200000,	0x00400000,	0x00800000,
					0x01000000,	0x02000000,	0x04000000, 	0x08000000,	
					0x10000000, 	0x20000000, 	0x40000000, 	0x80000000};

	out << "@NAME:TRAITS	VERSION:TYPED	FORMAT:TEXT\n";
	out << "@POS";
	for (size_t x=0; x<sites.sample_size() ; ++x){
		out << '\t' << "beta_" << x;
	}
	for (size_t x=0; x<sites.sample_size() ; ++x){
		out << '\t' << "delta_" << x;
	}
	out << "\n";

	uint32_t *P[2];
	uint32_t *P1=new uint32_t [sites.sample_size()];
	uint32_t *P2=new uint32_t [sites.sample_size()];

	P[0]=P1;
	P[1]=P2;

	sites.cache();
	sites.rewind();

	for (size_t y=0; y<traits.ntraits; y++){
		uint32_t z=traits.x2z[y];
		uint32_t k=traits.x2k[y];
		out << k*32+z;
		/* TODO Fix.*/
		//std::cerr << k << ", " << z << std::endl;
		sites.uncompress(P1, P2, k);
		const uint32_t mask=MASK[z];

		const double a=traits.W1[y][2];
		const double d=traits.W1[y][1];
		const double p=get_freq(P, sites.sample_size(), z);  
		const double q=1.-p;
		const double D=2*p*q-get_h(P, sites.sample_size(), z);  
		const double alpha=a+d*(q-p);
		double beta[3]={-p*2*alpha,(q-p)*alpha,2*q*alpha};
		double delta[3]={d*(-2*p*p+D),d*(2*p*q+D),d*(-2*q*q+D)};

		for (size_t i=0; i<sites.sample_size(); ++i){
			if ( mask2[i >> 5] & (1 << (i & 0x1F) ) ) 
				out << '\t' << beta[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
		}
		for (size_t i=0; i<sites.sample_size(); ++i){
			if ( mask2[i >> 5] & (1 << (i & 0x1F) ) ) 
				out << '\t' << delta[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
		}
		out << "\n";
	}
	sites.rewind();
	delete P1;
	delete P2;
}

/*print_bins is bad?*/
void 
print_bins(const Trait &trait, State &sites, std::ostream &out)
{
	out << "@NAME:TRAITS	VERSION:TYPED	FORMAT:TEXT\n";
	out << "@FREQ	A	D	N\n";
	double a[BINS];
	double d[BINS];
	double N[BINS];

	memset(a, 0, sizeof(double)*BINS);
	memset(d, 0, sizeof(double)*BINS);
	memset(N, 0, sizeof(double)*BINS);

	uint32_t *P[2];
	uint32_t *P1=new uint32_t [sites.sample_size()];
	uint32_t *P2=new uint32_t [sites.sample_size()];

	P[0]=P1;
	P[1]=P2;

	//sites.rewind();

	/*for (size_t K=0; K<2; ++K)
	{
	sites.uncompress(P[0], P[1]);
	for (size_t z=0; z<32; ++z)
		std::cerr << get_freq(P, 12, z) << std::endl;
	}*/

	sites.cache();
	sites.rewind();

	for (size_t x=0; x<trait.ntraits; ++x)
	{
		uint32_t z=trait.x2z[x];
		uint32_t k=trait.x2k[x];

		/*TODO fix*/
		//std::cerr << "uncompressing " << k << std::endl;
		sites.uncompress(P1, P2, k);
		
		double p=get_freq(P, sites.sample_size(), z);  

		//std::cerr << " byte " << k << "bit " << z << "	"<< p << "	" << trait.W1[x][0] << std::endl;

		uint32_t b=round(p*double(BINS-1) );

		a[b]+=trait.W1[x][0];
		d[b]+=trait.W1[x][1];
		N[b]+=1;
		
	}
	for (size_t x=0; x<BINS; ++x){
		if (N[x]>0.5)
			out << double(x)/double(BINS) << '\t' <<  a[x]/N[x] << '\t' << d[x]/N[x] << '\t' << N[x] << std::endl;
		else 
			out << double(x)/double(BINS) << '\t' <<  "NaN" << '\t' << "NaN" << '\t' << 0 << std::endl;
	}
	out << "@END_TABLE\n";
	sites.rewind();
	delete P1;
	delete P2;

}
///TODO This is becoming something serious, clean up the code!!


//TODO: Warning, not thread safe!!
std::map <uint32_t, double> harmonic_array;

double harmonic(const uint32_t &size)
{
	std::map <uint32_t, double>::iterator ind=harmonic_array.find(size);
	if ( ind == harmonic_array.end() ) {
		double sum=0;
		for (uint32_t x=1; x<(size+1); x++) sum+=1./double(x);
		harmonic_array[size]=sum;
		return sum;
	} else {
		return ind->second;
	}
	return -1;
}

/*print_stats is bad*/
void
print_stats(const Trait &trait, State &sites, const Individual *ind, std::ostream &out)
{
	out << "@NAME:TRAITS	VERSION:TYPED	FORMAT:TEXT\n";
	out << "@POS";
	for (size_t z=0; z<NUM_ENV; ++z)
	{
		out << "	PI_" << z;
	}
	for (size_t z=0; z<NUM_ENV; ++z)
	{
		out << "	THETA_" << z;
	}
	for (size_t w=0; w<NUM_ENV; ++w)
	{
		for (size_t z=w+1; z<NUM_ENV; ++z)
		{
			out << "	PI_" << w << "_" << z;
		}
	}
	for (size_t w=0; w<NUM_ENV; ++w)
	{
		for (size_t z=w+1; z<NUM_ENV; ++z)
		{
			out << "	FST_" << w << "_" << z;
		}
	}
	out << std::endl;
	uint32_t *P[2];
	uint32_t *P1=new uint32_t [sites.sample_size()];
	uint32_t *P2=new uint32_t [sites.sample_size()];

	P[0]=P1;
	P[1]=P2;

	size_t window=1000;

	double pi[NUM_ENV]={0}, t[NUM_ENV];
	double PI[NUM_ENV][NUM_ENV]={0}, T[NUM_ENV][NUM_ENV];
	double f_st[NUM_ENV][NUM_ENV]={0};
	double f_st_n[NUM_ENV][NUM_ENV]={0};
	uint32_t L[NUM_ENV]={0};
	size_t w=0;

	std::cerr << "print stats " << sites.sample_size() << std::endl;

	sites.cache();
	sites.rewind();

	for (size_t x=0; x<sites.genome_size(); ++x)
	{
		/*TODO Fix.*/
		sites.uncompress(P1, P2, x);
		for (size_t y=0; y<32; ++y)
		{
			for (size_t z=0; z<NUM_ENV; ++z)
			{
				t[z]=get_freq_if(P, sites.sample_size(), ind, z, y);
				t[z]=t[z]*(1.-t[z]);
				pi[z] +=t[z];
				L[z]+=(t[z]>0 && t[z]<1);
			}
			for (size_t w=0; w<NUM_ENV; ++w)
			{
				for (size_t z=w+1; z<NUM_ENV; ++z)
				{
					T[w][z]=get_freq_if_2(P, sites.sample_size(), ind, w, z, y);
					T[w][z]=T[w][z]*(1.-T[w][z]);
					PI[w][z]+=T[w][z];
					if (T[w][z]>0 )
					{
						f_st[w][z]+=(2.*T[w][z]-t[w]-t[z] )/(2.*T[w][z]);
						f_st_n[w][z]+=1.;
					}
				}
			}
			if (w==window)
			{
				out  <<  x*32+y; 
				for (size_t z=0; z<NUM_ENV; ++z)
				{
					out << '\t' << pi[z]/window;
					pi[z]=0;
				}
				for (size_t z=0; z<NUM_ENV; ++z)
				{
					double a=harmonic(count_e(ind, sites.sample_size(), z) );
					out << '\t' << L[z]/(2.*a)/window;
					L[z]=0;
				}
				for (size_t w=0; w<NUM_ENV; ++w)
				{
					for (size_t z=w+1; z<NUM_ENV; ++z)
					{
						out << '\t' << PI[w][z]/window;
						PI[w][z]=0;
					}
				}
				for (size_t w=0; w<NUM_ENV; ++w)
				{
					for (size_t z=w+1; z<NUM_ENV; ++z)
					{
						out << '\t' << f_st[w][z]/f_st_n[w][z];
						f_st[w][z]=0;
						f_st_n[w][z]=0;
					}
				}
				out << std::endl;
				w=0;
			}
			w++;
		}
	}
	out << "@END_TABLE\n";
	sites.rewind();
	delete P1;
	delete P2;
}

int 
main(int argc, char *argv[] )
{

	Environment env;
	
	std::string i_file;
	int mode=0;

	bool binary=false, p_pedigree=false, p_names=false, p_traits=false, p_geography=false, p_dist=false, p_header=false, p_matrix=false, noselect=false, noprint=false;
	bool p_stats=false;

	int N=1200, t=10, k=500, n_trait=50, skip=100, sub_sample_size=-1;

	uint32_t *mask;

	double r=5;			//recombination events per Individual per generation
	double pi=0.001;		//equilibrium diversity
	double var=0.01;		//trait var
	double theta=1;		//equilibrium diversity

	char type='g';			//geographic or random?

	env.set_name("pedigree_sim");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("A MPI capable forward time evolution simulation. Please direct questions to matthew.s.ackerman@gmail.com");

	env.optional_arg('N',"number", 	N,	"please provide a number.", "number of individuals in the populations.");
	env.optional_arg('g',"gen", 	t,	"please provide a number.", "number of generations the simulation runs.");
	env.optional_arg('s',"sites", 	k,	"please provide a number.", "number of sites/32.");
	env.optional_arg('v',"var", 	n_trait,"please provide a number.", "number of loci effecting the trait.");
	env.optional_arg('c',"chiasma",	r,	"please provide a number.", "number of recombinations events per individual per generation.");
	env.optional_arg('e',"equil", 	pi,	"please provide a number.", "equilibrium nucleotide diversity.");
	env.optional_arg('y',"type", 	type,	"please provide a number.", "type of mating < r | g > for random or geographic.");
	env.optional_arg('a',"theta", 	theta,	"please provide a number.", "angle separating environments");
	env.optional_arg('k',"skip", 	skip,	"please provide a number.", "number of generations to skip between printing files");
	env.optional_arg('i',"init", 	i_file,	"please provide a number.", "state file to load");
	env.optional_arg('S',"sigma", 	var,	"please provide a number.", "variance of trait per snp");
	env.optional_arg('Q',"mode", 	mode,	"please provide a number.", "undocumented mode.");
	env.optional_arg('1',"sample", 	sub_sample_size,	"please provide a number.", "sub sample.");

	env.flag(	'O',"opts", 	&env, 		&flag_options, "an error occurred while displaying the help message", "Prints a list of available options");		//DONE
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message", "Prints this message");				//DONE
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occurred while displaying the version message", "Prints the program version");		//DONE
	env.flag(	'b',"binary", 	&binary, 	&flag_set, 	"an error occurred while displaying the version message", "binary output");		//DONE
	env.flag(	'p',"pedigree", &p_pedigree,&flag_set, 	"an error occurred while displaying the version message", "prints the pedigree in the file 'pedigree.txt' ");		//DONE
	env.flag(	'n',"name", &p_names,	&flag_set, 	"an error occurred while displaying the version message", "prints the names of the samples.");		//DONE
	env.flag(	'G',"geography",&p_geography,&flag_set, 	"an error occurred while displaying the version message", "prints the geographic coordinates of each individual.");		//DONE
	env.flag(	't',"traits", &p_traits,	&flag_set, 	"an error occurred while displaying the version message", "prints the trait values every generation in files called t[GENERATION].txt");		//DONE
	env.flag(	'd',"dist", &p_dist,&flag_set, 	"an error occurred while displaying the version message", "prints the binned frequencies of a and d [Sorry guys, I'll figure out a way to explain this later.]");		//DONE
	env.flag(	'H',"header", &p_header,	&flag_set, 	"an error occurred while displaying the version message", "prints a header file.");		//DONE
	env.flag(	'm',"matrix", &p_matrix,	&flag_set, 	"an error occurred while displaying the version message", "prints a covariance matrix of additive and dominance deviations.");		//DONE
	env.flag(	'l',"lynch", &noselect,	&flag_set, 	"an error occurred while displaying the version message", "trait has no fitness effects.");		//DONE
	env.flag(	'o',"off", &noprint,	&flag_set, 	"an error occurred while displaying the version message", "don't print state file.");		//DONE
	env.flag(	'w',"window", &p_stats,	&flag_set, 	"an error occurred while displaying the version message", "print stats file.");		//DONE

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);


	mask=new uint32_t [(N >> 5)+1];
	memset(mask, 0xFFFFFFFF, sizeof(uint32_t)*( (N >> 5)+1) );
	if (sub_sample_size!=-1)
	{
		memset(mask, 0, sizeof(uint32_t)*( (N >> 5)+1) );
		std::vector <bool> v(N, 0);
		for (size_t x=0; x<sub_sample_size; ++x)
		{
			v[x]=1;
		}
		
		std::random_shuffle(v.begin(), v.end());
		for (size_t x=0; x<N; ++x)
		{
			if (v[x]) mask[x >> 5] |= (1 << (x & 0x1F) );
		}
	} else 
	{
		sub_sample_size=N;
	}
		
	std::mt19937 mt;

	mt.seed(time(NULL));

	if (p_header)
	{
		std::fstream header;
		header.open ("header.txt", std::fstream::out);
		print_head(header, k);
		header.close();
	}

	State Pstates(N);

	if (i_file!="")
	{
		Flat_file <State> state_file;
		state_file.open(i_file.c_str(), READ);
		Pstates=state_file.read_header();
		state_file.read(Pstates);
		state_file.close();
		N=Pstates.sample_size();
		k=Pstates.genome_size();
		std::cerr << N << ", " << k << std::endl;
		//states.clear();
	}

	State Ostates(N);

	r=r/k;

	Individual *parents=new Individual[N];
	Individual *children=new Individual[N];
	
	arr <int> data(12);
	make_map_data(data, N);
	Healpix_Map <int> map(data, RING);
	
	uint32_t **Pstate, **Ostate, **Results;
	uint32_t **Rstate;
	uint32_t **Recombination;

	Pstate=new uint32_t *[2];
	Pstate[0]=new uint32_t [N];
	Pstate[1]=new uint32_t [N];

	Ostate=new uint32_t *[2];
	Ostate[0]=new uint32_t[N];
	Ostate[1]=new uint32_t[N];

	Results=new uint32_t *[2];
	Results[0]=new uint32_t[N];
	Results[1]=new uint32_t[N];

	N_REC=0;	

#ifdef REC
	std::uniform_int_distribution<int> r32(0,31);
	std::uniform_int_distribution<int> r10(0,20);
	Gcounter gcounter(mt, r);

	Rstate=new uint32_t *[2];
	Rstate[0]=new uint32_t[N];
	Rstate[1]=new uint32_t[N];
#endif
	

#ifdef MUT
#ifndef REC
	std::uniform_int_distribution<int> r32(0,31);
#endif
	std::uniform_int_distribution<int> rN(0, N-1);
	std::poisson_distribution<int> poisson(N*pi/2./32.);
#endif
	//std::cerr << argc << std::endl;

//	check(geneology);
	Trait trait(k, n_trait, mt, var, mode);


//	std::exponential_distribution<double> distribution(3.5);


	if (i_file=="")
	{
		for (size_t K=0; K<k; ++K){
			for (size_t x=0; x<N; ++x)
			{
#ifdef MUT
				Pstate[0][x]=0x00000000;
				Pstate[1][x]=0x00000000;
#else
				if (x<N/2) 
				{
					Pstate[0][x]=0x00000000;
					Pstate[1][x]=0xFFFFFFFF;
//					Pstate[1][x]=0x00000000;
				} 
				else 
				{
//					Pstate[0][x]=0xFFFFFFFF;
					Pstate[0][x]=0x00000000;
					Pstate[1][x]=0xFFFFFFFF;
				}
#endif
			}
			Pstates.compress(Pstate[0], Pstate[1]);
		}
		Pstates.cache();
	} 
	{
	/*
	arr <int> data(12);
	make_map_data(data, 48);
	Healpix_Map <int> map2(data, RING);
	*/
	for (size_t x=0; x<N; ++x)
	{
		pointing p=map.pix2ang(x);
		//if (p.theta!=last_pheta)
		//	start_1=
		//	start=
		//int j=(map2.ang2pix(p)+start)%2;
		int j=int( ( (cos(p.theta)+1.)/2.)*(NUM_ENV) ); //^ (sin(p.phi)<0);
		parents[x].e=j;
		children[x].e=j;
	}
	}

//	std::cerr << __LINE__ << std::endl; 

	std::fstream pedigree;
	if (p_pedigree)
	{
		pedigree.open ("pedigree.txt", std::fstream::out);
		open_ped(pedigree);
	}


	for (size_t y=0; y<t; ++y){
//		#pragma omp parallel for
		std::cerr << __LINE__ << " y:" << y << std::endl; 
		iter(parents, children, N, type, map, mt);
#ifdef REC
		initR(Rstate, N, mt);
#endif
		for (size_t K=0; K<k; ++K)
		{
			Pstates.uncompress(Pstate[0], Pstate[1]);
#ifdef REC
/*
			Individual *c_it=children;
			uint32_t *O0_it=Ostate[0];
			uint32_t *O1_it=Ostate[1];
			uint32_t *R0_it=Rstate[0];
			uint32_t *R1_it=Rstate[1];
			Individual *c_end=children+N;
			while (c_it!=c_end)
			{
				half_update1(*c_it, *(++R0_it), Pstate, *(++O0_it) );
				half_update2(*(++c_it), *(++R1_it), Pstate, *(++O1_it) );
				
			}
*/

			for (size_t x=0; x<N; ++x)
			{
				children[x].update2( Rstate[0][x], Rstate[1][x], Pstate, Ostate, x);
			}
#else
			for (size_t x=0; x<N; ++x)
			{
				children[x].update(mt, Pstate, Ostate, x);
			}
#endif

#ifdef MUT
			mutate(Ostate[0], r32, rN, poisson, mt);
			mutate(Ostate[1], r32, rN, poisson, mt);
#endif
			if (!noselect) inc(children, trait, Ostate, N, K);

#ifdef REC
			recombine2((Rstate[0]), (Rstate[0]+N), r32, gcounter, mt);
			recombine2((Rstate[1]), (Rstate[1]+N), r32, gcounter, mt);
#endif
			Ostates.compress(Ostate[0], Ostate[1]);
		}
			
	
		Ostates.cache();
		std::cerr << "free :" << Ostates.get_free() << " ratio " << Ostates.compression_ratio() << std::endl;
		Pstates.clear();

		Pstates.swap(Ostates);
		std::swap(children, parents);
	
		for (size_t x=0; x<N; x++)
		{
			children[x].z=0;
		}

		if(p_pedigree){
			push_ped(parents, N, y, pedigree);
		}
		
		if (y%skip==0)
		{
		if(p_traits){
			std::fstream trait_file;
			trait_file.open("t"+std::to_string(y)+".txt", std::fstream::out);
			print_trait(N, mask, t, parents, trait_file);
			trait_file.close();
		}

		if(p_dist){
			std::fstream bin_file;
			bin_file.open("b"+std::to_string(y)+".txt", std::fstream::out);
			print_bins(trait, Pstates, bin_file);
			bin_file.close();
		}

		if(p_stats)
		{
			std::fstream stats_file;
			stats_file.open("s"+std::to_string(y)+".txt", std::fstream::out);
			print_stats(trait, Pstates, parents, stats_file);
			stats_file.close();
		}
		}
	}

	if(p_pedigree)
	{
		pedigree.close();
	}
	
	if (p_names)
	{
		std::fstream names;
		names.open ("name-file.txt", std::fstream::out);
		print_names(N, mask, t, names);
		names.close();
	}


	if(p_geography)
	{
		std::fstream map_file;
		map_file.open ("cord-file.txt", std::fstream::out);
		print_cord(N, mask, t, map, map_file);
		map_file.close();
	}

	if(p_matrix)
	{
		std::fstream cov_file;
		cov_file.open("cov.txt", std::fstream::out);
		print_cov(trait, Pstates, mask, cov_file);
		cov_file.close();
	}

	if (!noprint) 
	{
		Flat_file <State> state_file;
		Pstates.cache();
		if (binary)
			state_file.open( WRITE | BINARY );
		else 
			state_file.open( WRITE );
		State sub_states=sub_sample(Pstates, sub_sample_size, mask);
		state_file.write_header(sub_states);
		state_file.write(sub_states);
		state_file.close();
	}

	/*
	Pstates.rewind();
	for (size_t K=0; K<k; ++K)
	{
		Pstates.uncompress(Pstate[0], Pstate[1]);
		for (size_t z=0; z<32; ++z)
			std::cerr << get_freq(Pstate, N, z) << std::endl;  
	}

        std::fstream bin_file;
        bin_file.open("b"+std::to_string(-1)+".txt", std::fstream::out);
        print_bins(trait, Pstates, bin_file);
        bin_file.close();
	*/

	delete [] parents;
	delete [] children;

	delete [] Pstate[0];
	delete [] Pstate[1];
	delete [] Pstate;
	
	delete [] Ostate[0];
	delete [] Ostate[1];
	delete [] Ostate;

	delete [] Results[0];
	delete [] Results[1];
	delete [] Results;

#ifdef REC
	delete [] Rstate[0];
	delete [] Rstate[1];
	delete [] Rstate;
#endif 
	return 0;
}
