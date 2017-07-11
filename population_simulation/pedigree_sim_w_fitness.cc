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

//#include "interface.h"

#define REC
#define MUT

#define LOCI	32
//#define TRAITS	1000
class BGT 
{
private:
	size_t _size, _N;
public:
	uint32_t *chrm1, *chrm2;
	uint32_t blocks;
	BGT(uint32_t N, uint32_t k )
	{
		_size=N*k;
		_N=N;
		chrm1=new uint32_t [_size];
		chrm2=new uint32_t [_size];
	};
	void uncompress (uint32_t *a, uint32_t *b, const uint32_t &c)
	{
		memcpy(a,chrm1+c, sizeof(uint32_t)*_N );  
		memcpy(b,chrm2+c, sizeof(uint32_t)*_N );
	};
	void compress (const uint32_t *a, const uint32_t *b, const uint32_t &c)
	{
		memcpy(chrm1+c, a, sizeof(uint32_t)*_N );  
		memcpy(chrm2+c, b, sizeof(uint32_t)*_N );
	};
};

void random_mating(uint32_t *P1,  uint32_t *P2, const size_t &size)
{
	std::vector<uint32_t> m1(size), f1(size);

	std::iota(m1.begin(), m1.end(), 0);
	std::random_shuffle(m1.begin(), m1.end() );

	std::iota(f1.begin(), f1.end(), 0);
	std::random_shuffle(f1.begin(), f1.end() );

	for (size_t x=0; x<size; ++x)
	{
		P1[x]=m1[x];
		P2[x]=f1[x];	
	}
}

void random_sibs(uint32_t *P1,  uint32_t *P2, const size_t &size)
{
	//Place each number into a vector once. Randomize the vector, split into two halves, mate each pair twice.
	std::vector<uint32_t> all(size);
	std::iota(all.begin(), all.end(), 0);
	std::random_shuffle(all.begin(), all.end() );
	
	size_t oh=size/2;

	for (size_t x=0; x<oh; ++x)
	{
		P1[x]=all[x];
		P2[x]=all[x+oh];	
	}
	for (size_t x=oh; x<size; ++x)
	{
		P1[x]=all[x-oh];
		P2[x]=all[x];
	}
}

void random_consang(uint32_t *P1,  uint32_t *P2, const size_t &size)
{
	//Place each number into a vector once. Randomize the vector, split into two halves, mate each pair twice.
	std::vector<uint32_t> all(size);
	std::iota(all.begin(), all.end(), 0);
//	std::random_shuffle(all.begin(), all.end() );
	
	P1[0]=size-1;
	P2[0]=1;
	for (size_t x=1; x<size-1; ++x)
	{
		P1[x]=all[x-1];
		P2[x]=all[x+1];	
	}
	P1[size-1]=size-2;
	P2[size-1]=0;
}


/*
class coancestry {
{
	private:
	size_t _size, double *co;
	public:
	Slow_index(const size_t &size)
	{	
		_size=size;
		co=new double [size];
	}

	double
	operator () (const uint32_t &d1, const uint32_t &d2)	{return *(co+size*d1+d2);};

	double *
	operator () (const uint32_t &d1, const uint32_t &d2)	{return co+size*d1+d2;};
}
private:
public:
		
	return *(coancestry+x*size+y);
}

update_coancestry(co, x, y, P1x, P2x, P1y, P2y)
{
	if (x!=y)
		coancestry[x][y]=( co(P1x, P1y), co(P1x, P2y), co(P2x, P1y), co(P2x, P2y) )/4.;
	if (x==y)
		coancestry[x][y]=(1.+co(P1x, P2x) )/2.;
}

std::vector <std::vector <uint32_t> > make_mate_vector (double *coancestry, const size_t &size)
{
	std::vector <std::vector <uint32_t> > mate(size, vector <uint32_t > (5) );
	for (size_t x=0; x<size; ++x)
	{
		for (size_t y=x+1; y<size; ++y)
		{
			c=*(coancestry+x*size+y);
			if ( c >= 0.176 )
				mate[x][0].append(y);
			else if ( 0.088 <= c  && c < 0.176 )
				mate[x][1].append(y);
			else if ( 0.031 <= c  && c < 0.088 )
				mate[x][2].append(y);
			else if ( 0.008 <= c  && c < 0.031 )
				mate[x][3].append(y);
			else 
				mate[x][4].append(y);
		}
	}
	return mate;
}

void mating_vector(uint32_t *P1,  uint32_t *P2, double probs[5], double *coancestory, const size_t &size)
{
	std::discrete_distribution<int> distribution {self,sibs,halfsibs,1cousin,2cousin,3cousin,unrelated};
	std::vector <std::vector <uint32_t> >mate_vector=make_mate_vector(coancestory);

	for (size_t x=0; x<size; ++x)
	{
		P1=uniform(mt);
		N1=distribution(mt);
		SIZE=mate_vector[P1][N1].size();
		while (SIZE==0) SIZE=mate_vector[P1][N1++].size();
		R=uniform(mt)%SIZE;
		P2=mate_vector[P1][N1][R];
		P1[x]=P1;
		P2[x]=P2;
	}

	for (size_t x=0; x<size; ++x){
	for (size_t y=x+1; y<size; ++y){
		sum+=update_coancestry(x, y, P1[x], P2[x], P1[y], P2[y]);
	}
	}
	for (size_t x=0; x<size; ++x){
		co-=mean;
	}
}*/


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

	Individual() {
		P1=-1;
		P2=-1;
		z=1;
	};

	void
	update(std::mt19937 &mt, uint32_t *Pstate[], uint32_t *Ostate[], uint32_t x) {
		uint_fast32_t r1 = mt();
		uint_fast32_t r2 = mt();
		*(Ostate[0]+x) = *(Pstate[0]+P1) & r1 | *(Pstate[1]+P1) & ~r1; 
		*(Ostate[1]+x) = *(Pstate[0]+P2) & r2 | *(Pstate[1]+P2) & ~r2;
	};

	void
	update2(uint32_t &r1, uint32_t &r2, uint32_t *Pstate[], uint32_t *Ostate[], uint32_t x, std::uniform_int_distribution<int> &r32, std::uniform_int_distribution<int>&r10, std::mt19937 &mt ){
		*(Ostate[0]+x) = *(Pstate[0]+P1) & r1 | *(Pstate[1]+P1) & ~r1; 
		*(Ostate[1]+x) = *(Pstate[0]+P2) & r2 | *(Pstate[1]+P2) & ~r2;
	};
};

void disc_mating(uint32_t *P1,  uint32_t *P2, const size_t &size, const Healpix_Map <int> &map)
{
	std::vector <int> disc;
	int disc_size;
	double rad=0.050478006;
        for (size_t x=0; x<size; ++x)
        {
		map.query_disc( map.pix2ang(x), rad, disc);
		disc_size=disc.size();
                P1[x]=map[disc[rand() % disc_size]];
                P2[x]=map[disc[rand() % disc_size]];
        }
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
	double rad=0.050478006;

	//double rad=atan(2*sqrt(8.0/size) );
	double *zbreak=new double [ int(pow(tan(rad)/2, 2)*size*2) ];
        for (size_t x=0; x<size; ++x)
        {
		map.query_disc( map.pix2ang(x), rad, disc);
		disc_size=disc.size();
		disc_size=disc.size();
		zbreak[0]=exp(p[map[disc[0]]].z);
		for (size_t y=1; y<disc_size; ++y) {
			zbreak[y]=zbreak[y-1]+exp(p[map[disc[y]]].z);
		}	
                c[x].P1=map[disc[zfunc(zbreak, disc_size, mt)]];
                c[x].P2=map[disc[zfunc(zbreak, disc_size, mt)]];
        }
	delete zbreak;
}

void print_head(std::ostream &out, const int &length)
{
	out << "@HD	VN:0.4	SO:coordinate\n@SQ	SN:scaffold_1	LN:" << length*32 << "\n";
}

void print_names(int N, int t, std::ostream &out)
{
	
	out << "@NAME:FILES	VERSION:TYPED	FORMAT:TEXT\n";
	out << "@FILE_NAME	SAMPLE_NAME	...\n";
	out << "../analysis_files/mpileup.txt.gz";
	for (size_t x=0; x<N; ++x){
		out << '\t' << t*N+x;
	}
	out << "\n@END_TABLE\n";
}

void print_cord(int N, int t, const Healpix_Map <int> &map, std::ostream &out)
{
	
	out << "@NAME:TRAITS	VERSION:TYPED	FORMAT:TEXT\n";
	out << "@SAMPLE_NAME	LONG	LAT\n";
	for (size_t x=0; x<N; ++x){
		pointing p=map.pix2ang(x);
		out << t*N+x << '\t' <<  p.theta << '\t' << p.phi << std::endl;
	}
	out << "\n@END_TABLE\n";
}

void print_trait(int N, int t, Individual *ind, std::ostream &out)
{
	out << "@NAME:TRAITS	VERSION:TYPED	FORMAT:TEXT\n";
	out << "@SAMPLE_NAME	TRAIT1\n";
	for (size_t x=0; x<N; ++x){
		out << t*N+x << '\t' <<  ind[x].z << std::endl;
	}
	out << "\n@END_TABLE\n";
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
	for (size_t x=0; x<mut_num; ++x){
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

	const uint32_t MASK_BAD[32]={0x7ffffff,	0x7fffffff,	0x7ffffffd,	0x7ffffff9,
				0x7ffffff1,	0x7fffffe1,	0x7fffffc1,	0x7fffff81,
				0x7fffff01,	0x7ffffe01,	0x7ffffc01, 	0x7ffff801,	
				0x7ffff001, 	0x7fffe001, 	0x7fffc001, 	0x7fff8001,
				0x7fff0001,	0x7ffe0001,	0x7ffc0001,	0x7ff80001,
				0x7ff00001,	0x7fe00001,	0x7fc00001,	0x7f800001,
				0x7f000001,	0x7e000001,	0x7c000001, 	0x78000001,	
				0x70000001, 	0x60000001, 	0x40000001, 	0x00000001};

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
	uint32_t *k;
	double **W;
	std::map <uint32_t, uint32_t> k2x;
	Trait (const size_t &NLOCI, const size_t &NTRAITS, std::mt19937 &mt) {
		k=new uint32_t[NLOCI];
		W=new double *[NTRAITS];

		/*
		std::vector<int> v;
		v.reserve(NLOCI*32);
		int n(0);
		std::generate_n(std::back_inserter(v), NLOCI*32, [n]()mutable { return n++; }); 
		std::random_shuffle(v.begin(), v.end());
		*/
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
		uint8_t Z2=trait.k[k];
		uint32_t x=trait.k2x[k];
		for (uint8_t Z=0; Z<32; Z++)
		{
			if (Z && Z2 != 0)
			{
				uint32_t mask=MASK[Z];
				double W[3];
				std::memcpy(W, trait.W[x], sizeof(double)*3);
				for (size_t i=0; i<N; i++)
				{
					ind[i].z+=W[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
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
	{
		for (size_t y=0; y<generation_size; ++y){
			p+=( (Results[0][y] & MASK[Z])!=0)+( (Results[1][y] & MASK[Z])!=0);
		}
	}
	return p/double(generation_size*2.);
}


/*These were all initially cribbed from Cockerham and Weir 1984.*/
/*But I made some changes. Not sure if they reflect the right thing anymore*/
double get_genome_thingies(uint32_t **ind, 
	const size_t &generation_size, double a, double d)
{

	double p, q, f=0;
	double alpha;
	double delta, z_hat;

	double sigma_A=0., sigma_D=0., sigma_ADI=0., sigma_DI=0., sigma_z=0, f_hat=0, f2=0, N=0;

	for (uint8_t Z=0; Z<LOCI; Z++)
	{
		p=get_freq(ind, generation_size, Z);
		q=1.-p;
		if (p!=0 && q!=0){
			f=1.-get_h(ind, generation_size, Z)/(2.*p*q);
			alpha=a+d*(q-p);
		
			// sigma_A_2
			sigma_A+=2.*p*q*pow(alpha, 2)*(1+f);

			// sigma_D_2
			sigma_D+=4*pow(d,2)*p*q*(p*q+f*pow(1-2*p, 2) );

			// H_2
			sigma_ADI+=4*alpha*d*(pow(p,3)*q-pow(q,3)*p)*f;

			// D1_2
			sigma_DI+=pow(-2*d*p*q*f, 2);
			f_hat+=f;
			f_hat+=f;
			f2+=pow(f, 2);
			N+=1;
			z_hat=(-a)*(pow(p, 2)+p*q*f)+
				(d)*2*p*q*(1-f)+
				(a)*(pow(q, 2)+p*q*f);
		//	std::cerr << z_hat << ", " << -a-z_hat << ", " << d-z_hat << ", " << a-z_hat << std::endl;
			sigma_z+=pow(-a-z_hat, 2)*(pow(p, 2)+p*q*f)+
				 pow(d-z_hat, 2)*2*p*q*(1-f)+
				 pow(a-z_hat, 2)*(pow(q, 2)+p*q*f);
		}
	}
	std::cerr << "s_A, s_D, s_ADI, s_DI, s_z (e), f_hat, f_var" << std::endl;
	std::cerr << sigma_A << ", " << sigma_D << ", " << sigma_ADI << ", " << sigma_DI << ", " << sigma_z << "(" << sigma_A+sigma_D+2*sigma_ADI-sigma_DI << ")" <<  ", " << f_hat/N << ", " << (f2-pow(f_hat, 2)/N)/N<<std::endl;
	return 0;
}


void 
iter (Individual *parents, Individual *children,
	const size_t &generation_size, 
	const Healpix_Map <int> & map, std::mt19937 &mt)
{
	disc_mating(parents, children, generation_size, map, mt);
}

void 
print_ped(Individual *descendents[],
	const size_t &generation_size, 
	const size_t &simulation_length,
	std::ostream &out)
{
	Individual *ind=descendents[0];
	out << "Time" << '\t' << "Individual ID" << '\t' << "Paternal ID" << '\t' << "Maternal ID" << '\t' << "Sex\n"; //<< '\t' <<  "z" << std::endl;
	for (size_t g=0; g<simulation_length; ++g)
	{
		for (size_t x=0; x<generation_size; ++x)
		{
			out << g << '\t' << x+(g+1)*generation_size << '\t' << ind->P1+g*generation_size << '\t' << ind->P2+g*generation_size << '\t' << 3;
			//for (size_t tz=0; tz<TRAITS; ++tz) out << '\t' << ind->z[tz];
			out << std::endl;
			++ind;
		}
	}
}

void
print_results4(std::ostream &out, const int &N, uint32_t **Results)
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
//	out << "-------------------------------" << std::endl;
}

void
print_results_bin(std::ostream &out, const int &N, uint32_t **Results)
{ 
	for (size_t y=0; y<N; ++y){
		out.write((char *)(&(Results[0][y])), sizeof(uint32_t) );
		out.write((char *)(&(Results[1][y])), sizeof(uint32_t) );
	}
//	out << "-------------------------------" << std::endl;
}


/*void
print_results5(std::ostream &out, const int &N, uint_fast32_t **O, uint_fast32_t **P, uint_fast32_t **R)
{ 
	for (size_t b=0; b<32; ++b) {
		for (size_t y=0; y<N; ++y){
			std::bitset <32> O1(O[0][y]);
			std::bitset <32> O2(O[1][y]);

			std::bitset <32> P1(P[0][y]);
			std::bitset <32> P2(P[1][y]);

			std::bitset <32> R1(R[0][y]);
			std::bitset <32> R2(R[1][y]);

			R1[b] ? {p10b='|'; p11b=' ' } : {p10b=' '; p11b='|'};
			R2[b] ? {p20b='|'; p21b=' ' } : {p20b=' '; p21b='|'};

			P1[b] ? p10='0' : p10='1';
			P1[b] ? p11='0' : p11='1';

			P2[b] ? p20='0' : p20='1';
			P2[b] ? p20='0' : p21='1';

			out << '\t' << p10b << p10 << p10b << p11b << p11 << p11b << '\t' << p20b << p20 << p20b << p21b << p21 << p21b;
		}
		out << std::endl;
	}
}*/

void
print_results3(std::ostream &out, const int &k, const int &N, uint32_t **Results)
{ 
for (int K=0; K<k; ++K){
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
print_results2(std::ostream &out, const int &k, const int &N, uint32_t **Results)
{ 
for (int K=0; K<k; ++K){
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
print_results(std::ostream &out, const int &k, const int &N, uint32_t **Results)
{ 
	for (int K=0; K<k; ++K){
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

int 
main(int argc, char *argv[] )
{

	std::mt19937 mt;

	mt.seed(time(NULL));

	if (not(argc==8 ||argc==10 ||  argc==11)) {
		std::cerr << argv[0] << " <Number of Individual> <Number of generations> <Number of loci/8> <r> <pi> <mating_system> <out>\n";
		std::cerr << argv[0] << " <Number of Individual> <Number of generations> <Number of loci/8> <r> <pi> <mating_system> <out> <t1> <t2>\n";
		std::cerr << argv[0] << " <Number of Individual> <Number of generations> <Number of loci/8> <r> <pi> <mating_system> <out> <t1> <t2> <N2>\n";
		return 0;
	}


	size_t N=atoi(argv[1]);
	size_t t=atoi(argv[2]);
	size_t k=atoi(argv[3]);

	double r=atof(argv[4])/double(k*2);	//recombination events per Individual per generation
	double pi=atof(argv[5]);		//equalibrium diversity

	char type=argv[6][0];
	char type2=argv[7][0];
	bool binary=(type2=='b');

	size_t t1;//=atoi(argv[4]);
	size_t t2;//=atoi(argv[5]);
	size_t N2;//=ato

	if (argc==10) {
		t1=atoi(argv[8]);
		t2=atoi(argv[9]);
	}

	if ( argc==11) {
		t1=atoi(argv[8]);
		t2=atoi(argv[9]);
		N2=atoi(argv[10]);
	}

	std::fstream header;
	header.open ("header.txt", std::fstream::out);
	print_head(header, k);
	header.close();

	Individual *parents=new Individual[N];
	Individual *children=new Individual[N];
	
	BGT Pstates(N,k);
	BGT Ostates(N,k);

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
	std::poisson_distribution<int> poisson(double(N)*pi*4.);
#endif
	//std::cerr << argc << std::endl;

//	check(geneology);
	Trait trait(k, 500, mt);

	for (size_t x=0; x<N; ++x)
	{
		Pstate[0][x]=x;//0x00000000;
#ifdef MUT
		Pstate[1][x]=x;//0x00000000;
#else
		Pstate[1][x]=x;//0xFFFFFFFF;
#endif
	};
	for (size_t K=0; K<k; ++K){
		Pstates.compress(Pstate[0], Pstate[1], K);
	}
//	std::cerr << __LINE__ << std::endl; 
	for (size_t y=0; y<t; ++y){
//		#pragma omp parallel for
//		std::cerr << __LINE__ << " y:" << y << std::endl; 
		iter(parents, children, N, map, mt);
#ifdef REC
		initR(Rstate, N, mt);
#endif
		for (size_t K=0; K<k; ++K){
			//std::cerr << __LINE__ << " " << K << std::endl; 
			Pstates.uncompress(Pstate[0], Pstate[1], K);
			for (size_t x=0; x<N; ++x){
#ifdef REC
				children[x].update2( Rstate[0][x], Rstate[1][x], Pstate, Ostate, x, r32, r10, mt);
#else
				children[x].update(mt, Pstate, Ostate, x);
#endif
			}
			inc(children, trait, Ostate, N, K);
#ifdef MUT
			mutate(Ostate[0], r32, rN, poisson, mt);
			mutate(Ostate[1], r32, rN, poisson, mt);
#endif

#ifdef REC
			recombine2((Rstate[0]), (Rstate[0]+N), r32, gcounter, mt);
			recombine2((Rstate[1]), (Rstate[1]+N), r32, gcounter, mt);
#endif
			Ostates.compress(Ostate[0], Ostate[1], K);
		}
		std::swap(Pstates, Ostates);
		std::swap(children, parents);
	}

	for (size_t K=0; K<k; ++K){
		Ostates.uncompress(Results[0], Results[1], K);
		if (binary) print_results_bin(std::cout, N, Results);
		else print_results4(std::cout, N, Results);
	}
	std::fstream pedigree;

	pedigree.open ("name-file.txt", std::fstream::out);
	print_names(N, t, pedigree);
	pedigree.close();

	std::fstream map_file;
	map_file.open ("cord-file.txt", std::fstream::out);
	print_cord(N, t, map, map_file);
	map_file.close();

	map_file.open ("trait-file.txt", std::fstream::out);
	print_trait(N, t, parents, map_file);
	map_file.close();


// TODO WARNING WARNING WARNING!!! These are based off non-causitive loci, so they are at best a guess of the correct results....
//	get_genome_thingies(Results, N, a, d);


	delete [] Pstate[0];
	delete [] Pstate[1];
	delete [] Pstate;
	
	delete [] Ostate[0];
	delete [] Ostate[1];
	delete [] Ostate;

	delete [] Results[0];
	delete [] Results[1];
	delete [] Results;

	return 0;
}
