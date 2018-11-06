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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <ctime>

#include "interface.h"
#include "state.h"
#include "map_file.h"

#include "Eigen/Core"
#include "unsupported/Eigen/MatrixFunctions"

#define REC		//Define if you want linkage between sites, leave undefined for free recombination. 
#define MUT		//Define if you want to introduce mutation each generation.

#define CHECK_FIXED

#define LOCI	32	//Number of loci in a byte. Leave it at 32 unless you are me.
#define BINS	20

#define DEAD	-1

#define NUM_ENV		8
#define NUM_ENV_H	4

#define WORD	32

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
	std::cerr << N << "is an invalid number of Individuals for geometric structure. Try 12*i^2 for any i." << std::endl;
	exit(0);
}

class Individual
{
private:
	static uint32_t ID;
public:
	uint32_t P1, P2;
	double w, z, delta, beta, a, d;
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

void disc_mating(Individual *p, Individual *c, const size_t &size, const double &F, const Healpix_Map <int> &map, std::mt19937 &mt)
{
	std::vector <int> disc;
	int disc_size;

	//A five hundred mile dispersal disc
	//double rad=0.050478006;

	double rad=atan(2*sqrt(F/size) );
	//double rad=std::max(atan(2*sqrt(F/size) ), 1*0.050478006);
	double *zbreak=new double [ int(pow(tan(rad)/2, 2)*size*2) ];
	std::cerr << "size :" << int(pow(tan(rad)/2, 2)*size*2) << std::endl;  
	for (size_t x=0; x<size; ++x)
        {
		map.query_disc( map.pix2ang(x), rad, disc);
		disc_size=disc.size();
		disc_size=disc.size();
		//TODO These calls are a little expensive, cut them?
		zbreak[0]=exp(p[map[disc[0]]].w);
		for (size_t y=1; y<size_t (disc_size); ++y) {
			//TODO These calls are a little expensive, cut them?
			zbreak[y]=zbreak[y-1]+exp(p[map[disc[y]]].w);
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
	zbreak[0]=exp(p[0].w);
	for (size_t y=1; y<size; ++y) {
		zbreak[y]=zbreak[y-1]+exp(p[y].w);
	}	

        for (size_t x=0; x<size; ++x)
        {
                c[x].P1=zfunc(zbreak, size, mt);
                c[x].P2=zfunc(zbreak, size, mt);
        }
	delete [] zbreak;
}


void
rand_2pop(Individual *p, Individual *c, const size_t &size1, const size_t &size2, std::mt19937 &mt)
{

//	std::cerr << "2pop\n";

	double *zbreak1=new double [size1];
	double *zbreak2=new double [size2];

	zbreak1[0]=exp(p[0].w);
	zbreak2[0]=exp(p[size1].w);

	for (size_t y=1; y<size1; ++y) 
	{
		zbreak1[y]=zbreak1[y-1]+exp(p[y].w);
	}

	for (size_t y=size1+1; y<size1+size2; ++y) 
	{
		zbreak2[y-size1]=zbreak2[y-1-size1]+exp(p[y].w);
	}	

        for (size_t x=0; x<size1; ++x)
        {
                c[x].P1=zfunc(zbreak1, size1, mt);
                c[x].P2=zfunc(zbreak1, size1, mt);
//		std::cerr << "(1) "<< x << " " << c[x].P1 <<  ", " << c[x].P2 << " \n";
        }

        for (size_t x=size1; x<size1+size2; ++x)
	{
                c[x].P1=zfunc(zbreak2, size2, mt)+size1;
                c[x].P2=zfunc(zbreak2, size2, mt)+size1;
//		std::cerr << "(2) " << x << " " << c[x].P1 << ", " << c[x].P2 << " \n";
        }

	delete [] zbreak1;
	delete [] zbreak2;
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


void print_trait2(uint32_t N, uint32_t *mask, uint32_t t, Individual *ind, std::ostream &out)
{
	out << "@NAME:TRAITS	VERSION:TYPED	FORMAT:TEXT\n";
	out << "@SAMPLE_NAME	TRAIT1	BETA1	DELTA1	A1	D1\n";
	for (size_t x=0; x<N; ++x){
		if ( mask[x >> 5] & (1 << (x & 0x1F) ) ) 
			out << t*N+x << '\t' << ind[x].z << '\t' <<  ind[x].beta << '\t' <<  ind[x].delta << '\t' << ind[x].a << '\t' << ind[x].d << std::endl;
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

/*
class NTrait
{
private:
public:
	uint32_t *k, *x2k, *x2z;
	double **W;
	uint32_t ntraits;
	std::map <uint32_t, uint32_t> k2x;
	Trait (const size_t &NLOCI, const size_t &NTRAITS,const gsl_vector &mu, const gsl_matrix &L, std::mt19937 &mt) 
	{
		size_t c=mu?

		ntraits=NTRAITS;

		k=new uint32_t[NLOCI];
		memset(k, 0, sizeof(uint32_t)*NLOCI);

		x2k=new uint32_t[NTRAITS];
		x2z=new uint32_t[NTRAITS];

		W=new double *[NTRAITS];

		//CHOOSING LOCI TO BLAH
		std::vector<uint32_t> v(NLOCI*32);
		for (size_t x=0; x<NLOCI*32; x++) v[x]=x;
		std::random_shuffle(v.begin(), v.end());

		//TODO Potential memerror?
		std::vector<uint32_t> ks(v.begin(), v.begin()+NTRAITS);	
		//std::vector<uint32_t> ks(v.begin(), v.begin()+NTRAITS-1);	

		std::sort (ks.begin(), ks.end() );	

		const gsl_rng_type * T;

		gsl_rng *gslr;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		gslr = gsl_rng_alloc(T);

		for (size_t x=0; x<NTRAITS; ++x)
		{
			uint32_t K=ks[x] >> 5;
			uint32_t z=ks[x] & 0x0000001f;

			if ( k2x.find(K) == k2x.end() ) {	
				k2x[K]=x;
			}

			k[K]+=(1 << z);

			x2k[x]=K;
			x2z[x]=z;

			W[x]=new double pow(3,c);

			gsl_vector r;

			gsl_ran_multivariate_gaussian (gslr, mu, const gsl_matrix * L, r);

			for ()
			{
				W[x][?]=r[?];
			}
		}	
		gsl_rng_free (gslr);
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
	
}
*/


class Trait
{
private:
public:
	uint32_t *k, *x2k, *x2z;
	double **W1, **W2;
	uint32_t ntraits;
	std::map <uint32_t, uint32_t> k2x;
	Trait (const size_t &NLOCI, const size_t &NTRAITS, std::mt19937 &mt, const double &mu_a, const double &mu_d, const double &sigma_a, const double &sigma_d, const double &rho) 
	{

		ntraits=NTRAITS;

		k=new uint32_t[NLOCI];
		memset(k, 0, sizeof(uint32_t)*NLOCI);

		x2k=new uint32_t[NTRAITS];	//
		x2z=new uint32_t[NTRAITS];	//

		W1=new double *[NTRAITS];
		W2=new double *[NTRAITS];

		std::vector<uint32_t> v(NLOCI*WORD);
		for (size_t x=0; x<NLOCI*WORD; x++) v[x]=x;
		std::random_shuffle(v.begin(), v.end());

		//TODO Potential memerror?
		std::vector<uint32_t> ks(v.begin(), v.begin()+NTRAITS);	

		std::sort (ks.begin(), ks.end() );	

		const gsl_rng_type * T;
		gsl_rng *gslr;

		gsl_rng_env_setup();

		T = gsl_rng_default;
		gslr = gsl_rng_alloc (T);

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

			gsl_ran_bivariate_gaussian(gslr, sigma_a, sigma_d, rho, &W1[x][2], &W1[x][1]);

			//gsl_ran_bivariate_gaussian(gslr, sigma_a, sigma_d, rho, &W1[x][2], &W1[x][1]);

			W1[x][2]+=mu_a;
			W1[x][1]+=mu_d;
			W1[x][0]=-W1[x][2];
			
			W2[x]=new double[3];
			W2[x][0]=W1[x][0];
			W2[x][1]=W1[x][1];
			W2[x][2]=W1[x][2];
		}	
		gsl_rng_free (gslr);
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


class Interaction
{
private:
public:
	uint32_t *k1, *x2k1, *x2z1;
	uint32_t *k2, *x2k2, *x2z2;

	gsl_matrix **W1;

	uint32_t ntraits;

	std::map <uint32_t, uint32_t> k1_2x;
	std::map <uint32_t, uint32_t> k2_2x;

	Interaction (const size_t &NLOCI, const size_t &NTRAITS, std::mt19937 &mt, const Trait &trait, const Eigen::MatrixXd mu, const Eigen::MatrixXd L) 
	{
		srand(time(0));

		ntraits=NTRAITS;

		k1=new uint32_t[NLOCI];
		memset(k1, 0, sizeof(uint32_t)*NLOCI);

		k2=new uint32_t[NLOCI];
		memset(k2, 0, sizeof(uint32_t)*NLOCI);

		x2k1=new uint32_t[NTRAITS];	//
		x2z1=new uint32_t[NTRAITS];	//

		x2k2=new uint32_t[NTRAITS];	//
		x2z2=new uint32_t[NTRAITS];	//

		W1=new gsl_matrix *[NTRAITS];

		std::vector<uint32_t> v1(NTRAITS);
		for (size_t x=0; x<NTRAITS; x++) v1[x]=x;
		std::random_shuffle(v1.begin(), v1.end());

		std::vector<uint32_t> v2(NTRAITS);
		for (size_t x=0; x<NTRAITS; x++) v2[x]=x;
		std::random_shuffle(v2.begin(), v2.end());

		//TODO Potential memerror?
		std::vector < std::pair<uint32_t, uint32_t> > ks;

		for (int x=0; x<NTRAITS; x++)
		{
			if (v1[x]<v2[x])
				ks.push_back( std::pair<uint32_t, uint32_t>(v1[x], v2[x]) );
			else 
				ks.push_back( std::pair<uint32_t, uint32_t>(v2[x], v1[x]) );
			//ks[x].first  = v1[x];
			//ks[x].second = v2[x];
		}
	
		std::sort (ks.begin(), ks.end() );	

		Eigen::MatrixXd r=Eigen::MatrixXd::Zero(4,1);
		Eigen::MatrixXd X=Eigen::MatrixXd::Zero(4,1);

		Eigen::MatrixXd S00=L.block(0,0,4,4);
		Eigen::MatrixXd S40=L.block(4,0,4,4);
		Eigen::MatrixXd S04=L.block(0,4,4,4);
		Eigen::MatrixXd S44=L.block(4,4,4,4);
		Eigen::MatrixXd iS44 = S44.inverse();

		Eigen::MatrixXd mu_x = mu.block(0,0,4,1);
		Eigen::MatrixXd mu_y = mu.block(4,0,4,1);

		Eigen::MatrixXd SgivenX=S00-S04*iS44*S40;
		Eigen::MatrixXd sqrtmSgivenX = SgivenX.sqrt();

		std::normal_distribution<> norm{0,1};
		std::random_device rd{};
		std::mt19937 gen{rd()};

		for (size_t x=0; x<NTRAITS; ++x)
		{
			uint32_t x1 = ks[x].first;
			uint32_t x2 = ks[x].second;

			uint32_t K1 = trait.x2k[x1];
			uint32_t z1 = trait.x2z[x1];

			uint32_t K2 = trait.x2k[x2];
			uint32_t z2 = trait.x2z[x2];

			//std::cerr << "Trait " << ks[x] << " byte: " << K << ", bit" << z << std::endl;

			if ( k1_2x.find(K1) == k1_2x.end() ) {
				k1_2x[K1]=x;
			}

			if ( k2_2x.find(K2) == k2_2x.end() ) {
				k2_2x[K2]=x;
			}

			k1[K1]+=(1 << z1);
			k2[K2]+=(1 << z2);

			x2k1[x]=K1;
			x2z1[x]=z1;
	
			x2k2[x]=K2;
			x2z2[x]=z2;

			r(0,0) = norm(gen);
			r(1,0) = norm(gen);
			r(2,0) = norm(gen);
			r(3,0) = norm(gen);

			X(0,0)=trait.W1[x1][2];
			X(1,0)=trait.W1[x2][2];
			X(2,0)=trait.W1[x1][1];
			X(3,0)=trait.W1[x2][1];

			Eigen::MatrixXd mu_y_given_x = mu_y+S04*iS44*(X-mu_x);
			Eigen::MatrixXd y=sqrtmSgivenX*r+mu_y_given_x;	

			double MM = y(0, 0);
			double MH = y(1, 0);
			double HM = y(2, 0);
			double HH = y(3, 0);

			std::cerr << ks[x].first << ", " << ks[x].second << ", " << X(0,0) << ", " << X(1,0) << ", " << X(2,0) << ", " << X(3,0) << ", " << MM << ", " << MH << ", " << HM << ", " << HH << std::endl;

			W1[x]=gsl_matrix_alloc(3,3);
								// M1 H1 M2 H2
			gsl_matrix_set(W1[x], 0, 0, 0+MM);  	// -1  0 -1  0
			gsl_matrix_set(W1[x], 1, 0, 0-HM);  	//  0  1 -1  0
			gsl_matrix_set(W1[x], 2, 0, 0-MM);  	//  1  0 -1  0
			gsl_matrix_set(W1[x], 0, 1, 0-MH);  	// -1  0  0  1
			gsl_matrix_set(W1[x], 1, 1, 0+HH);  	//  0  1  0  1
			gsl_matrix_set(W1[x], 2, 1, 0+MH);  	//  1  0  0  1
			gsl_matrix_set(W1[x], 0, 2, 0-MM);  	// -1  0  1  0
			gsl_matrix_set(W1[x], 1, 2, 0+HM);  	//  0  1  1  0
			gsl_matrix_set(W1[x], 2, 2, 0+MM);  	//  1  0  1  0
			
		}
	};

	~Interaction (void)
	{
		//Implement	
	};
	
};
													     //const State
void inc(Individual *ind, Interaction &trait, uint32_t *P[], const size_t &N, const size_t &k, bool noselect, const State &state) 
{
	std::cerr << "INC!\n";
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
	State_stream stream;

	if(trait.k1[k]!=0)
	{
		uint32_t x = trait.k1_2x[k];
		uint32_t k1 = trait.x2k1[x];

		std::cerr << "x " << x << std::endl;
		std::cerr << "k " << trait.k1[k] << std::endl;

		uint32_t **Q=new uint32_t *[2];

		Q[0]=new uint32_t [N];
		Q[1]=new uint32_t [N];

		uint32_t *Q1=Q[0];
		uint32_t *Q2=Q[1];

		while (k1==k)
		{

			uint32_t z1 = trait.x2z1[x];
			uint32_t mask1=MASK[z1];

			std::cerr << "z1 " << z1 << std::endl;
			std::cerr << "mask1 " << mask1 << std::endl;

			uint32_t k2 = trait.x2k2[x];
			uint32_t z2 = trait.x2z2[x];

			uint32_t mask2=MASK[z2];
			
			gsl_matrix *W = gsl_matrix_alloc(3,3);
			gsl_matrix_memcpy(W, trait.W1[x]);

			state.set_stream(stream);
			state.uncompress(Q1, Q2, k2, stream);

			std::cerr << "k1, " << k1 << ", z1, " << z1 << ", k2, " << k2 << ", z2, " <<  z2 << ", ";// << std::endl;
			
			for (size_t i=0; i<N; i++)
			{
				std::cerr <<  ( ( (P1[i] & mask1)!=0)+( (P2[i] & mask1)!=0) ) << ( ( (Q1[i] & mask2)!=0)+( (Q2[i] & mask2)!=0) ) << " ";
				ind[i].z+=gsl_matrix_get(W, ( ( (P1[i] & mask1)!=0)+( (P2[i] & mask1)!=0) ), ( ( (Q1[i] & mask2)!=0)+( (Q2[i] & mask2)!=0) ) );
			}
			std::cerr << std::endl;
			x++;
			k1 = trait.x2k1[x];
		}

		delete Q[0];
		delete Q[1];
		delete Q;
	}

	
	if (!noselect)
	{
		for (size_t i=0; i<N; i++)
			ind[i].w=ind[i].z;
	}
}

void inc(Individual *ind, Trait &trait, uint32_t *P[], const size_t &N, const size_t &k, bool noselect) 
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
				double W1[3];//, W2[3];
				std::memcpy(W1, trait.W1[x], sizeof(double)*3);
				//std::memcpy(W2, trait.W2[x], sizeof(double)*3);
			
				for (size_t i=0; i<N; i++)
					ind[i].z+=W1[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
				x++;
			}
		}	
	}

	
	if (!noselect)
	{
		for (size_t i=0; i<N; i++)
			ind[i].w=ind[i].z;
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

long double get_freq(uint32_t **Results, 
	const size_t &generation_size,
	const uint32_t &Z)
{
	long double p=0;
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
	return p/(long double)(generation_size*2.);
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
	const size_t &generation_size, const char &type, const double &F, 
	const Healpix_Map <int> & map, std::mt19937 &mt, const uint32_t &t, const uint32_t &T)
{
	std::cerr << "Type: " << type << std::endl;
	if (type=='g')
	{
		if (t<T)
		{
			std::cerr << "rand" << std::endl;
			rand_mating(parents, children, generation_size, mt);
		}
		else
		{ 
			std::cerr << "disc" << std::endl;
			disc_mating(parents, children, generation_size, F, map, mt);
		}
	}
	else if (type=='r') 
		rand_mating(parents, children, generation_size, mt);
	else if (type=='a') 
	{
		if (t<T)
			rand_2pop(parents, children, generation_size/4, 3*generation_size/4, mt);
		else
			rand_mating(parents, children, generation_size, mt);
	}
	else 
	{
		std::cerr << "Invalid mating type. Valid options for -y are g, r, a." << std::endl;
		exit(0);
	}
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
push_ped2(Individual *ind,
	const size_t &generation_size, 
	const size_t &time,
	std::ostream &out)
{
	for (size_t x=0; x<generation_size; ++x)
	{
		out << std::setw(10) << time << std::setw(10) << x+(time+1)*generation_size << std::setw(10) << (ind->P1)+time*generation_size << std::setw(10) << (ind->P2)+time*generation_size << std::setw(10) << 3 << std::setw(10) << std::setprecision(2) << ind->z;
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

void 
set_beta(const Trait &traits, State &sites, Individual *ind, uint32_t *mask2, std::ostream &out) 
{
	const uint32_t MASK[LOCI]={0x00000001,	0x00000002,	0x00000004,	0x00000008,
					0x00000010,	0x00000020,	0x00000040,	0x00000080,
					0x00000100,	0x00000200,	0x00000400, 	0x00000800,	
					0x00001000, 	0x00002000, 	0x00004000, 	0x00008000,
					0x00010000,	0x00020000,	0x00040000,	0x00080000,
					0x00100000,	0x00200000,	0x00400000,	0x00800000,
					0x01000000,	0x02000000,	0x04000000, 	0x08000000,	
					0x10000000, 	0x20000000, 	0x40000000, 	0x80000000};

	uint32_t *P[2];
	uint32_t *P1=new uint32_t [sites.sample_size()];
	uint32_t *P2=new uint32_t [sites.sample_size()];

	P[0]=P1;
	P[1]=P2;

	sites.cache();
	sites.rewind();

	size_t l=traits.ntraits;
	size_t n=sites.sample_size();

	/*Eigen::MatrixXf B=Eigen::Zero(l,n);
	Eigen::MatrixXf D=Eigen::Zero(l,n);
	Eigen::MatrixXf C=makev(l);*/

	for (size_t i=0; i<sites.sample_size(); ++i)
	{
		ind[i].beta=0;
		ind[i].a=0;
		ind[i].delta=0;
		ind[i].d=0;
	}

	clock_t t0=clock();
	for (size_t y=0; y<traits.ntraits; y++)
	{
		uint32_t z=traits.x2z[y];
		uint32_t k=traits.x2k[y];
		/* TODO Fix.*/
		long double beta_sum=0, delta_sum=0;
		sites.uncompress(P1, P2, k);
		//std::cerr << k << "::" << z << "-" << P1[0] << std::endl;
		const uint32_t mask=MASK[z];

		const long double a=-traits.W1[y][0];
		const long double d=traits.W1[y][1];
		const long double p=get_freq(P, sites.sample_size(), z);  
		const long double q=1.-p;
		const long double D=2.*p*q-get_h(P, sites.sample_size(), z);  
		const long double alpha=a+d*(q-p);

		long double beta[3]={-p*2.*alpha,(q-p)*alpha,2.*q*alpha};
		long double a_[3]={-a,0,a};
		long double delta[3]={d*(-2.*p*p+D),d*(2.*p*q+D),d*(-2.*q*q+D)};
		long double d_[3]={ 0,d,0};

		for (size_t i=0; i<sites.sample_size(); ++i)
		{
			if ( mask2[i >> 5] & (1 << (i & 0x1F) ) ) 
			{
				ind[i].beta += beta[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
				ind[i].a += a_[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
			}
		}
		for (size_t i=0; i<sites.sample_size(); ++i)
		{
			if ( mask2[i >> 5] & (1 << (i & 0x1F) ) ) 
			{
				ind[i].delta += delta[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
				ind[i].d += d_[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
			}
		}
		if (y%1000==0) 
		{
			clock_t t1=clock();
			std::cerr << "set_beta: " << y << "|" << t1-t0 << " | " << traits.ntraits << std::endl;
			t0=clock();
		}
		
	}

	sites.rewind();

	delete P1;
	delete P2;
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
		if ( mask2[x >> 5] & (1 << (x & 0x1F) ) ) 
			out << '\t' << "beta_" << x;
	}
	for (size_t x=0; x<sites.sample_size() ; ++x){
		if ( mask2[x >> 5] & (1 << (x & 0x1F) ) ) 
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
		long double beta_sum=0, delta_sum=0;
		sites.uncompress(P1, P2, k);
		//std::cerr << k << "::" << z << "-" << P1[0] << std::endl;
		const uint32_t mask=MASK[z];

		const long double a=-traits.W1[y][0];
		const long double d=traits.W1[y][1];
		const long double p=get_freq(P, sites.sample_size(), z);  
		const long double q=1.-p;
		const long double D=2.*p*q-get_h(P, sites.sample_size(), z);  
		const long double alpha=a+d*(q-p);

		long double beta[3]={-p*2.*alpha,(q-p)*alpha,2.*q*alpha};
		long double delta[3]={d*(-2.*p*p+D),d*(2.*p*q+D),d*(-2.*q*q+D)};
	

		for (size_t i=0; i<sites.sample_size(); ++i){
			if ( mask2[i >> 5] & (1 << (i & 0x1F) ) ) 
				{
				out << '\t' << beta[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
				//std::cerr << int ( ((P1[i] & mask)!=0)+((P2[i] & mask)!=0) ) << '\t';
				//beta_sum+=beta[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
				}
		}
		for (size_t i=0; i<sites.sample_size(); ++i){
			if ( mask2[i >> 5] & (1 << (i & 0x1F) ) ) 
				{
				out << '\t' << delta[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
				//delta_sum+=delta[((P1[i] & mask)!=0)+((P2[i] & mask)!=0)];
				}
		}
		if (y%1000==0) std::cerr << y << "::" << traits.ntraits << std::endl;
		out << "\n";
	}
	sites.rewind();
	delete P1;
	delete P2;
}

/*print_cov is bad?*/
void 
print_cov2(const Trait &traits, State &sites, uint32_t *mask2, std::ostream &out) 
{
	std::cerr << "print_cov2:" << std::endl;
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
	out << '\t' << "a" << '\t' << "d";
	for (size_t x=0; x<sites.sample_size() ; ++x){
		if ( mask2[x >> 5] & (1 << (x & 0x1F) ) ) 
			out << '\t' << "m_" << x;
	}
	for (size_t x=0; x<sites.sample_size() ; ++x){
		if ( mask2[x >> 5] & (1 << (x & 0x1F) ) ) 
			out << '\t' << "h_" << x;
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
		sites.uncompress(P1, P2, k);
		if (y%1000==0) std::cerr << y << "::" << traits.ntraits << ", " << k << std::endl;
		const uint32_t mask=MASK[z];

		const long double a=-traits.W1[y][0];
		const long double d=traits.W1[y][1];

		out << '\t' << a << '\t' << d; 

		for (size_t i=0; i<sites.sample_size(); ++i){
			if ( mask2[i >> 5] & (1 << (i & 0x1F) ) ) 
				{
				out << '\t' << ((P1[i] & mask)!=0)+((P2[i] & mask)!=0);
				}
		}
		for (size_t i=0; i<sites.sample_size(); ++i){
			if ( mask2[i >> 5] & (1 << (i & 0x1F) ) ) 
				{
				out << '\t' << int(((P1[i] & mask)!=0)!=((P2[i] & mask)!=0));
				}
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

	bool binary=false, p_pedigree=false, p_names=false, p_traits=false, p_traits2=false, p_geography=false, p_dist=false, p_header=false, p_matrix2=false, p_matrix=false, noselect=false, noprint=false;
	bool p_stats=false;

	int N=1200, T=0, t=10, k=500, n_trait=50, skip=100, sub_sample_size=-1;

	uint32_t *mask;

	double r=36;			//recombination events per Individual per generation
	double pi=0.001;		//equilibrium diversity
	double var=0.01;		//trait var
	double theta=1;			//equilibrium diversity

	double mu_d=1;			//?
	double mu_a=1;			//?
	double sigma_d=1;		//?
	double sigma_a=1;		//?
	double rho=0;			//?
	double F=16;			//?

	char type='g';			//geographic or random?

	double sigma_i2=0.01, sigma_j2=0.01, sigma_k2=0.01;
	double sigma_ij=0.001, sigma_ik=0.001, sigma_jk=0.001, sigma_ai=0.001, sigma_aj=0.001;
	double sigma_ak=0.001, sigma_di=0.001, sigma_dj=0.001, sigma_dk=0.001;

	env.set_name("pedigree_sim");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("A MPI capable forward time evolution simulation. Please direct questions to matthew.s.ackerman@gmail.com");

	env.optional_arg('N',"number", 	N,	"please provide a number.", "number of individuals in the populations.");
	env.optional_arg('S',"sigma", 	var,	"please provide a number.", "variance of trait per snp");
	env.optional_arg('T',"time", 	T,	"please provide a number.", "time of switch.");
	env.optional_arg('a',"theta", 	theta,	"please provide a number.", "angle separating environments");
	env.optional_arg('c',"chiasma",	r,	"please provide a number.", "number of recombinations events per individual per generation.");
	env.optional_arg('e',"equil", 	pi,	"please provide a number.", "equilibrium nucleotide diversity.");
	env.optional_arg('f',"inbred", 	F,	"please provide a number.", "stable state inbreeding.");
	env.optional_arg('g',"gen", 	t,	"please provide a number.", "number of generations the simulation runs.");
	env.optional_arg('i',"init", 	i_file,	"please provide a number.", "state file to load");
	env.optional_arg('k',"skip", 	skip,	"please provide a number.", "number of generations to skip between printing files");
	env.optional_arg('s',"sites", 	k,	"please provide a number.", "number of sites/32.");
	env.optional_arg('v',"var", 	n_trait,"please provide a number.", "number of loci effecting the trait.");
	env.optional_arg('y',"type", 	type,	"please provide a number.", "type of mating < r | a | g > for random, admixed or geographic.");
	env.optional_arg(NULL, "mu_d", 	mu_d,	"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "mu_a", 	mu_a,	"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_d", sigma_d,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_a", sigma_a,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "rho", rho,"please provide a number.", "undocumented mode.");

	env.optional_arg(NULL, "sigma_ai", sigma_ai,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_aj", sigma_aj,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_ak", sigma_ak,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_di", sigma_di,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_dj", sigma_dj,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_dk", sigma_dk,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_i2", sigma_i2,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_ij", sigma_ij,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_ik", sigma_ik,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_j2", sigma_j2,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_jk", sigma_jk,"please provide a number.", "undocumented mode.");
	env.optional_arg(NULL, "sigma_k2", sigma_k2,"please provide a number.", "undocumented mode.");

	env.optional_arg('1',"sample", 	sub_sample_size,	"please provide a number.", "sub sample.");

	env.flag(	'G',"geography",&p_geography,&flag_set, 	"an error occurred while displaying the version message", "prints the geographic coordinates of each individual.");		//DONE
	env.flag(	'H',"header", &p_header,	&flag_set, 	"an error occurred while displaying the version message", "prints a header file.");		//DONE
	env.flag(	'M',"matrix2", &p_matrix2,	&flag_set, 	"an error occurred while displaying the version message", "prints a covariance matrix of additive and dominance actions.");		//DONE
	env.flag(	'O',"opts", 	&env, 		&flag_options, "an error occurred while displaying the help message", "Prints a list of available options");		//DONE
	env.flag(	'b',"binary", 	&binary, 	&flag_set, 	"an error occurred while displaying the version message", "binary output");		//DONE
	env.flag(	'd',"dist", &p_dist,&flag_set, 	"an error occurred while displaying the version message", "prints the binned frequencies of a and d [Sorry guys, I'll figure out a way to explain this later.]");		//DONE
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message", "Prints this message");				//DONE
	env.flag(	'l',"lynch", &noselect,	&flag_set, 	"an error occurred while displaying the version message", "trait has no fitness effects.");		//DONE
	env.flag(	'm',"matrix", &p_matrix,	&flag_set, 	"an error occurred while displaying the version message", "prints a covariance matrix of additive and dominance deviations.");		//DONE
//	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occurred while displaying the version message", "Prints the program version");		//DONE
	env.flag(	'p',"pedigree", &p_pedigree,&flag_set, 	"an error occurred while displaying the version message", "prints the pedigree in the file 'pedigree.txt' ");		//DONE
	env.flag(	'n',"name", &p_names,	&flag_set, 	"an error occurred while displaying the version message", "prints the names of the samples.");		//DONE
	env.flag(	'o',"off", &noprint,	&flag_set, 	"an error occurred while displaying the version message", "don't print state file.");		//DONE
	env.flag(	't',"traits", &p_traits,	&flag_set, 	"an error occurred while displaying the version message", "prints the trait values every generation in files called t[GENERATION].txt");		//DONE
	env.flag(	'u',"traits2", &p_traits2,	&flag_set, 	"an error occurred while displaying the version message", "prints the trait values every generation in files called t[GENERATION].txt");		//DONE
	env.flag(	'w',"window", &p_stats,	&flag_set, 	"an error occurred while displaying the version message", "print stats file.");		//DONE

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);

	std::mt19937 mt;

	//mt.seed(time(NULL));
	mt.seed(28);

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

	mask=new uint32_t [(N >> 5)+1];
	memset(mask, 0xFFFFFFFF, sizeof(uint32_t)*( (N >> 5)+1) );

	if (sub_sample_size>N)
	{
		std::cerr << "Sub sample is too large, reducing to " << N << ".\n";
		sub_sample_size==-1;
	}

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

	if (n_trait>32*k)
	{
		std::cerr << "More trait effecting loci than total loci, reducing -v\n";
		n_trait=k*32-1;
	}

	State Ostates(N);

	r=r/double(2*k);

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

//	check(geneology);
	Trait trait(k, n_trait, mt, mu_a, mu_d, sigma_a, sigma_d, rho);

	double sigma_ad=sigma_a*sigma_d*rho;
	double sigma_a2=sigma_a*sigma_a, sigma_d2=sigma_d*sigma_d;
	double mu_i=0, mu_j=0, mu_k=0;
#ifdef INTERACT
	Eigen::MatrixXd mu = Eigen::MatrixXd::Zero(8,1);
	Eigen::MatrixXd C  = Eigen::MatrixXd::Zero(8,8);

	mu <<   mu_a, mu_a, mu_d, mu_d, mu_i, mu_j, mu_j, mu_k;

	C  << 	sigma_a2 , 0. 	     , sigma_ad , 0. 	     , sigma_ai , sigma_aj , sigma_aj , sigma_ak ,
		0. 	 , sigma_a2 , 0.       , sigma_ad , sigma_ai , sigma_aj , sigma_aj , sigma_ak ,
		sigma_ad , 0.       , sigma_d2 , 0. 	     , sigma_di , sigma_dj , sigma_dj , sigma_dk ,
		0.	 , sigma_ad , 0.       , sigma_d2 , sigma_di , sigma_dj , sigma_dj , sigma_dk ,
		sigma_ai , sigma_ai , sigma_di , sigma_di , sigma_i2 , sigma_ij , sigma_ij , sigma_ik ,
		sigma_aj , sigma_aj , sigma_dj , sigma_dj , sigma_ij , sigma_j2 , 0.       , sigma_jk ,
		sigma_aj , sigma_aj , sigma_dj , sigma_dj , sigma_ij , 0.       , sigma_j2 , sigma_jk ,
		sigma_ak , sigma_ak , sigma_dk , sigma_dk , sigma_ik , sigma_jk , sigma_jk , sigma_k2;

	Interaction interaction(k, n_trait, mt, trait, mu, C);
#endif
	
#ifndef MUT
	std::bernoulli_distribution bernoulli(0.2);
#endif

	if (i_file=="")
	{
		for (size_t K=0; K<k; ++K){
#ifndef MUT
			uint32_t set=0;
			for (size_t bit=0; bit < 32; bit++)
			{
				if ( bernoulli(mt) )
					set |= 1UL << bit;
			}
 
#endif
			for (size_t x=0; x<N; ++x)
			{
#ifdef MUT
				Pstate[0][x]=0x00000000;
				Pstate[1][x]=0x00000000;
#else
				if (x<N/2) 
				{
					Pstate[0][x]=0x00000000;
//					Pstate[1][x]=set;
					Pstate[1][x]=0xFFFFFFFF;
				} 
				else 
				{
					Pstate[0][x]=0x00000000;
//					Pstate[1][x]=set;
					Pstate[1][x]=0xFFFFFFFF;
				}
#endif
			}
			Pstates.compress(Pstate[0], Pstate[1]);
		}
		Pstates.cache();
	} 

	{

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

	for (size_t x=0; x<N; x++)
	{
		parents[x].z=0;
		parents[x].w=0;
	}

	for (size_t y=0; y<t; ++y){
//		#pragma omp parallel for

		std::cerr << __LINE__ << " y:" << y << std::endl; 
		iter(parents, children, N, type, F, map, mt, y, T);

		for (size_t x=0; x<N; x++)
		{
			children[x].z=0;
			children[x].w=0;
		}

#ifdef REC
		initR(Rstate, N, mt);
#endif

#ifdef CHECK_FIXED
		uint32_t fixed=0xFFFFFFFF;
#endif
		for (size_t K=0; K<k; ++K)
		{
#ifdef CHECK_FIXED
			fixed=0xFFFFFFFF;
#endif

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
			if (!noselect || y%skip==skip-1 || y==t-1 )
			{
				inc(children, trait, Ostate, N, K, noselect);
#ifdef INTERACT
				inc(children, interaction, Ostate, N, K, noselect, Pstates);
#endif
			}

#ifdef REC
			recombine2((Rstate[0]), (Rstate[0]+N), r32, gcounter, mt);
			recombine2((Rstate[1]), (Rstate[1]+N), r32, gcounter, mt);
#endif

#ifdef CHECK_FIXED
			for (size_t x=0; x<N; ++x)
			{
				fixed &= Ostate[0][x];
				fixed &= Ostate[1][x];
			}
			if (fixed != 0)
			{
				fixed=~fixed;
				for (size_t x=0; x<N; ++x)
				{
					Ostate[0][x] &= fixed;
					Ostate[1][x] &= fixed;
				}
				fixed=~fixed;
			}
#endif 

//			check fixation
//			std::cerr << "Fixations ??? " << std::endl;

			Ostates.compress(Ostate[0], Ostate[1]);
		}
			
	
		Ostates.cache();
		std::cerr << "free :" << Ostates.get_free() << " ratio " << Ostates.compression_ratio() << std::endl;
		Pstates.clear();

		Pstates.swap(Ostates);
		std::swap(children, parents);
	

		if(p_pedigree){
			if (t >= 10 )
			{
				if(y > t-10)
					push_ped2(parents, N, y, pedigree);
			}
			else 
				push_ped2(parents, N, y, pedigree);
		}
		
		if (y%skip==skip-1)
		{
		if(p_traits){
			std::cerr << __LINE__ << std::endl;
			std::fstream trait_file;
			trait_file.open("t"+std::to_string(y)+".txt", std::fstream::out);
			print_trait(N, mask, t, parents, trait_file);
			trait_file.close();
		}

		if(p_dist){
			std::cerr << __LINE__ << std::endl;
			std::fstream bin_file;
			bin_file.open("b"+std::to_string(y)+".txt", std::fstream::out);
			print_bins(trait, Pstates, bin_file);
			bin_file.close();
		}

		if(p_stats)
		{
			std::cerr << __LINE__ << std::endl;
			std::fstream stats_file;
			stats_file.open("s"+std::to_string(y)+".txt", std::fstream::out);
			print_stats(trait, Pstates, parents, stats_file);
			stats_file.close();
		}
		}
	}

	if(p_pedigree)
	{
		std::cerr << __LINE__ << std::endl;
		pedigree.close();
	}
	
	if (p_names)
	{
		std::cerr << __LINE__ << std::endl;
		std::fstream names;
		names.open ("name-file.txt", std::fstream::out);
		print_names(N, mask, t, names);
		names.close();
	}

	if(p_traits)
	{
		std::cerr << __LINE__ << std::endl;
		std::fstream trait_file;
		trait_file.open("t_final.txt", std::fstream::out);
		print_trait(N, mask, t, parents, trait_file);
		trait_file.close();
	}

	if(p_traits2)
	{
		std::cerr << __LINE__ << std::endl;
		std::fstream trait_file;
		trait_file.open("beta_delta.txt", std::fstream::out);
		set_beta(trait, Pstates, parents, mask, trait_file);
		print_trait2(N, mask, t, parents, trait_file);
		trait_file.close();
	}

	if(p_geography)
	{
		std::cerr << __LINE__ << std::endl;
		std::fstream map_file;
		map_file.open ("cord-file.txt", std::fstream::out);
		print_cord(N, mask, t, map, map_file);
		map_file.close();
	}

	if(p_matrix)
	{
		std::cerr << __LINE__ << std::endl;
		std::fstream cov_file;
		cov_file.open("cov.txt", std::fstream::out);
		print_cov(trait, Pstates, mask, cov_file);
		cov_file.close();
	}

	if(p_matrix2)
	{
		std::cerr << __LINE__ << std::endl;
		std::fstream cov_file;
		cov_file.open("cov2.txt", std::fstream::out);
		print_cov2(trait, Pstates, mask, cov_file);
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
//		State sub_states=rand_drop(sub_sample(Pstates, sub_sample_size, mask), 0.25, mt);
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

	delete [] mask;

#ifdef REC
	delete [] Rstate[0];
	delete [] Rstate[1];
	delete [] Rstate;
#endif 
	return 0;
}
