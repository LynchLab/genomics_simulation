#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <bitset>
#include <time.h>
#include <cstring>
#include <omp.h>

#define REC

#define LOCI	32

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

void print_head(std::ostream &out, const int &length)
{
	out << "@HD	VN:0.4	SO:coordinate\n@SQ	SN:scaffold_1	LN:" << length*32 << "\n";
}

void print_names(int N, int t, std::ostream &out)
{
	
	out << "@NAME:FILES	VERSION:TYPED	FORMAT:TEXT\n";
	out << "@FILE_NAME	SAMPLE_NAME	...\n";
	out << "mpileup.txt";
	for (size_t x; x<N; ++x){
		out << '\t' << t*N+x;
	}
	out << "\n@END_TABLE\n";
}

/*
void print_ld(std::ostream &out, const int &length)
{
	out;
}
*/

/*
void print_relatedness(std::ostream &out, const ..)
{
	
}
*/

class individual
{
private:
	static uint32_t ID;
public:
	uint32_t P1, P2;
	double z;
//	uint32_t id;
//	uint32_t P1_id;
//	uint32_t P2_id;

/*
	individual(const uint32_t &p1, const uint32_t &p2) {
		P1 = p1;
		P2 = p2;
		z=1.0;
		id=++std::uniform_int_distribution<int> ri	P1_id=
		P2_id=
	};*/

	individual() {
		P1=-1;
		P2=-1;
		z=0.0;
//		id=++ID;
//		P1_id=0;
//		P2_id=0;
	};

	void
	update(std::mt19937 &mt, uint_fast32_t *Pstate[], uint_fast32_t *Ostate[], uint32_t x) {
		uint_fast32_t r1 = mt();
		uint_fast32_t r2 = mt();
		*(Ostate[0]+x) = *(Pstate[0]+P1) & r1 | *(Pstate[1]+P1) & ~r1; 
		*(Ostate[1]+x) = *(Pstate[0]+P2) & r2 | *(Pstate[1]+P2) & ~r2;
	};

	void
	update2(const uint_fast32_t &r1, const uint_fast32_t &r2, uint_fast32_t *Pstate[], uint_fast32_t *Ostate[], uint32_t x) {
		*(Ostate[0]+x) = *(Pstate[0]+P1) & r1 | *(Pstate[1]+P1) & ~r1; 
		*(Ostate[1]+x) = *(Pstate[0]+P2) & r2 | *(Pstate[1]+P2) & ~r2;
	}

};

uint32_t individual::ID=1;



void
recombine (uint_fast32_t &r, std::uniform_int_distribution<int> &ri, std::mt19937 &mt)
{
	//Switch based on the value of the lowest bit
	const uint_fast32_t MASK[32]={0x00000001,	0x00000003,	0x00000007,	0x0000000F,
					0x0000001F,	0x0000003F,	0x0000007F,	0x000000FF,
					0x000001FF,	0x000003FF,	0x000007FF, 	0x00000FFF,	
					0x00001FFF, 	0x00003FFF, 	0x00007FFF, 	0x0000FFFF,
					0x0001FFFF,	0x0003FFFF,	0x0007FFFF,	0x000FFFFF,
					0x001FFFFF,	0x003FFFFF,	0x007FFFFF,	0x00FFFFFF,
					0x01FFFFFF,	0x03FFFFFF,	0x07FFFFFF, 	0x0FFFFFFF,	
					0x1FFFFFFF, 	0x3FFFFFFF, 	0x7FFFFFFF, 	0xFFFFFFFF};
	if (r & 0x00000001) r=~MASK[ri(mt)];
	else r=MASK[ri(mt)];
}

void
rand_rec(uint_fast32_t &r, std::uniform_int_distribution<int> &r32, std::uniform_int_distribution<int> &r10, std::mt19937 &mt)
{
	if (!r10(mt)){
		recombine(r, r32, mt);
	} else {
		if (r & 0x00000001) r=0xFFFFFFFF;
		else r=0;	
	}
}

void 
initR (uint_fast32_t ***Rand,const size_t & N, const size_t & t, std::mt19937 &mt)
{
	uint_fast32_t *r1_ptr, *r2_ptr, *r1_end;
	uint_fast32_t ***r_out=Rand;
	uint_fast32_t ***rout_end=Rand+t;
	
	std::uniform_int_distribution<int> r32(0,32);
	std::uniform_int_distribution<int> r10(0,10);


	while (r_out!=rout_end)
	{
		r1_ptr=(*r_out)[0];
		r2_ptr=(*r_out)[1];
		r1_end=(*r_out)[0]+N;
		while (r1_ptr!=r1_end)
		{
//			std::cerr << "HI!\n";
			*r1_ptr=mt();
			*r2_ptr=mt();
			rand_rec(*r1_ptr, r32, r10, mt);
			rand_rec(*r2_ptr, r32, r10, mt);
			++r1_ptr;
			++r2_ptr;
		}	
//		std::cerr << "1\n";
		++r_out;
	} 
//	std::cerr << "BYE!\n";
}
/*These are all cribbed from Cockerham and Weir 1984.*/

double get_freq(uint_fast32_t **Results, 
	const size_t &generation_size,
	const uint32_t &Z)
{
	double p=0;
	const uint_fast32_t MASK[LOCI]={0x00000001,	0x00000002,	0x00000004,	0x00000008,
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


double get_genome_thingies(uint_fast32_t **ind, 
	const size_t &generation_size, double a, double d)
{

	double p, q, f=0;
	double alpha;
	double delta;

	double sigma_A=0., sigma_D=0., sigma_ADI=0., sigma_DI=0.;

	for (uint8_t Z=0; Z<LOCI; Z++)
	{
		p=get_freq(ind, generation_size, Z);
		q=1.-p;
		if (p!=0){
			alpha=a+d*(q-p);
//			alpha=q*(a+d*(q-p) );
			delta=2*p*q*d*f;
		
			// sigma_A_2
			sigma_A+=2.*p*q*pow(alpha, 2);

			// sigma_D_2
			sigma_D+=pow(2*p*q*d, 2);

			// H_2
			sigma_ADI+=4*alpha*d*(pow(p,3)*q-pow(q,3)*p);

			// D1_2
			sigma_DI+=4*pow(alpha,2)*(pow(q,2)*p+pow(p,2)*q);
		}
	}
	std::cerr << sigma_A << ", " << sigma_D << ", " << sigma_ADI << ", " << sigma_DI << std::endl;
	return 0;
}

double
inc(double &z, 
	uint_fast32_t &P1, 
	uint_fast32_t &P2,
	double a,
	double d)
{

	double W[3]={-a, d, a};

	const uint_fast32_t MASK[LOCI]={0x00000001,	0x00000002,	0x00000004,	0x00000008,
					0x00000010,	0x00000020,	0x00000040,	0x00000080,
					0x00000100,	0x00000200,	0x00000400, 	0x00000800,	
					0x00001000, 	0x00002000, 	0x00004000, 	0x00008000,
					0x00010000,	0x00020000,	0x00040000,	0x00080000,
					0x00100000,	0x00200000,	0x00400000,	0x00800000,
					0x01000000,	0x02000000,	0x04000000, 	0x08000000,	
					0x10000000, 	0x20000000, 	0x40000000, 	0x80000000};
	for (uint8_t Z=0; Z<LOCI; Z++)
	{
		z+=W[((P1 & MASK[Z])!=0)+((P2 & MASK[Z])!=0)];
	}
	return z;
}

void 
make (individual *descendents[], 
	const size_t &generation_size, 
	const size_t &simulation_length,
	std::mt19937 &mt)
{
	std::uniform_int_distribution<int> ri(0,generation_size-1);
	*descendents=new individual[generation_size*simulation_length];

	uint32_t *P1=new uint32_t[generation_size], *P2= new uint32_t[generation_size];

//	for (size_t x=0; x<generation_size;*simulation_length; ++x)
	for (size_t y=0; y<simulation_length; ++y)
	{
		random_sibs(P1, P2, generation_size);
		for (size_t x=0; x<generation_size; ++x)
		{
			(*descendents)[y*generation_size+x].P1=P1[x];
			(*descendents)[y*generation_size+x].P2=P2[x];
		}
	}
	delete P1;
	delete P2;
}

void 
make_subdivided (individual *descendents[], 
	const size_t &generation_size, 
	const size_t &gh,
	const size_t &simulation_length,
	const size_t &start,
	const size_t &stop,
	std::mt19937 &mt)
{
	std::uniform_int_distribution<int> ri(0,  generation_size-1);
	std::uniform_int_distribution<int> ri1(0,  gh-1);
	std::uniform_int_distribution<int> ri2(gh, generation_size-1);

	*descendents=new individual[generation_size*simulation_length];

	int P1, P2, gen=0;
	size_t i;

	for (size_t t=0; t<simulation_length; ++t)
	{
		for (size_t x=0; x<generation_size; ++x)
		{
			i=t*generation_size+x;
			if (start < t && t < stop){
				if (x < gh){
					P1=ri1(mt);
					P2=ri1(mt);
					(*descendents)[i].P1=P1;
					(*descendents)[i].P2=P2;
				} else {
					P1=ri2(mt);
					P2=ri2(mt);
					(*descendents)[i].P1=P1;
					(*descendents)[i].P2=P2;
				}
			} else {
				P1=ri(mt);
				P2=ri(mt);
				(*descendents)[i].P1=P1;
				(*descendents)[i].P2=P2;
			}
		}
	}
}
void 
print_ped(individual *descendents[],
	const size_t &generation_size, 
	const size_t &simulation_length,
	std::ostream &out)
{
	individual *ind=descendents[0];
	out << "Time" << '\t' << "Individual ID" << '\t' << "Paternal ID" << '\t' << "Maternal ID" << '\t' << "Sex" << '\t' <<  "z" << std::endl;
	for (size_t g=0; g<simulation_length; ++g)
	{
		for (size_t x=0; x<generation_size; ++x)
		{
			out << g << '\t' << x+(g+1)*generation_size << '\t' << ind->P1+g*generation_size << '\t' << ind->P2+g*generation_size << '\t' << 3 << '\t' << ind->z << std::endl;
			++ind;
		}
	}
}

void
print_results4(std::ostream &out, const int &N, uint_fast32_t **Results)
{ 
	for (size_t b=0; b<32; ++b) {
		for (size_t y=0; y<N; ++y){
			std::bitset <32> bits0(Results[0][y]);
			std::bitset <32> bits1(Results[1][y]);
			if (bits0[b]) {
				if (bits1[b]) {
					out << "	0	0";
				} else {
					out << "	0	1";
				}
			} else {
				if (bits1[b]) {
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
print_results3(std::ostream &out, const int &k, const int &N, uint_fast32_t **Results)
{ 
for (int K=0; K<k; ++K){
	for (size_t b=0; b<32; ++b) {
		for (size_t y=0; y<N; ++y){
			std::bitset <32> bits0(Results[0][N*K+y]);
			std::bitset <32> bits1(Results[1][N*K+y]);
			if (bits0[b]) {
				if (bits1[b]) {
					out << "	0	0";
				} else {
					out << "	0	1";
				}
			} else {
				if (bits1[b]) {
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
print_results2(std::ostream &out, const int &k, const int &N, uint_fast32_t **Results)
{ 
for (int K=0; K<k; ++K){
	for (size_t b=0; b<32; ++b) {
		out << "scaffold_1	" << b+K*32+1;
		for (size_t y=0; y<N; ++y){
			std::bitset <32> bits0(Results[0][N*K+y]);
			std::bitset <32> bits1(Results[1][N*K+y]);
			if (bits0[b]) {
				if (bits1[b]) {
					out << "	0";
				} else {
					out << "	1";
				}
			} else {
				if (bits1[b]) {
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
print_results(std::ostream &out, const int &k, const int &N, uint_fast32_t **Results)
{ 
	for (int K=0; K<k; ++K){
		for (size_t b=0; b<32; ++b) {
			out << "scaffold_1	" << b+K*32+1 << "	G";
			for (size_t y=0; y<N; ++y){
				std::bitset <32> bits0(Results[0][N*K+y]);
				std::bitset <32> bits1(Results[1][N*K+y]);
				if (bits0[b]) {
					if (bits1[b]) {
						out << "	12	TtTtTtTtTtTt	~~~~~~~~~~~";
					} else {
						out << "	12	TtTtTtAaAaAa	~~~~~~~~~~~";
					}
				} else {
					if (bits1[b]) {
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

	if (not(argc==6 || argc==9)) {
		std::cerr << argv[0] << " <Number of individual> <Number of generations> <Number of loci/8> <a> <d>\n";
		std::cerr << argv[0] << " <Number of individual> <Number of generations> <Number of loci/8> <a> <d> <t1> <t2> <N2>\n";
		return 0;
	}


	size_t N=atoi(argv[1]);
	size_t t=atoi(argv[2]);
	size_t k=atoi(argv[3]);

	size_t a=atof(argv[4]);
	size_t d=atof(argv[5]);

	size_t t1;//=atoi(argv[4]);
	size_t t2;//=atoi(argv[5]);
	size_t N2;//=atoi(argv[6]);

	if ( argc==9) {
		t1=atoi(argv[6]);
		t2=atoi(argv[7]);
		N2=atoi(argv[8]);
	}


	std::fstream header;
	header.open ("header.txt", std::fstream::out);
	print_head(header, k);
	header.close();

	individual *descendents;
	individual *this_generation;
		
	uint_fast32_t **Pstate, **Ostate, **Results;
	uint_fast32_t ***Rstate, ***thisR;	

	uint_fast32_t **Recombination;
//	std::bitset <32> br1(Results[0][N*K+y]);

	Pstate=new uint_fast32_t *[2];
	Pstate[0]=new uint_fast32_t [N];
	Pstate[1]=new uint_fast32_t [N];

	Ostate=new uint_fast32_t *[2];
	Ostate[0]=new uint_fast32_t[N];
	Ostate[1]=new uint_fast32_t[N];

	Results=new uint_fast32_t *[2];
//	Results[0]=new uint_fast32_t[N*k];
//	Results[1]=new uint_fast32_t[N*k];
	Results[0]=new uint_fast32_t[N];
	Results[1]=new uint_fast32_t[N];

	// [0]->N [1]->N ](t) **t,2,N
#ifdef REC
	Rstate=new uint_fast32_t **[t];
	for (size_t y=0; y<t; y++)
	{
		Rstate[y]=new uint_fast32_t *[2];
		Rstate[y][0]=new uint_fast32_t[N];
		Rstate[y][1]=new uint_fast32_t[N];
	}
#endif

	if (argc==9) make_subdivided(&descendents, N, N2, t, t1, t2, mt);
	else make(&descendents, N, t, mt);

#ifdef REC
	initR(Rstate, N, t, mt);
#endif
//	check(geneology);

	for (size_t K=0; K<k; ++K){	
		this_generation=descendents;

		for (size_t x=0; x<N; ++x)
		{
			Pstate[0][x]=0x00000000;
			Pstate[1][x]=0xFFFFFFFF;
		};

		thisR=Rstate;

		for (size_t y=0; y<t; ++y){
			#pragma omp parallel for
			for (size_t x=0; x<N; ++x){
				//array of xs and bp.
#ifdef REC
				this_generation[x].update2(*((*thisR)[0]+x), *((*thisR)[1]+x), Pstate, Ostate, x);
#else
				this_generation[x].update(mt, Pstate, Ostate, x);
#endif
				if (K<1) this_generation[x].z=inc(this_generation[x].z,*(Ostate[0]+x),*(Ostate[1]+x), a, d);
			}
		//	get_genome_thingies(Ostate, N);
			std::swap(Pstate,Ostate);
			this_generation+=N;
#ifdef REC
			++thisR;
#endif
		}
//		memcpy(Results[0]+N*K, Pstate[0], N*8);
//		memcpy(Results[1]+N*K, Pstate[1], N*8);
		memcpy(Results[0], Pstate[0], N*8);
		memcpy(Results[1], Pstate[1], N*8);
		print_results4(std::cout, N, Results);
	}

//	print_results(std::cout, k, N, Results);
//	print_results2(std::cout, k, N, Results);
//	print_results3(std::cout, k, N, Results);

	std::fstream pedigree;
	pedigree.open ("pedigree.txt", std::fstream::out);
	print_ped(&descendents, N, t, pedigree);
	pedigree.close();

	pedigree.open ("name-file.txt", std::fstream::out);
	print_names(N, t, pedigree);
	pedigree.close();

	this_generation-=N;
	
	get_genome_thingies(Results, N, a, d);
	
	delete [] descendents;

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
