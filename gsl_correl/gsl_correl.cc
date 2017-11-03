		#include "circular_list.h"
#include "interface.h"
#include "map_file.h"
#include "state.h"
#include "sample_name.h"
#include "correl_data.h"
#include "triu_index.h"

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cstring>
#include <sstream>
#include <tuple>
#include <map>
#include <fstream>

int main (int argc, char **argv){

	std::string names_file="", input_file="";
	int indX=-1, indY=-1;
	std::string namex="", namey="";

	Environment env;
	env.set_name("call_relatedness");
        env.set_version(VERSION);
        env.set_author("Matthew Ackerman");
        env.set_description("A POS relatedness caller. Please direct questions to matthew.s.ackerman@gmail.com");

        env.optional_arg('x',"namey",  namex,      "please provide a number.", "number of individuals in the populations.");
        env.optional_arg('y',"namex",  namey,      "please provide a number.", "number of individuals in the populations.");
        env.optional_arg('i',"input",  input_file,      "please provide a number.", "number of individuals in the populations.");
        env.positional_arg('n',"names",  names_file,      "please provide a number.", "number of individuals in the populations.");

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);

	State Pstates;
 	{
		Flat_file <State> state_file;
		if (input_file=="")
		{
			state_file.open(READ);
		}
		else
		{
			state_file.open(input_file.c_str(), READ);
		}
		Pstates=state_file.read_header();
		state_file.read(Pstates);
		state_file.close();
		std::cerr << Pstates.sample_size() << ", " << Pstates.genome_size() << std::endl;
	}


	gsl_matrix *beta = gsl_matrix_alloc (N, G);
	gsl_matrix *delta = gsl_matrix_alloc (N, G);

	for ?
	{
		
		t_beta[3]={alpha*(), alpha*(), alpha*()};
		t_delta[3]={d*(-2.*psqr+D), d*(2.*pq+D), d*(-2*qsqr+D)};
		for ? 
		{
		beta(k,l)=t_beta[?];
		delta(k,l)=t_delta[?];
		}	
	}
	//close();
}
