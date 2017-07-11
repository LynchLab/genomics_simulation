#include <iostream>
#include <cstdio>
#include <cstring>

#include "interface.h"
#include "stream-tools.h"
#include "fasta.h"

int btoi[256]={0};

const char base[4]={'A', 'C', 'G', 'T'};


void sub (std::string &string, int pos, int num)
{
	string[pos]=base[(btoi[string[pos]]+num)%4];
}

int main (int argc, char* argv[])
{
	btoi['C']=1;
	btoi['c']=1;
	btoi['G']=2;
	btoi['g']=2;
	btoi['T']=3;
	btoi['t']=3;

	Environment env;
	env.set_name("mutate");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman.");
	std::string var_name, ref_name, new_name;
	env.set_description("A program to simulate mutation on genomes. Please direct questions to matthew.s.ackerman@gmail.com");

	env.flag(	'O',"opts", 	&env, 		&flag_options, "an error occurred while displaying the help message", "Prints a list of available options");		//DONE
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message", "Prints this message");				//DONE
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occurred while displaying the version message", "Prints the program version");		//DONE
        env.optional_arg('v', "var",  var_name, "an error occurred while displaying the help message.", "variant filename (default stdin)");
        env.required_arg('r', "ref",  ref_name, "an error occurred while displaying the help message.", "reference filename (default none)");
        env.required_arg('n', "name",  new_name, "an error occurred while displaying the help message.", "new name for fasta (default none)");

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);
	
	std::istream *in=&std::cin;
	std::ostream *out=&std::cout;
	
	std::fstream fin;
	fin.open(ref_name);

	Fasta reference;
	reference.read_fasta(&fin);

	if (new_name.size()>0) reference.set_name(new_name);
	
	int pos, num;
	std::vector <std::string> fields;
	std::string line;

	while (std::getline(*in, line)){
		fields=split(line, ' ');
		pos=atoi(fields[1].c_str() );
		num=atoi(fields[2].c_str() );
		sub(reference.records[0].sequence, pos-1, num);
	}

	reference.write_fasta(out);
	return 0;
}
