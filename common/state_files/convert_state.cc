#include "interface.h"
#include <string.h>

int main (int argc, char **argv)
{
        Environment env;
	std::string filename="";
        env.set_name("convert_state");
        env.set_version(VERSION);
        env.set_author("Matthew Ackerman.");
        env.set_description("A program to convert state files. Please direct questions to matthew.s.ackerman@gmail.com");

	env.positional_arg(&filename, "an error has occured while setting filename.", "sets the input filename, default stdin");
        env.flag(       'h',"help",     &env,           &flag_help,     "an error occurred while displaying the help message", "Prints this message");                          //DONE
        env.flag(       'v',"version",  &env,           &flag_version,  "an error occurred while displaying the version message", "Prints the program version");                //DONE

        if ( parsargs(argc, argv, env) != 0 ) print_usage(env);	
        return 0;
	
}
