#include "mpi_fast.h"

#define BUFFER_SIZE 500

#ifndef NOGSL

#define MASTER 0 

size_t getx(const size_t &w, const size_t &size)
{
	size_t t=0;
	for (size_t x=0; x<size; x++){
		for (size_t y=x+1; y<size; y++){
			if (w==t) return x;
			t++;
		}
	}
	return 0;	
}

size_t gety(const size_t &w, const size_t &size)
{
	size_t t=0;
	for (size_t x=0; x<size; x++){
		for (size_t y=x+1; y<size; y++){
			if (w==t) return y;
			t++;
		}
	}
	return 0;
}

int estimateRel(int argc, char *argv[])
{
	int taskid,	/* task ID - also used as seed number */
	numtasks,	/* number of tasks */
	rc,             /* return code */
	i;		/* new comment */

	MPI_Status status;
	/* Obtain number of tasks and task ID */
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

#if DEBUG
	printf ("MPI task %d has started...\n", taskid);
#endif

	/* All the variables that can be set from the command line */

	std::string in_name="", rel_name="";

	Environment env;
	env.set_name("mapgd relatedness");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Uses called genotypes to estimate pairwise relatedness.");
	env.optional_arg('i', "input", 	in_name, "an error occurred while displaying the help message.", "input filename (default stdout)");
	env.optional_arg('o', "output", rel_name, "an error occurred while displaying the help message.", "output filename (default stdin)");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the help message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Indexed_file <Population> in; 		// Open the file with genotypic probabilities.
	Indexed_file <Population> in_mem; 	// the in memory gcf_file.

	std::stringstream file_buffer;

	Flat_file <Relatedness> rel_out; 	// Open a file for output.

	Relatedness relatedness;		//The class to read to
	Population genotype;			//The class to write from

	size_t file_buffer_size=0;
	char *char_buff;
	size_t sample_size;
	sample_names

	if(taskid==MASTER){
		if (in_name.size()!=0)
			in.open(in_name.c_str(), std::ios::in);
		else
			in.open(std::ios::in);

		if (rel_name.size()!=0)
			rel_out.open(rel_name.c_str(), std::ios::out);
		else
			rel_out.open(std::ios::out);
	
		genotype=in.read_header();			//This gives us the sample names.

		in_mem.open(&file_buffer, std::ios::out );
		in_mem.set_index(in.get_index() );
		in_mem.write_header(genotype);

		while(in.table_is_open() ){
			in.read(genotype);
			in_mem.write(genotype);
		}
		in_mem.close();
		rel_out.write_header(relatedness);
		sample_size=genotype.get_sample_names().size();
		file_buffer_size=file_buffer.str().size();
		char_buff=new char [file_buffer_size];
		memcpy(char_buff, file_buffer.str().c_str(), file_buffer_size);
	}
	MPI_Bcast(&file_buffer_size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&sample_size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	if (taskid!=MASTER) char_buff=new char [file_buffer_size];
	MPI_Bcast(char_buff, file_buffer_size, MPI_CHAR, MASTER, MPI_COMM_WORLD);
	file_buffer.str(char_buff);

	//TODO: This is uninitialized in most tasks.
	MPI_Datatype record;
	MPI_Type_contiguous(sizeof(relatedness), MPI_BYTE, &record);
	MPI_Type_commit(&record);

	//TODO replace std::map <Genotype_pair_tuple, size_t> hashed_genotypes;
	//TODO replace std::map <Genotype_pair_tuple, size_t> down_genotypes;
	
	Relatedness buffer_rel[BUFFER_SIZE];        
	std::fill_n(buffer_rel, BUFFER_SIZE, relatedness);

	//MPI_Send ( Data, count, data_type, root,	   0, communicator)
	//From process to root.
	//MPI_Recv ( Data, count, data_type, i, 0, comunicator, MPI_STATUS_IGNORE)

	size_t chunk=BUFFER_SIZE/numtasks;
	size_t used_size=chunk*numtasks;
	size_t done=0;
	size_t todo=sample_size*(sample_size-1)/2;
	size_t z=0;

	for (size_t x=0; x<sample_size; ++x){
		for (size_t y=x+1; y<sample_size; ++y){
			if (z>=used_size) 
			{
				if (taskid == MASTER) 
				{
					for (size_t i=1; i<numtasks; i++)
					{
						MPI_Recv (&buffer_rel[i*chunk], chunk, record, i, 0, MPI_COMM_WORLD, &status);
					}
					for (size_t w=0; w<used_size; w++)
					{
						size_t X=getx(w+done, sample_size);
						size_t Y=gety(w+done, sample_size);
						relatedness=buffer_rel[w];
						relatedness.set_X_name(name[X]);
						relatedness.set_Y_name(name[Y]);
						rel_out.write(relatedness);
					}
				} else {
					MPI_Send (&buffer_rel[taskid*chunk], chunk, record, MASTER, 0, MPI_COMM_WORLD);
				}
				done+=used_size;
				z-=used_size;
			} 
			if ( chunk*taskid <= z && z < chunk*(taskid+1) )
			{
				relatedness.zero();

				get_relatedness(relatedness, ?);

				buffer_rel[z]=relatedness;
			}
			z++;
		}
	}
	if (taskid == MASTER) 
	{
		for (size_t i=1; i<numtasks; i++) 
		{
			MPI_Recv (&buffer_rel[i*chunk], chunk, record, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		for (size_t w=0; w<todo-done; w++){
			size_t X=getx(w+done, sample_size);
			size_t Y=gety(w+done, sample_size);
			relatedness=buffer_rel[w];
			relatedness.set_X_name(name[X]);
			relatedness.set_Y_name(name[Y]);
			rel_out.write(relatedness);
		}
		rel_out.close();
	} else {
		MPI_Send (&buffer_rel[taskid*chunk], chunk, record, MASTER, 0, MPI_COMM_WORLD);
	}
#if DEBUG
	printf ("MPI task %d has ended...\n", taskid);
#endif
	MPI_Finalize();
	return 0;					//Since everything worked, return 0!.
}

#else 

int 
estimateRel(int argc, char *argv[])
{
	std::cerr << "This command depends on gsl, which could not be found. Please e-mail matthew.s.ackerman@gmail.com for help.\n";
	return 0;
}

#endif 

