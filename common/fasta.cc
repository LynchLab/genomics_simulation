#include "fasta.h"

Fasta_record::Fasta_record (const std::string &line)
{
	desc=line.substr(1);
}

void Fasta::new_record(const std::string &line)
{
	records.push_back(Fasta_record (line) );
}

void Fasta::extend_record(const std::string &line){
	records.back().sequence.append(line);
}
void Fasta::read_fasta (std::istream *in)
{
	std::string line;
	while (std::getline (*in, line) )
	{
		if (line[0]=='>' )
		{
			new_record(line);
		}
		else {
			extend_record(line);
		}
	}	
}

void Fasta::set_name(const std::string &line)
{
	records[0].desc=line;
}
void Fasta::write_fasta(std::ostream *out)
{
	std::vector<Fasta_record>::iterator end=records.end();
	std::vector<Fasta_record>::iterator it=records.begin();
	size_t start, string_size;
	while (it!=end)
	{
		*out << '>' << it->desc << std::endl;
		start=0;
		string_size=it->sequence.size();
		for (;start<string_size;start+=LINE_LENGTH)
		{
			*out << it->sequence.substr(start, LINE_LENGTH ) << std::endl;
		}
		++it; 
	}
}
