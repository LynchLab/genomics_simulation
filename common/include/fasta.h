#ifndef __FASTA_H__
#define __FASTA_H__

#include <vector>
#include <string>
#include <iostream>

#define LINE_LENGTH 90

struct Fasta_record
{
	Fasta_record(const std::string &);
	std::string desc;
	std::string sequence;
};

class Fasta 
{
private:
public:
	std::vector <Fasta_record> records;
	void new_record(const std::string &);
	void extend_record(const std::string &);
	void read_fasta (std::istream *);
	void write_fasta (std::ostream *);
	void set_name (const std::string &);
};

#endif
