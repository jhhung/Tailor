#include "../include/includes.hpp"
#include <boost/filesystem.hpp>

using namespace std;
int main (int argc, char** argv) 
{
	string usage = R"(

*********************************************************************************
tailor

  A software using BWT to perform genomic mapping with ability to detect tailing. 

  All hits will be reported as SAM format, equivalent to -v 0 -a mode of bowtie,
except for the ability to detect unmatched sequence in the 3' end.

Usage:

 1. building index

>	tailor build
  
 2. mapping using fastq as input
 
>	tailor map

*********************************************************************************

)";
	if (argc < 2) {
		cerr << usage << endl;
		exit (1);
	}
	if (strcmp (argv[1], "build") == 0) {
		tailor_build::main (argc-1, argv+1);
	}
	else
	if (strcmp (argv[1], "map") == 0) {
		tailor_map::main (argc-1, argv+1);
	}
}
