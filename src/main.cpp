#include "../include/includes.hpp"
#include <boost/filesystem.hpp>

using namespace std;
int main (int argc, char** argv) 
{
	string usage = R"(

*********************************************************************************
+------+
|Tailor|
+------+
  A software using BWT to perform genomic mapping with ability to detect 
untemplated addition of nucleotide to the 3' end of small RNA (tailing). 
  All hits will be aligned to a reference sequence with exact match. Any unmapped
sequences at the 3' end are considered "tail". The exact matching process is 
equivalent to -v 0 -a mode of bowtie.
  Reports will be in SAM format. Tails will be described as "soft-clip" in CIGAR
and the sequences are reported under "TL:Z:" in the optional fields.

Usage:

 1. building index of the genome
>	tailor build --help
  
 2. mapping fastq to the index
>	tailor map --help
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
