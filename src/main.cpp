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
  Tailor uses BWT to perform genomic mapping with ability to detect non-templated
addition of nucleotide to the 3' end of the query sequence (tailing).
  All hits will be aligned to a reference sequence with exact match. Any unmapped
sequences at the 3' end are considered "tail". The exact matching process is 
equivalent to -v 0 -a mode of bowtie.
    Tailor also offer to allow mismatches in the middle of the query string. But
this is not the default behavior.
    Reports will be in SAM format. Tails will be described as "soft-clip" in CIGAR
and the sequences are reported under "TL:Z:" in the optional fields. Mismatches, if
allowed, will be reported in the "MD" tag.

    Tailor is freely avaible on github: jhhung.github.com/Tailor

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
	} else {
		cerr << "Error: unrecognized option " << argv[1] << endl;
		cerr << usage << endl;
		exit (1);
	}
}
