#ifndef TAILOR_BUILD_HPP_
#define TAILOR_BUILD_HPP_

#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "tailer.hpp"

namespace tailor_build {

int main (int argc, char** argv) {
	std::string usage = R"(

*********************************************************************************
tailor

  A software using BWT to perform genomic mapping with ability to detect tailing. 

  All hits will be reported as SAM format, equivalent to -v 0 -a mode of bowtie,
except for the ability to detect unmatched sequence in the 3' end.

  tailor build
    build index from a fasta file. 

*********************************************************************************

)";

	/** option variables **/
	std::string inputFasta {};
	std::string indexPrefix {};
	bool overwrite {false};
	boost::program_options::options_description opts {usage};
	try {
		opts.add_options ()
			("help,h", "display this help message and exit")
			("input,i", boost::program_options::value<std::string>(&inputFasta)->required(), "The input fasta file.")
			("prefix,p", boost::program_options::value<std::string>(&indexPrefix)->required(), "Prefix of index file to generate.")
			("force,f", boost::program_options::bool_switch(&overwrite)->default_value(false), "Overwrite the existing infex files if they already exist.")
		;
		boost::program_options::variables_map vm;
		boost::program_options::store (boost::program_options::parse_command_line(argc, argv, opts), vm);
		boost::program_options::notify(vm);
	} catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << opts << std::endl;
        exit (1);
    } catch (...) {
        std::cerr << "Unknown error!" << std::endl;
        std::cerr << opts << std::endl;
        exit (1);
    }
	/* test whether input file exist */
	if (!boost::filesystem::exists (inputFasta)) {
		std::cerr << "Error: Input fasta file " << inputFasta << " does not exist! Please double check. Existing..." << std::endl;
		exit (1);
	}

	if (indexPrefix.back () != '.') {
		indexPrefix += '.';
	}
	/* check whether index already exist */
	if (checkIndexIntact (indexPrefix) && !overwrite) {
		std::cerr << "Error: index files already exist. If you want to overwrite them, please run it again with option -f.\nExisting..." << std::endl;
		exit (2);
	}
	/* executing buildBWT */
	buildBWT2 (inputFasta, indexPrefix);
}
}

#endif
