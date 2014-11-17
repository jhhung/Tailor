/*
# Tailor, a BWT-based aligner for non-templated RNA tailing
# Copyright (C) 2014 Min-Te Chou, Bo W Han, Jui-Hung Hung
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

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

# To generate index files from a fasta file.

>  tailor build

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
			("force,f", boost::program_options::bool_switch(&overwrite)->default_value(false), "Overwrite the existing index files if they already exist.")
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
