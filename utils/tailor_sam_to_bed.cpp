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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "htslib/sam.h"
#include "boost/program_options.hpp"

using namespace std;
int main(int argc, char** argv)
{
    const string usage = R"(
This program convert sam/bam output of tailor to bed2 format.
Please do not use it for other purpose.
    )";
    string input_sam_file;
    string output_bed_file;
    int max_tail_len;
    boost::program_options::variables_map vm;
    boost::program_options::options_description opts {usage};
    try {
        opts.add_options ()
        ("help,h",		"display this help message and exit")
        ("input,i",		boost::program_options::value<std::string>(&input_sam_file )->default_value(string{"-"} ), "input SAM file")
        ("output,o",	boost::program_options::value<std::string>(&output_bed_file)->default_value(string{"-"}), "output BED file")
        ("maxtail,m",	boost::program_options::value<int>(&max_tail_len)->default_value(6), "maximal length allowed for the tail")
        ;
        boost::program_options::store (boost::program_options::parse_command_line(argc, argv, opts), vm);
        boost::program_options::notify(vm);
        if (vm.count("help"))	{ std::cerr << opts << std::endl; exit (1); }
    }
    catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << opts << std::endl;
        exit (1);
    } catch (...) {
        std::cerr << "Unknown error!" << std::endl;
        std::cerr << opts << std::endl;
        exit (1);
    }

		samFile* in1 = sam_open(input_sam_file.c_str(), "r");
    if (in1 == NULL) {
				cerr << "error openning " << input_sam_file << endl;
        return EXIT_FAILURE;
    }
    bam_hdr_t* header = sam_hdr_read(in1);
    bam1_t* aln = bam_init1();
    int exit_code = 0;
    uint8_t* NH_s;
    uint8_t* TL_s;
    uint8_t* MD_s;
    while (sam_read1(in1, header, aln) >= 0) {
        NH_s = bam_aux_get(aln, "NH");
        TL_s = bam_aux_get(aln, "TL");
        MD_s = bam_aux_get(aln, "MD");
//        uint32_t *cigar = bam_get_cigar(aln);
//        for (int i = 0; i < aln->core.n_cigar; ++i) {
//            cout << bam_cigar_oplen(cigar[i]) << bam_cigar_opchr(cigar[i]);
//        }
        fprintf(stdout, "%s\t%d\t%d\t%s\t%d\t%c\t%s\t%s\n",
            header->target_name[aln->core.tid], // query name
            aln->core.pos , // position
            bam_endpos(aln), // aln->core.pos + aln->core.l_qseq + aln->core.qual - 255
            bam_get_qname(aln),
            NH_s == NULL ? 0 : bam_aux2i(NH_s),
            aln->core.flag & 16 ? '-' : '+',
            TL_s == NULL ? "*" : bam_aux2Z(TL_s),
            MD_s == NULL ? "*" : bam_aux2Z(MD_s)
        );
    }
    bam_destroy1(aln);
    bam_hdr_destroy(header);
    sam_close(in1);
    return exit_code;
}
