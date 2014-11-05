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

#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

using namespace std;

int main (int argc, char** argv) {
    std::string usage = R"(
    
    This script does two things
    1. filtering the input fastq file by the phred quality
    2. pool the species with the same sequence together. As soon as they pass the phred threshold, the quality string is dumped.
    The number of reads is writting in the header of the fastq.
    
    )";
    /** option variables **/
    std::string	input_fastq_file {};
    std::string	output_fastq_file {};
    int minimal_phred {};
    int	offset {};
    boost::program_options::options_description opts {usage};
    boost::program_options::variables_map vm;
    try {
        opts.add_options ()
        ("help,h",		"display this help message and exit")
        ("input,i",		boost::program_options::value<std::string>(&input_fastq_file )->default_value(string{"stdin"} ), "input fastq file")
        ("output,o",	boost::program_options::value<std::string>(&output_fastq_file)->default_value(string{"stdout"}), "output fastq file, only with reads passed the threshold. The reads with same sequences are pooled and the number of time been sequenced is in the header")
        ("minphred,q",	boost::program_options::value<int>(&minimal_phred)->default_value(25), "minimal phred score allowed, choose between 0 - 40")
        ("offset,s",	boost::program_options::value<int>(&offset)->default_value(33), "offset between ascii string and the phred. use 33 for sanger and 64 for illumina")
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
    if (minimal_phred > 40) {
        cerr << "error: phred score is too hight, please use 0 - 40" << endl;
        exit (1);
    }
    
    std::ios::sync_with_stdio(false);
    istream* p_ist_in {&std::cin};
    if (input_fastq_file != "stdin" && input_fastq_file != "-" ) {
        p_ist_in = new std::ifstream {input_fastq_file};
    }
    boost::iostreams::filtering_istream in;
    char magic_number[4] = "\0\0\0";
    p_ist_in->get(magic_number, 3);
    if (
           magic_number[0] == '\037'
        && magic_number[1] == (char)'\213'
        && magic_number[2] == '\0'
        ) {
        in.push(boost::iostreams::gzip_decompressor());
    } else
    if (
           magic_number[0] == 'B'
           && magic_number[1] == 'Z'
       )
    {
        in.push(boost::iostreams::bzip2_decompressor());
    } else
    if (magic_number[0] == '@')
    {
//        plain text
    } else {
        cerr << "unknown format" << endl;
        return 1;
    }
    p_ist_in->seekg(0, ios::beg);
    in.push(*p_ist_in);
    std::ostream* out {&std::cout};
    if (output_fastq_file!="stdout" && output_fastq_file!="-") {
        out = new std::ofstream {output_fastq_file};
    }
    char c = in.peek ();
    string sequence, quality;
    unordered_map<string, uint64_t> counter;
    while (in.good ()) {
        /** testing header **/
        if (c!='@') {
            cerr << "illega char " << c << "\nthe input might not be fastq. please check again" << endl;
            return 1;
        }
        /** reading header **/
        in.ignore (numeric_limits<streamsize>::max(), '\n');
        /** reading sequence **/
        getline (in, sequence);
        /** reading third line **/
        in.ignore (numeric_limits<streamsize>::max(), '\n');
        /** reading and filtering quality **/
        bool passed = true;
        for ( in.get(c); c!='\n'; in.get(c) ) {
            if (c < minimal_phred + offset) {
                passed = false;
                break;
            }
        }
        /** incrementing the counter **/
        if (passed) {
            counter[sequence] += 1;
        }
        /** peeking next char, invoke eof **/
        c = in.peek ();
    }
    /** reporting **/
    for (const auto& s : counter) {
        *out << '@' << s.second << '\n'
        << s.first << '\n'
        << "+\n";
        for (int i = 0 ; i < s.first.size (); ++ i)  {
            *out << 'I';
        }
        *out << '\n';
    }
    if (out != &std::cout) static_cast<ofstream*>(out)->close();
    return 0;
}
