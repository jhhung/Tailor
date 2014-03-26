#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <boost/program_options.hpp>

using namespace std;

string reverse_complement (const string&);

int main (int argc, char** argv) {
	std::string usage = R"(
****************************
this script filer the SAM output of Tailor and convert 
it to BED2 format. It is specifically designed to deal
with Tailor's SAM output thus misses a lot of format/safty check. 
Please don't use it for other purposes.

****************************

)";
	/** option variables **/
	std::string	input_sam_file	{};
	std::string	output_bed_file	{};
	int		max_tail_len		{};
	boost::program_options::variables_map vm;
	boost::program_options::options_description opts {usage};
	try {
		opts.add_options ()
		("help,h",		"display this help message and exit")
		("input,i",		boost::program_options::value<std::string>(&input_sam_file )->default_value(string{"stdin"} ), "input SAM file")
		("output,o",	boost::program_options::value<std::string>(&output_bed_file)->default_value(string{"stdout"}), "output BED file")
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
	std::istream* in {&std::cin};
	if (input_sam_file != "stdin" && input_sam_file != "-" ) {
		in = new std::ifstream {input_sam_file};
	}
	std::ostream* out {&std::cout};
	if (output_bed_file!="stdout" && output_bed_file!="-") {
		out = new std::ofstream {output_bed_file};
	}
	/** reading header **/
	string line;
	char c {in->peek ()};
	while (c == '@' && in->good ()) {
		getline (*in, line);
		// *out << line << '\n';
		c = in->peek ();
	}
	/** reading SAM **/
	
	/** for sam **/
	string qname;	// chromosome
	int flag;		// flag
	string rname;	// reads name
	int pos;		// mapping position, 1-based
	int mapq;		// mapping quality, 255 - tail length
	string cigar;	// cigar string
	string rnext;	// "*" for tailor
	int pnext;		// 0 for tailor
	int tlen;		// 0 for tailor
	string seq;		// sequence
	string qual;	// quality
	string NH_str;			// value of HN tag
	string TL_str;		// value of TL tag
	/** for bed **/
	char strand;
	
	while (in->peek () && in->good ()) {
		*in >> qname;
		*in >> flag;
		*in >> rname;
		*in >> pos;
		*in >> mapq;
		*in >> cigar;
		*in >> rnext;
		*in >> pnext;
		*in >> tlen;
		*in >> seq;
		*in >> qual;
		*in >> NH_str;
		if (mapq != 255) {
			*in >> TL_str;
		}

		if ( (255 - mapq) <= max_tail_len) {
			*out <<	rname << '\t'
				<<	pos - 1 << '\t'
				<<	pos + seq.size () + mapq - 255 - 1 << '\t'
				<<	qname << '\t'
				<<	NH_str.substr (NH_str.find_last_of (":")+1) << '\t';
				
			if (flag & 0x10) { 
				*out << "-\t"
					<< reverse_complement (seq) << '\t';
			} else {
				*out << "+\t"
					<<  seq << '\t';
			}
			if (mapq!=255) {
				*out <<	TL_str.substr (TL_str.find_last_of (":")+1) << '\t';
			} else {
				*out << '*' << '\t';
			}
			*out <<	255 - mapq << '\n';	
		}
	}
	if (in != &std::cin) {
		static_cast<std::ifstream*> (in)->close ();
		delete in;
	}
	if (out != &std::cout) {
		static_cast<std::ofstream*>(out)->close ();
		delete out;
	}
	return 0;
}

string reverse_complement (const string& sense) {
	string reverse {sense.rbegin (), sense.rend ()};
	for (char& c : reverse) {
		switch (c) {
			case 'A':
			case 'a':
			c = 'T'; break;
			case 'T':
			case 't':
			c = 'A'; break;
			case 'C':
			case 'c':
			c = 'G'; break;
			case 'G':
			case 'g':
			c = 'C'; break;
		}
	}
	return reverse;
}
