#ifndef BOWHAN_HPP_
#define BOWHAN_HPP_
#include "boost/filesystem.hpp"
#include "boost/thread.hpp"
#include "abwt_thread.hpp"
void buildBWT2 (const std::string& fileName, const std::string& prefixName) {
	/* read input fasta file */
	std::ifstream in {fileName};
	/* string to store the sense + reverse complementary of the genome seq */
	std::string seq, seqRC {};
	/* running accumulator recording the length of each chr */
	INTTYPE tempLen {0}, accumulatedLength {0};
	/* for concatenated seq */
	std::map <INTTYPE, INTTYPE> NPosLen { };

	/* file to store which regions has which chr*/
	std::ofstream chrStartPos {prefixName + "chrStart"};
	/* file to store the length of each chr */
	std::ofstream chrLen {prefixName + "chrLen"};
	/* read in each fasta and make two string */
	while (in.good ()) {
		Fasta<std::vector> fa {in};
		/* store start position of each chr */
		chrStartPos << fa.getName () << '\t' << accumulatedLength << '\n';
		/* get chr length */
		tempLen = fa.getLengthNoN ();
		/* store chr length */
		chrLen << fa.getName () << '\t' << tempLen << '\n';
		/* update accumulated length */
		accumulatedLength += tempLen;
		/* update NPosLen */
		fa.updateNpos (NPosLen);
		seq += fa.getSeqNoN ();
	}
	chrStartPos.close ();
	chrLen.close ();
	/* resize to enough space for the reverse complemetary sequence and a $ sign */
	seq.resize (seq.size () * 2 + 1); // TODO: resize does mallocating the extra space and also initialization, the later is not necessary
	auto iter = seq.begin ();
	std::advance (iter, (seq.size ()-1)/2);
	auto iter2 = iter;
	--iter2;
	do {
		switch (*iter2) {
		case 'A':
			*iter = 'T'; break;
		case 'T':
			*iter = 'A'; break;
		case 'G':
			*iter = 'C'; break;
		case 'C':
			*iter = 'G'; break;
		}
		++iter;
	} while (iter2-- != seq.begin ());
	*iter = '$';
	/* writing NPosLen to file */
	{
		boost::iostreams::filtering_ostream fos;
		fos.push (boost::iostreams::zlib_compressor());
		fos.push (boost::iostreams::file_sink (prefixName + "NposLen.z"));
		boost::archive::binary_oarchive oa (fos);
		oa << NPosLen;
	}
	{
		ABSequence<std::string> x ( seq );
		ABWT<ABSequence<std::string>> y (x, 512, 64, prefixName);
	}
}

bool checkIndexIntact (const std::string& prefixName) {
	if (boost::filesystem::exists (prefixName + "t_bwt.bwt") &&
	boost::filesystem::exists (prefixName + "t_table.bwt") &&
	boost::filesystem::exists (prefixName + "t_seq.bwt") &&
	boost::filesystem::exists (prefixName + "NposLen.z") &&
	boost::filesystem::exists (prefixName + "chrStart") &&
	boost::filesystem::exists (prefixName + "chrLen") ) {
		return true;
	} else {
		return false;
	}
}
// read dual-BWT
ABWT_table loadBWT2 (const std::string& prefixName, std::ostream* out) {
	ABWT_table abwtt;
	ABSequence<std::string> seq;
	abwtt.readBWT(prefixName + "t_bwt.bwt");
	abwtt.readTable(prefixName + "t_table.bwt");
	abwtt.readSEQ(prefixName + "t_seq.bwt", seq);
	abwtt.readNPosLen (prefixName + "NposLen.z");
	abwtt.readChrStartPos (prefixName + "chrStart");
	abwtt.readChrLen (prefixName + "chrLen");
//	abwtt.using_jbwt();
	/* continue to write header */
	for (const auto& chrSizes : abwtt.chr_length) {
		*out << "@SQ\tSN:" << chrSizes.first << "\tLN:" << chrSizes.second << '\n';
	}
	return abwtt;
}


// tailing searching with dual strand
void searchBWT_tail2 (ABWT_table&& abwtt, std::string fileName, std::size_t nthreads, std::ostream* out, int minLen, int allow_mismatch) {
	std::ifstream in {fileName};
/*	
	boost::mt19937 rng;
	rng.seed (static_cast<unsigned int>(std::time(0) + getpid ()));
	boost::uniform_int<> uinInt (1,std::numeric_limits<int>::max());
	boost::variate_generator<boost::mt19937&, boost::uniform_int<> > vg (rng, uinInt);
	std::string randPrefix = std::to_string ( vg ());
*/
	boost::thread* threads [nthreads];
	for (int i = 0 ; i < nthreads; ++i) {
		threads[i] = new boost::thread {ABWT_threads<ABWT_table> {abwtt, &in, out, minLen, allow_mismatch}};
	}
	for (int i = 0 ; i < nthreads; ++i) {
		if (threads[i]->joinable ())
			threads[i]->join ();
		delete threads[i];
	}
}

// tailing version for dual BWT
void tailing2 (const std::string prefixName, const std::string fastqName, std::ostream* out, std::size_t nthread, int minLen, int allow_mismatch) {
	/* writting sam header */
	*out << "@HD" << '\t' << "VN:1.0" << '\t' << "SO:unsorted\n";
	searchBWT_tail2 (loadBWT2 (prefixName, out), fastqName, nthread, out, minLen, allow_mismatch);
}

#endif /* BOWHAN_HPP_ */
