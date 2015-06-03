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

#ifndef ABWT_FORMAT_HPP_
#define ABWT_FORMAT_HPP_

/* predeclaration of Fasta class for Segment class */
template <template <typename, typename ...> class>
class Fasta;

/* Segmen class */
// a class that stores Ns and ACGTs separately, one Segment is basically: NNN...NNNACGT...ACGT
// there could be one segment that does not have N, then _offset == 0
// there could be one segment that has only N, then _len == 0
// one Fasta can contain many Segments
struct Segment {
public:
	INTTYPE _offset {0};
	INTTYPE _len {0};
	std::string _seq {};
	Segment& operator = (const Segment&);
public:
	explicit Segment (std::istream& is) {
		bool finish = false;
		char c;

		while (is.peek() !='>' && is.good () && !finish) {
			c = is.get ();
			switch (c) {
				case 'A': case 'a': case 'C': case 'c':
				case 'G': case 'g': case 'T': case 't': 
					is.putback(c);
					finish = true;
					break;
				case '\r': case '\n': case ' ':
					break;
				case 'N': case 'n':
					++ _offset; break;
				default:
					std::cerr << "unexpected char: " << c << "; turning into N" << std::endl;
					++ _offset; break;
			}
		}

		finish = false;

		while (is.peek () != '>' && is.good() && !finish) {
			c = is.get ();
			switch(c) {
				case 'A': case 'a': case 'C': case 'c':
				case 'G': case 'g': case 'T': case 't':
					++ _len;
					_seq += std::toupper (c);
					break;
				case '\r': case '\n': case ' ':
					break;
				case 'N': case 'n':
					is.putback (c);
					finish = true;
					break;
				default:
					std::cerr << "unexpected char: " << c << "; turning into N" << std::endl;
					is.putback ('N');
					finish = true;
					break;
			}
		}
		//std::cerr << "_len:\t" << _len << "\t_offset:\t" << _offset << std::endl;
	}

	Segment (const Segment&) = delete ;

	Segment (Segment&& other):
		_offset (other._offset),
		_len (other._len),
		_seq {std::move (other._seq)}
	{}

	Segment& operator = (Segment&& other) {
		if (this != &other) {
			_offset = other._offset ;
			_len = other._len;
			_seq.swap (other._seq);
		}
		return *this;
	}
}; /* end of class Segment definition */

/* definition of Fasta class */
template <template <typename, typename ...> class CONTAINER = std::vector>
class Fasta {
private:
	std::string _name {};
	CONTAINER<Segment> _sequences {}; /// each segment is N...NACGT...ACGT until next N
	INTTYPE _length {0};
	INTTYPE _lengthNoN {0};
	Fasta& operator = (const Fasta&);
public:
	class badFasta {};
	Fasta () = default ;
	explicit Fasta (std::istream & is) {
		std::string tmp {};
		char c;
		if (is.peek () == '>') {
			is.ignore (); /// consume the '>'
			// getline (is, _name); /// read the rest as the name
			while (is.get (c)) {
				if ( c == '\t' || c ==' ') {
					is.ignore (100000, '\n');
					break;
				}
				else if ( c == '\n' )
					break;
				else
					_name += c;
			}
			while (is.peek () != '>' && is.good ()) { /// if not ready for reading next Fasta and the stream is good
				_sequences.emplace_back (is); /// continue to read fragments
				auto& _seg = _sequences.back ();
				_length += _seg._len + _seg._offset ;
				_lengthNoN += _seg._len;
			}
		} else
			throw badFasta (); /// if next char is not '>', then throw
	}
	Fasta (const Fasta&) = delete;
	Fasta (Fasta&& other):
		_name (std::move (other._name)),
		_sequences (std::move (other._sequences))
	{}
	Fasta& operator = (Fasta&& other) {
		if (this != &other) {
			_name.swap (other._name);
			_sequences.swap (other._sequences);
		}
		return *this;
	}

	void updateNpos (std::map <INTTYPE, INTTYPE>& NposLen) const {
		auto seg = _sequences.begin ();
		auto stopIter = _sequences.end(); advance (stopIter, -1);
		std::pair <std::map <INTTYPE, INTTYPE>::iterator, bool> lastPos2  {};
	// adding first segment
		/* first chromosome */
		if (NposLen.empty ()) {
			lastPos2 = NposLen.insert ( std::make_pair (0, seg->_offset) ); /// seg->_offset is consumed in current cycle
//			std::cerr << "Just inserted:\t" << (lastPos2.first)->first <<'\t' << (lastPos2.first)->second << '\n';
		}
		/* non first chromosome */
		else {
			auto lastPos = NposLen.rbegin ();  /// get the last entry, update its N
			lastPos->second = seg->_offset; /// the ->first of last entry is gonna used by this new fasta, but we won't get the information of the N until now, hereby we update it
			lastPos2 = NposLen.insert ( std::make_pair (lastPos->first, seg->_offset) ); /// insert a new one, update its _offset. (but _len won't get updated until next segment...)
		}
	// adding middle segment, every time we meet a new segment, we immediately use its _offset. But its _len won't get used until next circle
		while (seg != stopIter) {
//			std::cerr << "(lastPos2.first)->first:\t" << (lastPos2.first)->first << '\n';
//			std::cerr << "seg->_len:\t" << seg->_len << '\n';
//			std::cerr << "seg->_offset:\t"<< seg->_offset<<'\n';
			auto tmp = (lastPos2.first)->first + seg->_len;
			lastPos2 = NposLen.insert (std::make_pair ( tmp , (lastPos2.first)->second + (++seg)->_offset  ));
//			std::cerr << "Just inserted:\t" << (lastPos2.first)->first <<'\t' << (lastPos2.first)->second << '\n';
		}
		if (seg->_len != 0) { /// if last segment has ATGC, we have to update their _len
			NposLen.insert (std::make_pair ( (lastPos2.first)->first + seg->_len,  (lastPos2.first)->second )); // this will be updated in next chr
		}

//		for (const auto x : NposLen) {
//			std::cerr << x.first << '\t' << x.second << std::endl;
//		}
	}

	friend std::ostream& operator << (std::ostream& os, const Fasta& fasta) {
		os << '>' << fasta._name << '\n';
		for (const Segment& seg : fasta._sequences) {
			for (int i = 0; i < seg._offset; ++i) {
				os << 'N';
			}
			os << seg._seq;
		}
		os << '\n';
		return os;
	}
	std::string getName () const {
		return _name;
	}
	std::string getSeq () const {
		std::string _sequence {};
		for (const Segment& seg : _sequences) {
			for (int i = 0; i < seg._offset; ++i) {
				_sequence += 'N';
			}
			_sequence += seg._seq;
		}
		return _sequence;
	}
	std::string getSeqNoN () const {
		std::string _sequence;
		for (const auto& seg : _sequences) {
			_sequence += seg._seq;
		}
		return _sequence;
	}
	std::string getReverseSeq () const {
		std::string tmp {this->getSeqNoN()};
		return std::string {tmp.crbegin(), tmp.crend()};
	}

	std::string getReverseComplementSeq () const {
		std::string _sequence {this->getReverseSeq()};
		for (auto iter = _sequence.begin(); iter!= _sequence.end(); ++ iter) {
			switch (*iter) {
			case 'A':
				*iter = 'T'; break;
			case 'T': 
				*iter = 'A'; break;
			case 'G':
				*iter = 'C'; break;
			case 'C':
				*iter = 'G'; break;
			case 'N':
				break;
			default:
				throw badFasta ();
			}
		}
		return _sequence;
	}
	INTTYPE getLength () const {
		return _length;
	}
	INTTYPE getLengthNoN () const {
		return _lengthNoN;
	}
}; /* end of class Fasta definition */

/* definition for class Fastq */
class Fastq {
private:
	std::string _name {};
	std::string _sequence {};
	std::string _quality {};
	Fastq& operator = (const Fastq&);
public:
	class badFastq {};
	Fastq () = default;
	explicit Fastq (std::istream& is) {
		char c;
		if (is.peek() == '@') {
			is.ignore ();
//			std::getline (is, _name);
			while (is.get (c)) {
				if ( c == '\t' || c ==' ') {
					is.ignore (100000, '\n');
					break;
				}
				else if ( c == '\n' )
					break;
				else
					_name += c;
			}
			std::getline (is, _sequence);
			boost::to_upper (_sequence);
			is.ignore (100000, '\n');
			std::getline (is, _quality);
			is.peek ();
		} else if (is.peek() == '>') {
			/// to support fasta
			is.ignore();
			while (is.get (c)) {
				if ( c == '\t' || c ==' ') {
					is.ignore (100000, '\n');
					break;
				}
				else if ( c == '\n' )
					break;
				else
					_name += c;
			}
			std::string temp;
			while(is.peek() != '>' && is.good()) {
				std::getline(is, temp);
				_sequence += temp;
			}
			_quality.resize(_sequence.size(), 'I');
//			is.peek();
		} else {
			throw badFastq ();
		}
	}

	Fastq (const Fastq& other):
		_name {other._name},
		_sequence {other._sequence},
		_quality {other._quality}
		{}
	Fastq (Fastq&& other):
		_name {std::move (other._name)},
		_sequence {std::move (other._sequence)},
		_quality {std::move (other._quality)}
		{}
	Fastq& operator=(Fastq&& other) {
		if (this != &other) {
			_name.swap (other._name);
			_sequence.swap (other._sequence);
			_quality.swap (other._quality);
			}
			return *this;
		}
		int seq_size() const
		{
			return _sequence.size();
		}
		std::string getName () const {
			return _name;
		}
		std::string getSeq () const {
			return _sequence;
		}
		std::string getQuality () const {
			return _quality;
		}
		std::string getRevQuality () const {
			return std::string {_quality.rbegin(), _quality.rend()};
		}
		friend std::ostream& operator << (std::ostream& os, const Fastq& fastq) {
			os << '@' << fastq._name << '\n' << fastq._sequence << "\n+\n" << fastq._quality << '\n';
		}
}; /* end of class Fastq definition */

/* definition of class Sam */
class Sam {
public:
	enum SAM_FLAG {
		MAPPED = 0,
		PAIRED_END = 1,
		EACH_END_ALIGNED = 2,
		UNMAPPED = 4,
		NEXT_UNMAPPED = 8,
		REVERSE_COMPLEMENTED = 16,
		NEXT_REVERSE_COMPLEMENTED = 32,
		FIRST_SEG = 64,
		SECOND_SEG = 128,
		SECONDARY_ALIGNMENT = 256,
		NOT_PASSING_QUALITY = 512,
		PCR_DUP = 1024,
		FLAG_SIZE = 11
	};
private:
	std::string QNAME = "";
	std::bitset<SAM_FLAG::FLAG_SIZE> FLAG;
	std::string RNAME = "";
	INTTYPE POS = 0;
	int MAPQ = 255;
	std::string CIGAR = "";
	std::string RNEXT = "";
	INTTYPE PNEXT = 0;
	int64_t TLEN = 0;
	std::string SEQ = "";
	std::string QUAL = "";
	std::string _MD = "";
	INTTYPE _NH = 0;
	std::string _tailSeq = "";
//	std::unordered_map<std::string, std::string> OPTIONAL_FIELDS {};

public:
	/**! default constructor **/
	Sam () = default ;
	/** ctor from individuals*/
	Sam (std::string&& _QNAME, SAM_FLAG _FLAG, std::string&& _RNAME, INTTYPE _POS, int _MAPQ, std::string&& _CIGAR, std::string&& _RNEXT, INTTYPE _PNEXT, int64_t _TLEN, const std::string& _SEQ, std::string&& _QUAL, INTTYPE NH, std::string&& tailSeq = "", std::string&& MD=""):
		QNAME {_QNAME},
		FLAG {_FLAG},
		RNAME {_RNAME},
		POS {_POS},
		MAPQ {_MAPQ},
		CIGAR {_CIGAR},
		RNEXT {_RNEXT},
		PNEXT {_PNEXT},
		TLEN {_TLEN},
		SEQ {_SEQ},
		QUAL {_QUAL},
		_NH {NH},
		_tailSeq {tailSeq},
		_MD {MD}
		{ }
		/**! copy ctor **/
		Sam (const Sam&) = default;
		/**! move ctor **/
		Sam (Sam&&) = default;
		//	/**! lvalue assignment operator **/
		//	Sam& operator=(const Sam& other) {
		//		if (this == &other) return *this;
		//		QNAME = other.QNAME;
		//		FLAG  = other.FLAG;
		//		RNAME = other.RNAME;
		//		POS   = other.POS;
		//		MAPQ  = other.MAPQ;
		//		CIGAR = other.CIGAR;
		//		RNEXT = other.RNEXT;
		//		PNEXT = other.PNEXT;
		//		TLEN  = other.TLEN;
		//		SEQ   = other.SEQ;
		//		QUAL  = other.QUAL;
		//		OPTIONAL_FIELDS = other.OPTIONAL_FIELDS;
		//		return *this;
		//	}
		/**! rvalue assignment operator **/
		Sam& operator=(Sam&& other) {
			if (this != &other) {
				QNAME.swap (other.QNAME);
				FLAG = other.FLAG;
				RNAME.swap (other.RNAME);
				POS   = other.POS;
				MAPQ  = other.MAPQ;
				CIGAR.swap (other.CIGAR);
				RNEXT.swap (other.RNEXT);
				PNEXT = other.PNEXT;
				TLEN  = other.TLEN;
				SEQ.swap (other.SEQ);
				QUAL.swap (other.QUAL);
				_NH = other._NH;
				_tailSeq.swap (other._tailSeq);
				_MD.swap(other._MD);
				other.~Sam ();
			}
			return *this;
		}
		/**! destructor **/
		~Sam () = default;

		template <typename T>
		bool checkFlag (T f) const {
			return FLAG.test (f);
		}

		template <typename T, typename ... Args>
		bool checkFlag(T f, Args ... rest) const {
			if (sizeof ...(rest))
				return FLAG.test(f) && checkFlag (rest...);
			return FLAG.test (f);
		}
		// print out the sequence of SEQ
		int get_size () const {
			return SEQ.size ();
		}
		// output
		friend std::ostream& operator<< (std::ostream& os, const Sam& sam) {
			os << sam.QNAME << '\t'
					<< sam.FLAG.to_ulong() << '\t'
					<< sam.RNAME << '\t'
					<< sam.POS << '\t'
					<< sam.MAPQ << '\t'
					<< sam.CIGAR << '\t'
					<< sam.RNEXT << '\t'
					<< sam.PNEXT << '\t'
					<< sam.TLEN << '\t'
					<< sam.SEQ << '\t'
					<< sam.QUAL << '\t'
					//<< "MD:Z:" << sam._MD << '\t'
					<< "NH:i:" << sam._NH;
			if (!sam._MD.empty ())
				os << "\tMD:Z:" << sam._MD;
			if (!sam._tailSeq.empty ())
				os << "\tTL:Z:" << sam._tailSeq;
			os << '\n';
			return os;
		}
};

#endif /* ABWT_FORMAT_HPP_ */
