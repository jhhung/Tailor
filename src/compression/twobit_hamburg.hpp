///@file twobit.hpp
///@brief define two-bit compression scheme for obtained sequence. 
#ifndef TWOBIT_HPP_
#define TWOBIT_HPP_

#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <boost/utility/binary.hpp>


class TwoBitSequence
{
    public:
	TwoBitSequence () : sequence (0) { }  /// default constructor
	TwoBitSequence (int length) : sequence (length) { }  /// constructor from a fixed length
	TwoBitSequence (const std::string&); /// constructor from a DNA string
	TwoBitSequence (const TwoBitSequence&); /// copy constructor
	TwoBitSequence (const boost::dynamic_bitset<>&);	
	TwoBitSequence (boost::dynamic_bitset<>&&);
	TwoBitSequence (TwoBitSequence&&); ///move constructor
	TwoBitSequence& operator = (const TwoBitSequence&); /// assignment operator
	TwoBitSequence& operator = (TwoBitSequence&&); /// assignment operator
	TwoBitSequence& reverse (); 	/// reverse sequence
	TwoBitSequence& complement ();	/// complement sequence
	TwoBitSequence& reverse_complement (); /// reverse complement sequence
	TwoBitSequence reverse_copy () ; 	/// reverse sequence
	TwoBitSequence complement_copy () ;	/// complement sequence
	TwoBitSequence reverse_complement_copy () ; /// reverse complement sequence
	std::string to_string () const;
	std::string to_DNA () const;
	char operator [] (uint64_t) const;
	uint64_t size () const;
	std::string substr (uint64_t, uint64_t) const;
    private:
	boost::dynamic_bitset<> sequence;
};

TwoBitSequence::TwoBitSequence (const std::string& str_seq)
{
    uint64_t l = str_seq.size ();
    uint64_t pos {};
    sequence.resize (l<<1);
    for (uint64_t i = 0; i < l ; ++ i)
    {
	pos = (i << 1);
	switch (str_seq[i])
	{
	    case 'A':
	    case 'a':
		sequence[pos]   = 0;
		sequence[pos+1] = 0;
		break;
	    case 'C':
	    case 'c':
		sequence[pos]   = 1;
		sequence[pos+1] = 0;
		break;
	    case 'G':
	    case 'g':
		sequence[pos]   = 0;
		sequence[pos+1] = 1;
		break;
	    case 'T':
	    case 't':
	    case 'U':
	    case 'u':
		sequence[pos]   = 1;
		sequence[pos+1] = 1;
		break;
	    case 'N':
	    case 'n':
		break;
	    default:
		std::cerr << "[Error] : illegal sequence " << str_seq[i] << " in sequence\n";
		exit (1);
		break;
	}
    }
}

std::string TwoBitSequence::to_string () const
{
    std::string buffer;
    buffer.reserve (sequence.size());
    boost::to_string (sequence, buffer);
    return buffer;
}

std::string TwoBitSequence::to_DNA () const
{
    return this->substr (0, sequence.size()>>1);
}

char TwoBitSequence::operator [] (uint64_t p) const
{

    uint64_t pos {p<<1};
    if (pos >= sequence.size())
	throw std::out_of_range ("[Error] : TwoBitSequence : index longer than size");
    bool first = sequence [pos];
    bool second = sequence [pos+1];
    if (first && second)
	return 'T';
    else if (first && !second)
	return 'C';
    else if (!first && second)
	return 'G';
    else if (!first && !second)
	return 'A';
}

uint64_t TwoBitSequence::size () const
{
    return sequence.size ();
}

std::string TwoBitSequence::substr (uint64_t start, uint64_t len) const
{
    std::string sstr;
    sstr.resize (len);
    for (uint64_t i = 0; i < len; ++ i)
    {
	sstr[i] = this->operator [](start+i);
    }
    return sstr;
}


TwoBitSequence::TwoBitSequence (const TwoBitSequence& other)
{
    sequence = other.sequence;
}

TwoBitSequence::TwoBitSequence (TwoBitSequence&& other)
{
    sequence.clear ();
    sequence.swap (other.sequence);
}

TwoBitSequence& TwoBitSequence::operator = (const TwoBitSequence& other)
{
    if (this == & other)
	return *this;
    sequence.clear ();
    sequence = other.sequence;
    return *this;
}

TwoBitSequence::TwoBitSequence (const boost::dynamic_bitset<>& db)
{
    sequence = db;
}

TwoBitSequence::TwoBitSequence (boost::dynamic_bitset<>&& db)
{
    sequence.clear ();
    sequence.swap (db);
}



TwoBitSequence& TwoBitSequence::operator = (TwoBitSequence&& other)
{
    assert (this != & other);
    sequence.clear ();
    sequence.swap (other.sequence);
    other.~TwoBitSequence ();
    return *this;
}

TwoBitSequence& TwoBitSequence::reverse ()
{
    uint64_t len1  { sequence.size () };
    decltype(sequence) temp {len1};	/// to store the result
    decltype(sequence) mask {len1};	/// to mask certain bits
    mask[0] = mask[1] = 1;			/// initialized
    int64_t len = reinterpret_cast<int64_t&> (len1);	/// will need the minus value
    /// human genome is 3095693983nt (10111000100001001000101010011111) 32bits,
    /// thus 64 bits is more than enough and the sign bit will not be used (but uint32_t won't be enough)
    len-=2;
    while (mask.any())	/// mask 11000...000 will be the last mask
    {
	if (len >= 0)
	    temp |= ( (mask&sequence) << len);
	else
	    temp |= ( (mask&sequence) >> -len);
	//		std::cout << "seq  is :" << sequence << '\n';
	//		std::cout << "mask is :" << mask << '\n';
	//		std::cout << "temp is :" << temp << '\n';
	//		std::cout << "len  is :" << len << '\n';
	mask<<=2;	/// left shift to retrieve the next nt
	len-=4;		/// two bits shifted from each side
    }
    sequence.swap (temp);
    return *this;
}

TwoBitSequence& TwoBitSequence::complement ()
{
    sequence.flip ();	/// in-place flip
    return *this;
}

TwoBitSequence& TwoBitSequence::reverse_complement ()
{
    this->reverse ();
    this->complement ();
    return *this;
}

TwoBitSequence TwoBitSequence::reverse_copy ()
{
    uint64_t len1  { sequence.size () };
    decltype(sequence) temp {len1};	/// to store the result
    decltype(sequence) mask {len1};	/// to mask certain bits
    mask[0] = mask[1] = 1;			/// initialized
    int64_t len = reinterpret_cast<int64_t&> (len1);	/// will need the minus part
    /// human genome is 3095693983nt (10111000100001001000101010011111) 32bits,
    /// thus 64 bits is more than enough and the sign bit will not be used (but uint32_t won't be enough)
    len-=2;
    while (mask.any())	/// mask 11000...000 will be the last mask
    {
	if (len >= 0)
	    temp |= ( (mask&sequence) << len);
	else
	    temp |= ( (mask&sequence) >> -len);
	//		std::cout << "seq  is :" << sequence << '\n';
	//		std::cout << "mask is :" << mask << '\n';
	//		std::cout << "temp is :" << temp << '\n';
	//		std::cout << "len  is :" << len << '\n';
	mask<<=2;	/// left shift to retrieve the next nt
	len-=4;		/// two bits shifted from each side
    }
    return TwoBitSequence (std::move(temp));
}

TwoBitSequence TwoBitSequence::complement_copy ()
{
    TwoBitSequence temp {sequence.flip()};
    sequence.flip ();	/// restore the state of sequence, things would be much easier if m_bits is not private...
    return temp;
}

TwoBitSequence TwoBitSequence::reverse_complement_copy ()
{
    TwoBitSequence temp {*this};
    temp.reverse ();
    temp.complement ();
    return temp;
}

#endif
