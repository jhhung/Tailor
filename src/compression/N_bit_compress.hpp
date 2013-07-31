/// @file n_bit_compress.hpp
/// @brief define a data structure for achieving N bit encode/compression scheme and accordingly saving the encoded result.  Overall, two kinds of classes, i.e. CompressRuleTrait and NBitCompress are provided. \n The CompressRuleTrait is served as a template parameter for the encode/compression core class of NBitCompress, wherein the CompressRuleTrait ACTUALLY conveys the information of how many bits and the actual encode/compression rule that the encode/compression sceme takes. \n The NBitCompress takes the responsibility to carry out the encode/compression operation. 
/// @author JHH Corp.
#ifndef N_BIT_COMPRESS_HPP_
#define N_BIT_COMPRESS_HPP_
#include <iostream>
#include <unordered_map>
#include <string>
#include <cstdint>
#include <deque>
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/deque.hpp"
#include "boost/serialization/string.hpp"
#include "boost/serialization/utility.hpp"
#include "boost/dynamic_bitset.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
//#include "../constant_def.hpp"
#include "compressrule_trait.hpp"

	uint64_t Total=0;

namespace boost 
{	//since we want to access boost::dynamic_bitset<>'s private member of m_bits & m_num_bits, by providing an easy way to serialize 
	//boost::dynamic_bitset<> object, we make use of a boost::dynamic_bitset<>'s friend function, template void to_block_range, to
	//gain access to the aforementioned private members.  We then take a full specialization of that function, to achieve the access
	//of that two priavte members
	template <>
	inline void
	to_block_range<>( const boost::dynamic_bitset<>& b,
					std::vector<boost::dynamic_bitset<>::block_type, boost::dynamic_bitset<>::allocator_type>** a) 
	{
		*a = (std::vector<boost::dynamic_bitset<>::block_type, boost::dynamic_bitset<>::allocator_type>*)(&b.m_bits);
	}

	template <>
	inline void
	to_block_range<>( const boost::dynamic_bitset<>& b, size_t** c ) 
	{
		*c = (size_t*)(&b.m_num_bits);
	}
}

template <int v>
struct Int2Type
{
	enum { value = v};
};

/// @class NBitCompress
/// @brief encode/compression achieves the main encode/compression operation based on the rule designated by the CompressRuleTrait, taken as a template parameter T. 
/// @tparam a type parameter T, i.e. the CompressRuleTrait
template <typename T>
class NBitCompress
{
private:
/// @memberof NBitCompess
/// @brief the data structure holding the encoded binary result
	boost::dynamic_bitset<> sequence;	
/// @memberof NBitCompress
/// @brief the data structure holding the except characters and its position, so that a later decode/decompression operation of the encoded/compressed data can be achieved.
	std::deque < std::pair < uint32_t, std::pair <char, uint32_t> > > ExceptDeq;	
/// @memberof NBitCompress
/// @brief the data structure holding the mask characters and its position, so that a later decode/decompression operation of the encoded/compressed data can be achieved.
	std::deque < std::pair < uint32_t, NBitCompress <CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > > > MaskDeq;
public:
	mutable uint32_t Normal_count, Except_count, Mask_count;
//private functions
/// @memberof NBitCompress
/// @brief taking care of some detail jobe for the operator[] overloading
	char index_impl (uint64_t p) const  //Implement detail for operator[]
	{
		uint64_t pos {p*T::bit_count};
		if (pos >= sequence.size())
			throw std::out_of_range ("[Error] : TWBImpl : index longer than size");
        for (auto i=TraitObject.table.begin(); i!=TraitObject.table.end(); ++i)
        {
			int th = 0;
			while (th < T::bit_count)
			{
				if (i->second.second[th]==sequence[pos+th])
					++th;
				else 
					break;
			}
            if ( (i->second.first == true) && (th==T::bit_count) )
				return i->first;
		}
	}
/// @memberof NBitCompress
/// @brief taking care of some detail jobe for the constructor taking a std::string while the CompressType carried by the inputted CompressRuleTrait is equal to CompressType::N_BITS
	void ConsImp (const std::string& str_seq, Int2Type<CompressType::LOWERCASE>)
	{
		uint32_t N_len = 0;
		uint32_t N_start = 0;
		uint64_t pos = 0;
		char CharTemp = '0';
		sequence.resize (str_seq.size()*T::bit_count);
		uint32_t str_length = str_seq.size();
		for (uint32_t i = 0; i < str_length ; ++ i, pos+=T::bit_count)
		{
			if ( CharTemp != '0' )
			{
				if ( str_seq[i] == CharTemp )
				{
					++N_len;
					continue;
				}
				else
				{
					ExceptDeq.push_back (std::make_pair(N_start, std::make_pair (CharTemp, N_len)));
					N_start = 0, N_len = 0;
					CharTemp = '0';
				}
			}
			if ( CharNormal.find (str_seq[i]) != std::string::npos )
			{
				auto cr = TraitObject.table.find (str_seq[i]);
				for (auto ii = 0; ii!=T::bit_count; ++ii)
					sequence[pos+ii] = cr->second.second [ii];
				CharTemp = '0';
				continue;
			}
			else 
			{
				CharTemp = str_seq[i];
				if (CharExcept.find ( CharTemp ) == std::string::npos)
					CharExcept += CharTemp;
				N_start = i;
				++N_len;
			}
		}
		if (N_len!=0)
		{
			auto N_interv = std::make_pair ( N_start, std::make_pair (CharTemp, N_len) );
			ExceptDeq.push_back ( N_interv );
		}
	};
	void ConsImp (const std::string& str_seq, Int2Type<CompressType::N_BITS>)	//Implement for constructor taking a std::string
	{
		uint32_t N_len = 0;
		uint32_t N_start = 0;
		uint64_t pos = 0;;
		char CharTemp = '0';
		sequence.resize (str_seq.size()*T::bit_count);
		uint32_t str_length = str_seq.size();
		for (uint32_t i = 0; i < str_length ; ++ i, pos+=T::bit_count)
		{
			if ( CharTemp != '0' )
			{
				if ( str_seq[i] == CharTemp )
				{
					++N_len;
					continue;
				}
				else
				{
					ExceptDeq.push_back (std::make_pair(N_start, std::make_pair (CharTemp, N_len)));
					N_start = 0, N_len = 0;
					CharTemp = '0';
				}
			}
			if ( CharNormal.find (str_seq[i]) != std::string::npos )
			{
				auto cr = TraitObject.table.find (str_seq[i]);
				for (auto ii = 0; ii!=T::bit_count; ++ii)
					sequence[pos+ii] = cr->second.second [ii];
				CharTemp = '0';
				continue;
			}
			else 
			{
				CharTemp = str_seq[i];
				if (CharExcept.find ( CharTemp ) == std::string::npos)
					CharExcept += CharTemp;
				N_start = i;
				++N_len;
			}
		}
		if (N_len!=0)
		{
			auto N_interv = std::make_pair ( N_start, std::make_pair (CharTemp, N_len) );
			ExceptDeq.push_back ( N_interv );
		}
	};
/// @memberof NBitCompress
/// @brief taking care of some detail jobe for the constructor taking a std::string while the CompressType carried by the inputted CompressRuleTrait is equal to CompressType::N_BITS_MASK
	void ConsImp (const std::string& str_seq, Int2Type<CompressType::N_BITS_MASK>)	//Implement for constructor taking a std::string
	{
		uint32_t N_len = 0;
		uint32_t N_start = 0;
		uint64_t pos {};
		char CharTemp = '0';
		sequence.resize (str_seq.size()*T::bit_count);
		std::locale loc;
		uint32_t str_length = str_seq.size();
		for (uint32_t i = 0; i < str_length ; ++ i, pos+=T::bit_count)
		{
			if ( CharTemp != '0' )
			{
				if ( str_seq[i] == CharTemp )
				{
					++N_len;
					continue;
				}
				else if (std::islower (str_seq[i], loc) && std::islower (CharTemp, loc) )
				{
					++N_len;
					continue;
				}
				else if ( std::isupper ( CharTemp, loc ) )
				{
					ExceptDeq.push_back (std::make_pair(N_start, std::make_pair (CharTemp, N_len)));
					N_start = 0, N_len = 0;
					CharTemp = '0';
				}
				else 
				{
					MaskDeq.push_back (std::make_pair (N_start, 
						NBitCompress< CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > (str_seq.substr (N_start, N_len))));
					N_start = 0, N_len = 0;
					CharTemp = '0';
				}
			}
			if ( CharNormal.find (str_seq[i]) != std::string::npos )// )
			{
				auto cr = TraitObject.table.find (str_seq[i]);
				for (int ii = 0; ii!=T::bit_count; ++ii)
				{
					sequence[pos+ii] = cr -> second.second [ii];
				}
				CharTemp = '0';
				continue;
			}
			else 
			{
				CharTemp = str_seq[i];
				if ( std::isupper ( CharTemp, loc ) ) 
				{
					if (CharExcept.find ( CharTemp ) == std::string::npos) 
						CharExcept += CharTemp; 
				}
				N_start = i;
				++N_len;
			}
		}
		if (N_len!=0)
		{
			if ( std::isupper (CharTemp, loc ) )
				ExceptDeq.push_back (std::make_pair (N_start, std::make_pair (CharTemp, N_len)));
			else
			{
				MaskDeq.push_back (std::make_pair (N_start, 
						NBitCompress< CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > (str_seq.substr (N_start, N_len))));

			}
		}
std::cerr<<"sequence's sizeof "<<sizeof(sequence)<<std::endl;
std::cerr<<"ExceptDeq's sizeof "<<sizeof(ExceptDeq)<<std::endl;
std::cerr<<"MaskDeq's sizeof "<<sizeof(MaskDeq)<<std::endl;
	}
/// @memberof NBitCompress
/// @brief taking care of some detail jobe for the reverse function
	void reverse_impl (boost::dynamic_bitset<>& temp, 
			std::deque<std::pair<uint32_t, std::pair<char, uint32_t> > >& temp_ExceptDeq, 
			std::deque<std::pair<uint32_t, NBitCompress< CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > > >& temp_MaskDeq)
	{
		uint64_t len_bseq { sequence.size () };
		uint64_t len_char = len_bseq / T::bit_count;
		decltype(sequence) mask {sequence.size()}; // to mask certain bits
		for (int i=0; i!=T::bit_count; ++i)
			mask[i] = 1;
		int64_t len = reinterpret_cast<int64_t&> (len_bseq); // will need the minus value
		// human genome is 3095693983nt (10111000100001001000101010011111) 32bits,
		// thus 64 bits is more than enough and the sign bit will not be used (but uint32_t won't be enough)
		len-=T::bit_count;
		while (mask.any())  // mask 11000...000 will be the last mask
		{
			if (len >= 0)
				temp |= ( (mask&sequence) << len);
			else
				temp |= ( (mask&sequence) >> -len);
			mask<<=T::bit_count;
			len-=2*T::bit_count;
		}
		std::for_each ( ExceptDeq.begin(), ExceptDeq.end(),
						[&] (std::pair< uint32_t, std::pair<char, uint32_t> >& Q)
						{Q.first = len_char - Q.first - Q.second.second;
						 temp_ExceptDeq.push_front (Q);		}
						);
		std::for_each ( MaskDeq.begin(), MaskDeq.end(),
						[&] (std::pair< uint32_t, NBitCompress< CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > >& Q)
						{Q.first = len_char-Q.first-Q.second.GetSize()/T::bit_count;
						Q.second.reverse ();
						temp_MaskDeq.push_front (Q);		}
						);
	}
public:
/// @memberof NBitCompress
/// @brief holding the alphabet and conversion rule conveyed in the template parameter of CompressRuleTrait 
	T TraitObject;
/// @memberof NBitCompress
/// @brief holding the list of target character for the current encode/compression scheme
	std::string CharNormal;
/// @memberof NBitCompress
/// @brief holding the list of except character for the current encode/compression operation
	std::string CharExcept;
// Constructors	
	/// @brief default constructor
	/// @see TEST(NBitCompress, default_constructor_Normal) and
	/// @see TEST(NBitCompress, default_constructor_Mask) unit test examples
	NBitCompress ()	//default constructor
		: TraitObject (T())
		, CharNormal (TraitObject.alphabet), CharExcept ("")
		, sequence (0), ExceptDeq (), MaskDeq ()
	{
	}
	/// @brief constructor taking length information
	/// @see TEST(NBitCompress, constructor_with_length_Normal) and
	/// @see TEST(NBitCompress, constructor_with_length_Mask) unit test examples
	NBitCompress (uint64_t length)   // constructor from a fixed length
		: TraitObject (T())
		, CharNormal (TraitObject.alphabet), CharExcept ("")
		, sequence (length), ExceptDeq (), MaskDeq ()
	{}
	/// @brief constructor taking a std::string, prabably the most frequently used constructor
	/// @see TEST(NBitCompress, constructor_with_string_Normal) and
	/// @see TEST(NBitCompress, constructor_with_string_Mask) unit test examples
	NBitCompress (const std::string& str_seq)	//constructor with a string
		: TraitObject (T())
		, CharNormal (TraitObject.alphabet), CharExcept ("")
		, sequence (0), ExceptDeq (), MaskDeq ()
	{
		ConsImp (str_seq, Int2Type<T::type>() );
	}
	/// @brief copy constructor
	/// @see TEST(NBitCompress, copy_constructor_Normal) and
	/// @see TEST(NBitCompress, copy_constructor_Mask) unit test examples
	NBitCompress (const NBitCompress& other)	// copy constructor
		: TraitObject (other.TraitObject)
		, CharNormal (other.CharNormal), CharExcept (other.CharExcept)
		, sequence (other.sequence), ExceptDeq (other.ExceptDeq), MaskDeq (other.MaskDeq)
	{}
	/// @brief constructor with two elements
	/// @see TEST(NBitCompress, constructor_with_elements_Normal) unit test example
	NBitCompress (const boost::dynamic_bitset<>& db,	//constructor with two elements 
					const std::deque<std::pair<uint32_t, std::pair<char, uint32_t> > >& Ed)
		: TraitObject (T())
		, CharNormal (TraitObject.alphabet), CharExcept("")
		, sequence (db), ExceptDeq (Ed), MaskDeq ()
	{
		std::for_each ( ExceptDeq.begin(), ExceptDeq.end(),
						[&] (const std::pair < uint32_t, std::pair<char, uint32_t> > Q)
						{ if (CharExcept.find(Q.second.first)==std::string::npos)
							CharExcept += Q.second.first;}
						);
	}
	/// @brief constructor with two rvalue elements
	/// @see TEST(NBitCompress, constructor_with_rvalue_elements_Normal) unit test example
	NBitCompress (boost::dynamic_bitset<>&& db,	//constructor with two rvalue elements 
					std::deque<std::pair<uint32_t, std::pair<char, uint32_t> > >&& Ed)
		: TraitObject (T())
		, CharNormal (TraitObject.alphabet), CharExcept("")
		, MaskDeq ()
	{
		sequence.clear ();
		sequence.swap (db);
		ExceptDeq.clear();
		ExceptDeq.swap (Ed);
		std::for_each ( ExceptDeq.begin(), ExceptDeq.end(),
						[&] (const std::pair < uint32_t, std::pair<char, uint32_t> > Q)
						{ if (CharExcept.find(Q.second.first)==std::string::npos)
							CharExcept += Q.second.first;}
						);
	}
	/// @brief constructor with three elements
	/// @see TEST(NBitCompress, constructor_with_rvalue_elements_Mask) unit test example
	NBitCompress (const boost::dynamic_bitset<>& db,	//constructor with three elements 
					const std::deque<std::pair<uint32_t, std::pair<char, uint32_t> > >& Ed,
					const std::deque<std::pair<uint32_t, NBitCompress< CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > > >& Md)//std::string> >& Md)
		: TraitObject (T())
		, CharNormal (TraitObject.alphabet), CharExcept("")
		, sequence (db), ExceptDeq (Ed), MaskDeq (Md)
	{
		std::for_each ( ExceptDeq.begin(), ExceptDeq.end(),
						[&] (const std::pair < uint32_t, std::pair<char, uint32_t> > Q)
						{ if (CharExcept.find(Q.second.first)==std::string::npos)
							CharExcept += Q.second.first;}
						);
	}
	/// @brief constructor with three rvalue elements
	/// @see TEST(NBitCompress, constructor_with_rvalue_elements_Mask) unit test examples
	NBitCompress (boost::dynamic_bitset<>&& db,	//constructor with three rvalue elements 
					std::deque<std::pair<uint32_t, std::pair<char, uint32_t> > >&& Ed,
					std::deque<std::pair<uint32_t, NBitCompress< CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > > >&& Md)//std::string> >&& Md)
		: TraitObject (T())
		, CharNormal (TraitObject.alphabet), CharExcept("")
	{
		sequence.clear ();
		sequence.swap (db);
		ExceptDeq.clear (); 
		ExceptDeq.swap (Ed);
		MaskDeq.clear ();
		MaskDeq.swap (Md);
		std::for_each ( ExceptDeq.begin(), ExceptDeq.end(),
						[&] (const std::pair < uint32_t, std::pair<char, uint32_t> > Q)
						{ if (CharExcept.find(Q.second.first)==std::string::npos)
							CharExcept += Q.second.first;}
						);
	}
/// @breif move constructor
/// @see TEST(NBitCompress, move_constructor_Normal) and
/// @see TEST(NBitCompress, move_constructor_Mask) unit test examples
	NBitCompress (NBitCompress&& other)//move constructor
		: TraitObject (T())
		, CharNormal (TraitObject.alphabet)
	{
		assert (this != & other);
		CharExcept.swap (other.CharExcept);
		sequence.clear ();
		sequence.swap (other.sequence);
		ExceptDeq.clear ();
		ExceptDeq.swap (other.ExceptDeq);
		MaskDeq.clear ();
		MaskDeq.swap (other.MaskDeq);
	}
/// @breif move constructor
/// @see TEST(NBitCompress, Move_assignment_Normal) and
/// @see TEST(NBitCompress, Move_assignment_Mask) unit test examples
	NBitCompress& operator = (NBitCompress&& other) /// move assignment
	{
		assert (this != & other);
		sequence.clear ();
		sequence.swap (other.sequence);
		ExceptDeq.clear ();
		ExceptDeq.swap (other.ExceptDeq);
		MaskDeq.clear ();
		MaskDeq.swap (other.MaskDeq);
		CharNormal.swap (other.CharNormal);
		CharExcept.swap (other.CharExcept);
	}
/// @breif assignment
/// @see TEST(NBitCompress, assignment_Normal) and
/// @see TEST(NBitCompress, assignment_Mask) unit test examples
	NBitCompress& operator = (const NBitCompress& other) /// assignment operator
	{
		if (this == & other)
			return *this;
		TraitObject = other.TraitObject;
		CharNormal = other.CharNormal;
		CharExcept = other.CharExcept;
		sequence = other.sequence;
		ExceptDeq = other.ExceptDeq;
		MaskDeq = other.MaskDeq;
		return *this;
	}
// Get_functions
	uint64_t GetSize () const
	{
		return sequence.size ();
	}
	boost::dynamic_bitset<> GetSeq ()
	{
		return sequence;
	}
	std::deque<std::pair<uint32_t, std::pair <char, uint32_t> > > GetExceptDeq ()
	{ 
		return ExceptDeq;
	}
	std::deque<std::pair<uint32_t, NBitCompress< CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > > > GetMaskDeq ()
	{   
		return MaskDeq;
	}
// basic functions
	void Printer (void)
	{
//		std::cerr<<"convert sequence "<< this->MakeSeqString()<<std::endl;
		NBitPrinter();
	}
	void NBitPrinter (void)
	{
//		std::cerr<<"converted Nbit code "<<this->MakeNBitString()<<std::endl;
		std::cerr<<"ExceptDeq info "<<ExceptDeq.size()<<'\t'<<sizeof(ExceptDeq)<<std::endl;
		std::for_each ( ExceptDeq.begin(),
						ExceptDeq.end(),
						[] (const std::pair< uint32_t, std::pair <char, uint32_t> > Q)
						{std::cerr<<Q.first<<'\t'<<Q.second.first<<'\t'<<Q.second.second<<std::endl;}
						);
		if (MaskDeq.size()!=0)
		{
			std::cerr<<"MaskDeq info"<<MaskDeq.size()<<'\t'<<(char)sizeof(MaskDeq)<<std::endl;
			std::for_each ( MaskDeq.begin(),
							MaskDeq.end(),
							[] (std::pair< uint32_t, NBitCompress< CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > >& Q)//std::string> Q)
							{
//								std::cerr<<Q.first<<std::endl;
								std::cerr<<"infor for each element for a NBitCompress MaskDeq "<<sizeof(Q)<<std::endl;
								std::cerr<<Q.first<<std::endl;			
								Q.second.Printer();});
		}
	}
	/// @brief provide a conversion from the stored binary codes into its alphabetical characters form
	/// @see TEST(NBitCompress, MakeString_Normal) and
	/// @see TEST(NBitCompress, MakeString_Mask) unit test examples
	std::string MakeSeqString () //const
	{
		return GetSubStr (0, sequence.size() / T::bit_count);
	}
	void MakeSeqString (std::string& sstr)
	{
		this->GetSubStr (sstr, 0, sequence.size() / T::bit_count);//(0, sequence.size()>>1);
	}
	/// @brief provide a conversion from the stored binary codes into std::string form, with each of its bit converted into character form
	/// @see TEST(NBitCompress, MakeString_Normal) and
	/// @see TEST(NBitCompress, MakeString_Mask) unit test examples
	std::string MakeNBitString () const  
	{
		std::string buffer;
		buffer.reserve (sequence.size());
		boost::to_string (sequence, buffer);
		return buffer;
	}
	/// @brief provide the operator[] functionalitiy, which is substantially the same as the operator[] function of std::string object
	/// @see TEST(NBitCompress, middle_bracket_Normal) and
	/// @see TEST(NBitCompress, middle_bracket_Mask) unit test examples
	char operator [] (uint64_t p) const
	{
		//firstly execute operator[] to get first character and accordingly updates the value of the
		//Mask_count, Except_count, or Normal_count.  
		//When sstr[i] falls in CharNormal, the Normal_count is set to indicate 
		//how many characters belong to CharNormal is on the way.
		//When sstr[i] falls in CharExcept or mask rule, the Except_count and Mask_count are respectively set to indicate
		//how many except characters or how many mask characters are on the way.
		auto y = std::lower_bound ( ExceptDeq.begin(),
									ExceptDeq.end(), 
									p, 
									[](const std::pair<uint32_t, std::pair<char, uint32_t> >&pr, uint32_t oprnd)
									{ return (pr.first + pr.second.second -1) < oprnd;}
									);
//	if(! (sequence[p*T::bit_count] && sequence[p*T::bit_count+1]))
//	{
		for (auto i = y; i != y-2 ; --i)
		{
			if ( y-ExceptDeq.begin() == ExceptDeq.size() )
				break;
			if ( ( p >= i->first) && ( p < (i->first + i->second.second)  ) )
			{
				Except_count = i->first + i->second.second - p - 1;
				return i->second.first;
			}
			if (i == ExceptDeq.begin())
				break;
		}
		if (MaskDeq.size()!=0)
		{
			auto Y = std::lower_bound ( MaskDeq.begin(),
										MaskDeq.end(),
										p,
										[](const std::pair<uint32_t, NBitCompress< CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > >&pr, uint32_t oprnd)
										{ return (pr.first + pr.second.GetSize()/T::bit_count -1) < oprnd;}
										);
			for (auto j = Y; j != Y-2 ; --j)
			{
				if ( Y-MaskDeq.begin() == MaskDeq.size() )
					break;
				if ( ( p >= j->first) && ( p < (j->first + j->second.GetSize()/T::bit_count)  ) )
				{
					Mask_count = j->first + j->second.GetSize()/T::bit_count - p - 1;
					return j->second[p-(j->first)];
				}
				if (j == MaskDeq.begin())
					break;
			}
			if ( (Y-MaskDeq.begin() != MaskDeq.size())
				 && (y-ExceptDeq.begin() != ExceptDeq.size()) )
			{
				auto a = y->first-p-1;
				auto b = Y->first-p-1;
				if (a<b)
					Normal_count = a;//ReturnCount = a;
				else
					Normal_count = b;//ReturnCount = b;
			}
		}
//	}
		else
			Normal_count = y->first-p-1;
		return index_impl (p);
	}
	/// @brief provide the substr functionalitiy, which is substantially the same as the substr function of std::string object
	/// @see TEST(NBitCompress, Substr_Normal) and
	/// @see TEST(NBitCompress, Substr_Mask) unit test examples
	std::string GetSubStr (uint64_t start, uint64_t len) //const
	{
		std::string sstr;
		GetSubStr (sstr, start, len);
		return sstr;
	}
	void GetSubStr (std::string& sstr, uint64_t start, uint64_t len) //const
	{							   
		sstr.reserve(len);
		std::string mask_content;   
		sstr.resize (len);		  
		Mask_count = 0;
		Except_count = 0;
		Normal_count = 0;	
		//have the three counts initialized to zero whenever the GetSubStr fucntion is called
		for (uint64_t i = 0; i < len; ++ i)
		{   
			if (Normal_count != 0)	
			{   
			//indicating that a previous character belong to CharNormal is obtained, and the Normal_count accordingly indicates that
			//current sstr[i] is also belong to CharNormal and can be parsed without calling the full operat[] function
				sstr[i] = index_impl ( start + i );
				--Normal_count;
			}
			else if (Except_count != 0)	
			{
			//indicating that a previous character belong to CharNormal is obtained, and the Normal_count accordingly indicates that
			//indicating that a previous character belong to CharExcept is obtained, and the Except_count accordingly indicates that
			//current sstr[i] is also belong to CharExcept and can be parsed without calling the full operat[] function
				sstr[i] = sstr[i-1];
				--Except_count;
			}   
			else if ( (Mask_count!=0) && (mask_content.size()!=0) )	
			{		   
			// indicating that a previous character included by mask rule is obtained, and the Mask_count accordingly indicates that
			// current sstr[i] is also included by mask rule and can be parsed without calling the full operat[] function
			// direct value assign from the mask_content, a std::string object for storing the mask characters, is executed
				sstr[i] = mask_content[0];
				mask_content.erase(mask_content.begin());
				--Mask_count;   
			}   
			else if ( (Mask_count!=0) && (mask_content.size()==0) )
			{				   
			// indicating that a previous character included by mask rule is obtained, and the Mask_count accordingly indicates that
			// current sstr[i] is also included by mask rule and can be parsed without calling the full operat[] function
			// this is the second character coming right after that previous character, and an operation for obtaining the mask_content
			// is also executed.  Thus, the next comming mask character can be parsed by direct valu assign from the mask_content
			auto Y = std::lower_bound ( MaskDeq.begin(), MaskDeq.end(),	start+i,
				[](std::pair<uint32_t, NBitCompress< CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > >&pr, uint32_t oprnd)
					{ return (pr.first + pr.second.GetSize()/T::bit_count - 1) < oprnd;}	);
				if ( (i+start > Y->first) && (i+start< Y->first + Y->second.GetSize()/T::bit_count) )
					mask_content = (Y->second.MakeSeqString()).substr(i+start-Y->first);
				else if (i+start == Y->first-1)
					mask_content = (Y->second.MakeSeqString()).substr(1);
				std::cerr<<mask_content<<std::endl;
				sstr[i] = mask_content[0];
				mask_content.erase(mask_content.begin());
				--Mask_count;
			}
			else
			{
			// in the situation that all the three counts not able to give any parse hint, direct call the operator[] to parse the current character
				sstr[i] = this->operator[] ( start + i );								
			}												
		}													
	}
	friend class boost::serialization::access;
	template <typename Archive> 
	void serialize(Archive& ar, const unsigned int) 
	{
//		std::vector<boost::dynamic_bitset<>::block_type, boost::dynamic_bitset<>::allocator_type>* q;
//		to_block_range(sequence, &q);
//		size_t* qq;
//		to_block_range(sequence, &qq);
//		ar & (*q) & (*qq) & CharNormal & ExceptDeq;
//		if (T::type==CompressType::N_BITS_MASK)
//			ar & MaskDeq;
//		std::cerr<<"current T::type "<<T::type<<std::endl;
        serialize_base<Archive> (ar, Int2Type<T::type>());
    }
    
    template <typename Archive>
    void serialize_base (Archive& ar, Int2Type<CompressType::N_BITS_MASK> )
    {   
		std::vector<boost::dynamic_bitset<>::block_type, boost::dynamic_bitset<>::allocator_type>* q;
		to_block_range(sequence, &q);
		size_t* qq;
		to_block_range(sequence, &qq);
//		std::cerr<<"serializing bitset & exceptdeq & maskdeq "<<std::endl;
        ar & (*q) & (*qq) & CharNormal & ExceptDeq  & MaskDeq & TraitObject & ExceptDeq;
//		std::cerr<<"total seq length "<<Total<<std::endl;
		Total =0;
    }   
	template <typename Archive>
    void serialize_base (Archive& ar, Int2Type<CompressType::LOWERCASE> ) 
    {   
		std::vector<boost::dynamic_bitset<>::block_type, boost::dynamic_bitset<>::allocator_type>* q;
		to_block_range(sequence, &q);
		size_t* qq;
		to_block_range(sequence, &qq);
//		std::cerr<<"maskdeq in NBitCompress format"<<std::endl;
        ar & (*q) & (*qq) & ExceptDeq;
//		std::cerr<<"done a MaksDeq.second 's serializing "<<std::endl;
		Total+= sequence.size() / T::bit_count;
    }   
	
	template <typename Archive>
    void serialize_base (Archive& ar, Int2Type<CompressType::N_BITS> ) 
    {   
		std::vector<boost::dynamic_bitset<>::block_type, boost::dynamic_bitset<>::allocator_type>* q;
		to_block_range(sequence, &q);
		size_t* qq;
		to_block_range(sequence, &qq);
        ar & (*q) & (*qq)  & CharNormal & ExceptDeq & TraitObject;
    }   
// other functions
///	@biref Porivde the functionality for having the stored data revered in order
/// @see TEST(NBitCompress, reverse_Normal) and
/// @see TEST(NBitCompress, reverse_Mask) unit test examples
	NBitCompress& reverse ()
	{
		decltype(sequence) temp {sequence.size()}; /// to store the result
		std::deque < std::pair<uint32_t, std::pair <char, uint32_t> > > temp_ExceptDeq;
		std::deque < std::pair<uint32_t, NBitCompress< CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > > > temp_MaskDeq;//std::string> > temp_MaskDeq;
		reverse_impl ( temp, temp_ExceptDeq, temp_MaskDeq ); 
		sequence.swap (temp);
		ExceptDeq.swap (temp_ExceptDeq);
		if (MaskDeq.size()!=0)
			MaskDeq.swap (temp_MaskDeq);
		return *this;
	}
///	@biref Porivde the functionality for having each of the characters, in the stored data, complemented into its complement characters, e.g. A->T, G->C, and vice versa.
/// @see TEST(NBitCompress, complement_Normal) and
/// @see TEST(NBitCompress, complement_Mask) unit test examples
	NBitCompress& complement () // complement sequence
	{
		sequence.flip ();   // in-place flip
		return *this;
	}
///	@biref Porivde the functionality for both reverse and complement
/// @see TEST(NBitCompress, reverse_complement_Normal) and
/// @see TEST(NBitCompress, reverse_complement_Mask) unit test examples
	NBitCompress& reverse_complement ()/// reverse complement sequence
	{
		this->reverse ();
		this->complement ();
		return *this;
	}		   
///	@biref Porivde the functionality for provide an NBitCompress object, with its content equal to the reversed content of the current NBitCompress object
/// @see TEST(NBitCompress, reverse_copy_Normal) and
/// @see TEST(NBitCompress, reverse_copy_Mask) unit test examples
	NBitCompress reverse_copy ()
	{
		decltype(sequence) temp {sequence.size()}; //to store the result
		std::deque < std::pair<uint32_t, std::pair <char, uint32_t> > > temp_ExceptDeq;
		std::deque < std::pair<uint32_t, NBitCompress< CompressRuleTrait <T::bit_count, CompressType::LOWERCASE> > > > temp_MaskDeq;//std::string> > temp_MaskDeq;
		reverse_impl ( temp, temp_ExceptDeq, temp_MaskDeq );
		return NBitCompress (std::move(temp), std::move(temp_ExceptDeq), std::move(temp_MaskDeq) );
	}
///	@biref Porivde the functionality for provide an NBitCompress object, with its content equal to the complemented content of the current NBitCompress object
/// @see TEST(NBitCompress, complement_copy_Normal) and
/// @see TEST(NBitCompress, complement_copy_Mask) unit test examples
	NBitCompress complement_copy () // complement sequence
	{
		NBitCompress temp (sequence.flip(), ExceptDeq, MaskDeq);//{sequence.flip()};
		sequence.flip ();   // restore the state of sequence, things would be much easier if m_bits is not private...
		return temp;
	}   
///	@biref Porivde the functionality for provide an NBitCompress object, with its content equal to the reversed and complemented content of the current NBitCompress object
/// @see TEST(NBitCompress, reverse_complement_copy_Normal) and
/// @see TEST(NBitCompress, reverse_complement_copy_Mask) unit test examples
	NBitCompress reverse_complement_copy () // reverse complement sequence
	{   
		NBitCompress temp {*this};
		temp.reverse ();
		temp.complement ();
		return temp;
	}
	template <typename U>
	friend bool operator== (const NBitCompress<U>& q, const NBitCompress<U>& qq) ;
};

	template <typename T>
	bool operator== (const NBitCompress<T>& q, const NBitCompress<T>& qq) 
	{
		if ( (q.sequence==qq.sequence) && (q.ExceptDeq==qq.ExceptDeq) )
		//&& (q.MaskDeq==qq.MaskDeq) )
			return true;
		else
			return false;	
	}
#endif
