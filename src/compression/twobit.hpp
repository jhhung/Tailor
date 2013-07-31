///@file twobit.hpp
///@brief define twobit compression scheme 
///@author JHH Corp.
#ifndef TWOBIT_HPP_
#define TWOBIT_HPP_
#include <map>
#include <locale>
#include <algorithm>
#include <bitset>
#include <iostream>
#include <deque>
#include <string>
#include "boost/dynamic_bitset.hpp"
#include "boost/utility/binary.hpp"
//#include "../constant_def.hpp"

/*enum TwoBit_types
{
	Normal,
	Mask
};*/
class TWBImpl 
{
protected:
	std::string CharExcept;
	boost::dynamic_bitset<> sequence;
	std::deque < std::pair < uint32_t, std::pair <char, uint32_t> > > ExceptDeq;

public:
// constructor
	TWBImpl () 	/// default constructor
		: CharExcept (""), sequence (0), ExceptDeq ()
	{}//std::cerr<<"\ndefault constructing"<<std::endl;}  

	TWBImpl (uint64_t length) 	/// constructor from a fixed length
		: CharExcept (), sequence (length), ExceptDeq ()
	{}//std::cerr<<"\nconstructing with length"<<std::endl;}  

	TWBImpl (const TWBImpl& other)	/// copy constructor
		: sequence (other.sequence), ExceptDeq (other.ExceptDeq)
	{}//std::cerr<<"\ncopy constructing"<<std::endl;};

	TWBImpl (const boost::dynamic_bitset<>& db, const std::deque<std::pair<uint32_t, std::pair<char, uint32_t> > >& Ed)
		: sequence (db), ExceptDeq (Ed)
	{}//std::cerr<<"constructing with a bitset and a deque"<<std::endl;}

	TWBImpl (boost::dynamic_bitset<>&& db, std::deque<std::pair<uint32_t, std::pair<char, uint32_t> > >&& Ed)
	{
//		std::cerr<<"move constructing"<<std::endl;
		sequence.clear ();
		sequence.swap (db);
		ExceptDeq.clear();
		ExceptDeq.swap (Ed);
	}

	TWBImpl (TWBImpl&& other) ///move constructor
	{
//		std::cerr<<"move constructing with other twbseq"<<std::endl;
		assert (this != & other);
		sequence.clear ();
		sequence.swap (other.sequence);
		ExceptDeq.swap (other.ExceptDeq);
	}

	TWBImpl& operator = (TWBImpl&& other) /// assignment operator
	{
//		std::cerr<<"move assign"<<std::endl;
		assert (this != & other);
		sequence.clear ();
		sequence.swap (other.sequence);
		ExceptDeq.swap (other.ExceptDeq);
	}

	TWBImpl& operator = (const TWBImpl& other) /// assignment operator
	{
//		std::cerr<<"assign"<<std::endl;
		if (this == & other)
			return *this;
		sequence.clear ();
		sequence = other.sequence;
		ExceptDeq.clear();
		ExceptDeq = other.ExceptDeq;
		return *this;
	}

//basic functions
	void TWBPrinter (void)
	{
		std::cerr<<"converted twobit code "<<this->MakeTWBString()<<std::endl;
		std::cerr<<"ExceptDeq info"<<std::endl;
		std::for_each ( ExceptDeq.begin(), 
						ExceptDeq.end(),
						[] (const std::pair< uint32_t, std::pair <char, uint32_t> > Q)
						{std::cerr<<Q.first<<'\n'<<Q.second.first<<'\t'<<Q.second.second<<std::endl;}
						);
	}

	std::string MakeTWBString () const
	{
		std::string buffer;
		buffer.reserve (sequence.size()); 
		boost::to_string (sequence, buffer);
		return buffer;
	}

	char index_impl (uint64_t p) const  //Implement detail for operator[]
	{
		uint64_t pos {p<<1};
		if (pos >= sequence.size())
			throw std::out_of_range ("[Error] : TWBImpl : index longer than size");
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

	uint64_t GetSize () const
	{
		return sequence.size ();
	}

// Get_functions
	boost::dynamic_bitset<> GetSeq ()
	{
		return sequence;
	}

	std::deque<std::pair<uint32_t, std::pair <char, uint32_t> > > GetExceptDeq ()
	{
		return ExceptDeq;
	}
};

/// @brief define data type of TwoBitSequence, i.e. having the characters of 'A', 'C', 'G', and 'T'/'U' of the sequence encoded into two bit binary codes: \n 00, 01, 10, and 11.  Besides, data structure, such as ExceptDeque, is employed to keep track of the non-ACGT characters, such as N, which occationally appears in the sequence.
/// @tparam none-type parameter of int / enum of TwoBit_types, indicating whether mask functionality is involved   
template <int>
class TwoBitSequence
{};

/// @brief specialized form of the TwoBitSequence class, with the TwoBit_types specialized to be Normal, i.e. value of 0, indicating that no mask functionality is involved.
template <>
class TwoBitSequence < TwoBit_types::Normal >
	: public TWBImpl
{
private:
	std::string CharNormal;
	mutable uint32_t Normal_count, Except_count;
	void FillSeq ( char ChrIn, uint64_t pos )
	{
		switch (ChrIn)
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
		}	
	}
public:
// constructor
	TwoBitSequence () /// default constructor
		: TWBImpl (), CharNormal ("AaCcGgTtUu")
	{}//std::cerr<<"\ndefault constructing"<<std::endl;}  

	TwoBitSequence ( uint64_t length ) /// constructor from a fixed length
		: TWBImpl (length), CharNormal ("AaCcGgTtUu")
	{}//std::cerr<<"\nconstructing with length"<<std::endl;}  

	TwoBitSequence (const std::string& str_seq) /// constructor from a DNA string
		: TWBImpl (), CharNormal ("AaCcGgTtUu")
	{
//		std::cerr<<"\nconstructing with a string"<<std::endl;
		uint32_t len_seq = str_seq.size();
		uint64_t len_twb (len_seq << 1);
		uint32_t N_len = 0;
		uint32_t N_start = 0;
		uint64_t pos {};
		char CharTemp = '0';
		sequence.resize (len_twb);
		std::locale loc;
		for (uint32_t i = 0; i < len_seq ; ++ i)
		{
			pos = (i << 1);
			if ( CharTemp != '0' )
			{
				if ( str_seq[i] == CharTemp )
				{
					++N_len;
					continue;
				}
 				else //if ( std::isupper ( CharTemp, loc ) )
				{
					auto N_interv = std::make_pair ( N_start, std::make_pair ( CharTemp, N_len ) );
					ExceptDeq.push_back (N_interv);
					N_start = 0;
					N_len = 0;
					CharTemp = '0';
				}
			}
			if ( CharNormal.find (str_seq[i]) != std::string::npos )
			{
				FillSeq ( str_seq[i], pos );
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
			//std::cerr<<"str_seq[i] != CharTemp && CharTemp is UpperCase"<<std::endl;
			auto N_interv = std::make_pair ( N_start, std::make_pair (CharTemp, N_len) );
			ExceptDeq.push_back ( N_interv );
//std::cerr<<"Nmap updated with ExceptDeq["<<CharTemp<<"] push_backed with pair "<<N_start<<'\t'<<N_len<<" and has size of "<<ExceptDeq[CharTemp].size()<<std::endl;
		}
	}

	TwoBitSequence (const TwoBitSequence& other) /// copy constructor
		: CharNormal (other.CharNormal)
	{
		sequence = other.sequence;
		CharExcept = other.CharExcept;
		ExceptDeq = other.ExceptDeq;
	}//std::cerr<<"\ncopy constructing"<<std::endl;};

	TwoBitSequence (const boost::dynamic_bitset<>& db, const std::deque<std::pair<uint32_t, std::pair<char, uint32_t> > >& Ed)
		: TWBImpl ( db, Ed )
	{}//std::cerr<<"constructing with a bitset and a deque"<<std::endl;}

	TwoBitSequence (boost::dynamic_bitset<>&& db, std::deque<std::pair<uint32_t, std::pair<char, uint32_t> > >&& Ed)
		: TWBImpl ( db, Ed )
	{}//std::cerr<<"move constructing"<<std::endl;

	TwoBitSequence (TwoBitSequence&& other)//move constructor
	{
		assert (this != & other);
		sequence.clear ();
		sequence.swap (other.sequence);
		ExceptDeq.swap (other.ExceptDeq);
		CharExcept.swap (other.CharExcept);
		CharNormal.swap (other.CharNormal);
	}

	TwoBitSequence& operator = (TwoBitSequence&& other) /// assignment operator
	{
		assert (this != & other);
		sequence.clear ();
		sequence.swap (other.sequence);
		ExceptDeq.swap (other.ExceptDeq);
		CharExcept.swap (other.CharExcept);
		CharNormal.swap (other.CharNormal);
	}

	TwoBitSequence& operator = (const TwoBitSequence& other) /// assignment operator
	{
		if (this == & other)
			return *this;
		sequence.clear ();
		sequence = other.sequence;
		ExceptDeq.clear();
		ExceptDeq = other.ExceptDeq;
		CharExcept = other.CharExcept;
		CharNormal = other.CharNormal;
		return *this;
	}

//basic functions
	void Printer (void)
	{
		std::cerr<<"convert sequence "<< this->MakeSeqString()<<std::endl;
		std::cerr<<"ExceptDeq info "<<std::endl;
		std::for_each ( ExceptDeq.begin(), ExceptDeq.end(), 
						[] (const std::pair<uint32_t, std::pair < char, uint32_t> >& Q)
						{ std::cerr<<Q.second.first<<'\t'<<Q.first<<'\t'<<Q.second.second<<std::endl;}
						);
	}

	std::string MakeSeqString () const
	{
		return this->GetSubStr (0, sequence.size()>>1);
	}

	char operator [] (uint64_t p) const
	{
		auto y = std::lower_bound ( ExceptDeq.begin(), 
									ExceptDeq.end(), 
									p, 
									[](const std::pair<uint32_t, std::pair<char, uint32_t> >&pr, uint32_t oprnd)
									{ return (pr.first + pr.second.second -1) < oprnd;}
									);
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
		if (y-ExceptDeq.begin() != ExceptDeq.size()) 
			Normal_count = y->first-p-1;
		return index_impl (p);
	}

	std::string GetSubStr (uint64_t start, uint64_t len) const
	{
		std::string sstr;
		std::string mask_content;
		sstr.resize (len);
		uint32_t Except_count = 0;
		uint32_t Normal_count = 0;
		std::locale loc;
		for (uint64_t i = 0; i < len; ++ i)
		{
			if (Normal_count != 0)
			{
				sstr[i] = index_impl ( start + i );
				--Normal_count;
			}
			else if (Except_count != 0)
			{
				sstr[i] = sstr[i-1];
				--Except_count;
			}
			else
			{
				auto z = this->operator[] ( start + i );
				sstr[i] = z;
			}
		}
		return sstr;
	}

//other functions
	TwoBitSequence& reverse ()
	{
		uint64_t len_twb { sequence.size () };
		uint32_t len_seq = len_twb>>1;
		decltype(sequence) temp {len_twb}; /// to store the result
		decltype(sequence) mask {len_twb}; /// to mask certain bits
		mask[0] = mask[1] = 1;		  /// initialized
		int64_t len = reinterpret_cast<int64_t&> (len_twb);	/// will need the minus value
		/// human genome is 3095693983nt (10111000100001001000101010011111) 32bits,
		/// thus 64 bits is more than enough and the sign bit will not be used (but uint32_t won't be enough)
		len-=2;
		while (mask.any())  /// mask 11000...000 will be the last mask
		{
			if (len >= 0)
				temp |= ( (mask&sequence) << len);
			else
				temp |= ( (mask&sequence) >> -len);
			mask<<=2;   /// left shift to retrieve the next nt
			len-=4;	 /// two bits shifted from each side
		}
		sequence.swap (temp);

		std::deque < std::pair<uint32_t, std::pair <char, uint32_t> > > temp_ExceptDeq;
		std::for_each ( ExceptDeq.begin(), ExceptDeq.end(),
						[&] (std::pair< uint32_t, std::pair<char, uint32_t> >& Q)
						{Q.first = len_seq - Q.first - Q.second.second;
						 temp_ExceptDeq.push_front (Q);}
						);
		ExceptDeq.swap ( temp_ExceptDeq );
		return *this;
	}

	TwoBitSequence& complement () /// complement sequence
	{
		sequence.flip ();   /// in-place flip
		return *this;
	}
	
	TwoBitSequence& reverse_complement () /// reverse complement sequence
	{
		this->reverse ();
		this->complement ();
		return *this;
	}

	TwoBitSequence reverse_copy ()	/// reverse sequence
	{
		uint64_t len_twb  { sequence.size () };
		uint32_t len_seq = len_twb>>1;
		decltype(sequence) temp {len_twb}; /// to store the result
		decltype(sequence) mask {len_twb}; /// to mask certain bits
		mask[0] = mask[1] = 1;		  /// initialized
		int64_t len = reinterpret_cast<int64_t&> (len_twb);	/// will need the minus part
		/// human genome is 3095693983nt (10111000100001001000101010011111) 32bits,
		/// thus 64 bits is more than enough and the sign bit will not be used (but uint32_t won't be enough)
		len-=2;
		while (mask.any())  /// mask 11000...000 will be the last mask
		{
			if (len >= 0)
				temp |= ( (mask&sequence) << len);
			else
				temp |= ( (mask&sequence) >> -len);
			mask<<=2;   /// left shift to retrieve the next nt
			len-=4;	 /// two bits shifted from each side
		}

		std::deque < std::pair<uint32_t, std::pair <char, uint32_t> > > temp_ExceptDeq;
		std::for_each ( ExceptDeq.begin(), ExceptDeq.end(),
						[&] (std::pair< uint32_t, std::pair<char, uint32_t> >& Q)
						{Q.first = len_seq - Q.first - Q.second.second;
						 temp_ExceptDeq.push_front (Q);}
						);
		return TwoBitSequence (std::move(temp), std::move(temp_ExceptDeq) );
	}

	TwoBitSequence complement_copy () /// complement sequence
	{
		TwoBitSequence temp (sequence.flip(), ExceptDeq);//{sequence.flip()};
		sequence.flip ();   /// restore the state of sequence, things would be much easier if m_bits is not private...
		return temp;
	}

	TwoBitSequence reverse_complement_copy () /// reverse complement sequence
	{
		TwoBitSequence temp {*this};
		temp.reverse ();
		temp.complement ();
		return temp;
	}
};

/// @brief specialized form of the TwoBitSequence class, with the TwoBit_types specialized to be Mask, i.e. value of 0, indicating that current TwoBitSequence object includes mask handling capability.
template <>
class TwoBitSequence < TwoBit_types::Mask >
	: public TWBImpl
{
private:
	std::string CharNormal;
	std::deque < std::pair < uint32_t, std::string > > MaskDeq;
	mutable uint32_t Normal_count, Except_count, Mask_count;//skip_count, NonN_count;
	void FillSeq ( char ChrIn, uint64_t pos )
	{
		switch (ChrIn)
		{
			case 'A':
			//case 'a':
				sequence[pos]   = 0;
				sequence[pos+1] = 0;
				break;

			case 'C':
			//case 'c':
				sequence[pos]   = 1;
				sequence[pos+1] = 0;
				break;

			case 'G':
			//case 'g':
				sequence[pos]   = 0;
				sequence[pos+1] = 1;
				break;

			case 'T':
			//case 't':
			case 'U':
			//case 'u':
				sequence[pos]   = 1;
				sequence[pos+1] = 1;
				break;
		}	
	}

public:
//Constructors
	TwoBitSequence () /// default constructor
		: TWBImpl (), CharNormal ("ACGTU")
	{}//std::cerr<<"\ndefault constructing"<<std::endl;}  

	TwoBitSequence ( uint64_t length )
		: TWBImpl (length) 	/// constructor from a fixed length
		, CharNormal ("ACGTU")
	{}//std::cerr<<"\nconstructing with length"<<std::endl;}  

	TwoBitSequence (const std::string& str_seq) /// constructor from a DNA string
		: TWBImpl () 
		, CharNormal ("ACGTU")
		, MaskDeq ()
		, Normal_count ()
		, Except_count ()
		, Mask_count ()
	{
		uint32_t len_seq = str_seq.size();
		uint64_t len_twb (len_seq << 1);
		uint32_t N_len = 0;
		uint32_t N_start = 0;
		uint64_t pos {};
		char CharTemp = '0';
		sequence.resize (len_twb);
		std::locale loc;
		for (uint32_t i = 0; i < len_seq ; ++ i)
		{
			pos = (i << 1);
			if ( CharTemp != '0' )
			{
				if ( str_seq[i] == CharTemp )
				{
					++N_len;
					continue;
				}
				else if ( std::islower (str_seq[i], loc) && std::islower (CharTemp, loc) )
				{
					++N_len;
					continue;
				}
 				else if ( std::isupper ( CharTemp, loc ) )
				{
					auto N_interv = std::make_pair ( N_start, std::make_pair ( CharTemp, N_len ) );
					ExceptDeq.push_back (N_interv);
					N_start = 0;
					N_len = 0;
					CharTemp = '0';
				}
				else
				{
					std::string mask_element (str_seq.substr(N_start, N_len));//(sequence[N_start<<1], N_len);
					auto M_interv = std::make_pair ( N_start, mask_element );
					MaskDeq.push_back ( M_interv ); 
					N_start = 0;
					N_len = 0;
					CharTemp = '0';
				}
			}
			if ( CharNormal.find (str_seq[i]) != std::string::npos )
			{
				FillSeq ( str_seq[i], pos );
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
			{	
				auto N_interv = std::make_pair ( N_start, std::make_pair (CharTemp, N_len) );
				ExceptDeq.push_back ( N_interv );
			}
			else
			{
				std::string mask_element (str_seq.substr(N_start, N_len));//(sequence[N_start<<1], N_len);
				auto M_interv = std::make_pair ( N_start, mask_element );
				MaskDeq.push_back ( M_interv ); 
			}
		}
	}

	TwoBitSequence (const TwoBitSequence& other) /// copy constructor
		: CharNormal (other.CharNormal)
		, MaskDeq (other.MaskDeq)
		, Normal_count () 
		, Except_count ()
		, Mask_count ()
	{
		sequence = other.sequence;
		CharExcept = other.CharExcept;
		ExceptDeq = other.ExceptDeq;
	}//std::cerr<<"\ncopy constructing"<<std::endl;};

	TwoBitSequence (const boost::dynamic_bitset<>& db, 
					const std::deque<std::pair<uint32_t, std::pair<char, uint32_t> > >& Ed,// = std::deque<std::pair<uint32_t, int> > (0)  )
					const std::deque<std::pair<uint32_t, std::string> >& Md)
		: TWBImpl ( db, Ed )
		, MaskDeq(Md)
		, Normal_count ()
		, Except_count ()
		, Mask_count ()
	{}//std::cerr<<"constructing with a bitset and a deque"<<std::endl;}

	TwoBitSequence (boost::dynamic_bitset<>&& db, 
					std::deque<std::pair<uint32_t, std::pair<char, uint32_t> > >&& Ed,
					std::deque<std::pair<uint32_t, std::string> >&& Md)
		: TWBImpl ( db, Ed )
	{
		MaskDeq.swap (Md);
	}//std::cerr<<"move constructing"<<std::endl;

	TwoBitSequence (TwoBitSequence&& other)//move constructor
	{
//		std::cerr<<"move constructing with other twbseq"<<std::endl;
		assert (this != & other);
		sequence.swap (other.sequence);
		ExceptDeq.swap (other.ExceptDeq);
		MaskDeq.swap (other.MaskDeq);
		CharExcept.swap (other.CharExcept);
		CharNormal.swap (other.CharNormal);
	}

	TwoBitSequence& operator = (TwoBitSequence&& other) /// assignment operator
	{
//		std::cerr<<"move assign"<<std::endl;
		assert (this != & other);
		sequence.swap (other.sequence);
		ExceptDeq.swap (other.ExceptDeq);
		MaskDeq.swap (other.MaskDeq);
		CharExcept.swap (other.CharExcept);
		CharNormal.swap (other.CharNormal);
	}

	TwoBitSequence& operator = (const TwoBitSequence& other) /// assignment operator
	{
//		std::cerr<<"assign"<<std::endl;
		if (this == & other)
			return *this;
		sequence = other.sequence;
		ExceptDeq = other.ExceptDeq;
		MaskDeq = other.MaskDeq;
		CharExcept = other.CharExcept;
		CharNormal = other.CharNormal;
		return *this;
	}

//basic functions
	void Printer (void)
	{
		std::cerr<<"convert sequence "<< this->MakeSeqString()<<std::endl;
		std::cerr<<"MaskDeq info"<<std::endl;
		std::for_each ( MaskDeq.begin(),
						MaskDeq.end(),
						[] (const std::pair<uint32_t, std::string>& Q)
						{std::cerr<<Q.first<<'\t'<<Q.second<<std::endl;}
						);
		std::cerr<<"ExceptDeq info"<<std::endl;
		std::for_each ( ExceptDeq.begin(), 
						ExceptDeq.end(),
						[] (const std::pair< uint32_t, std::pair <char, uint32_t> > Q)
						{std::cerr<<Q.first<<'\n'<<Q.second.first<<'\t'<<Q.second.second<<std::endl;}
						);
	}

	std::string MakeSeqString () const
	{
		return this->GetSubStr (0, sequence.size()>>1);
	}

	char operator [] (uint64_t p) const
	{
		auto y = std::lower_bound ( ExceptDeq.begin(), 
									ExceptDeq.end(), 
									p, 
									[](const std::pair<uint32_t, std::pair<char, uint32_t> >&pr, uint32_t oprnd)
									{ return (pr.first + pr.second.second -1) < oprnd;}
									);
		//		std::cerr<<"\ny-ExceptDeq.begin()"<<y-ExceptDeq.begin()<<std::endl;
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
		auto Y = std::lower_bound ( MaskDeq.begin(), 
									MaskDeq.end(), 
									p, 
									[](const std::pair<uint32_t, std::string >&pr, uint32_t oprnd)
									{ return (pr.first+pr.second.size()-1) < oprnd;}
									);
		for (auto j = Y; j != Y-2 ; --j)
		{
			if ( Y-MaskDeq.begin() == MaskDeq.size() )
				break;
			if ( ( p >= j->first) && ( p < (j->first + j->second.size())  ) )
			{
				Mask_count = j->first + j->second.size() - p - 1;
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
		return index_impl (p);
	}

	std::string GetSubStr (uint64_t start, uint64_t len) const
	{
		std::string sstr;
		std::string mask_content;
		sstr.resize (len);
		uint32_t Mask_count = 0;
		uint32_t Except_count = 0;
		uint32_t Normal_count = 0;
		std::locale loc;
		for (uint64_t i = 0; i < len; ++ i)
		{
			if (Normal_count != 0)
			{
				sstr[i] = index_impl ( start + i );
				--Normal_count;
			}
			else if (Except_count != 0)
			{
				sstr[i] = sstr[i-1];
				--Except_count;
			}
			else if ( (Mask_count!=0) && (mask_content.size()!=0) )
			{			
				sstr[i] = mask_content[0];
				mask_content.erase(mask_content.begin());
				--Mask_count;
			}
			else if ( (Mask_count!=0) && (mask_content.size()==0) )
			{
				std::for_each ( MaskDeq.begin(), 
								MaskDeq.end(),
								[&] (const std::pair <uint32_t, std::string>& Q)
								{ if ( i-1 == Q.first )  mask_content = Q.second.substr (1); }
								);
				sstr[i] = mask_content[0];
				mask_content.erase(mask_content.begin());
				--Mask_count;
			}
			else
			{
				auto z = this->operator[] ( start + i );
				sstr[i] = z;
			}
		}
		return sstr;
	}

	std::deque<std::pair<uint32_t, std::string> > GetMaskDeq ()
	{
		return MaskDeq;
	}

// other functions
	TwoBitSequence& reverse ()
	{
		uint64_t len_twb { sequence.size () };
		uint32_t len_seq = len_twb>>1;
		decltype(sequence) temp {len_twb}; /// to store the result
		decltype(sequence) mask {len_twb}; /// to mask certain bits
		mask[0] = mask[1] = 1;		  /// initialized
		int64_t len = reinterpret_cast<int64_t&> (len_twb);	/// will need the minus value
		/// human genome is 3095693983nt (10111000100001001000101010011111) 32bits,
		/// thus 64 bits is more than enough and the sign bit will not be used (but uint32_t won't be enough)
		len-=2;
		while (mask.any())  /// mask 11000...000 will be the last mask
		{
			if (len >= 0)
				temp |= ( (mask&sequence) << len);
			else
				temp |= ( (mask&sequence) >> -len);
			mask<<=2;   /// left shift to retrieve the next nt
			len-=4;	 /// two bits shifted from each side
		}
		sequence.swap (temp);

		std::deque < std::pair<uint32_t, std::pair <char, uint32_t> > > temp_ExceptDeq;
		std::for_each ( ExceptDeq.begin(), ExceptDeq.end(),
						[&] (std::pair< uint32_t, std::pair<char, uint32_t> >& Q)
						{Q.first = len_seq - Q.first - Q.second.second;
						 temp_ExceptDeq.push_front (Q);}
						);
		ExceptDeq.swap ( temp_ExceptDeq );

		std::deque < std::pair<uint32_t, std::string> > temp_MaskDeq;
		std::for_each ( MaskDeq.begin(), MaskDeq.end(),
						[&] (std::pair< uint32_t, std::string>& Q)
						{Q.first = len_seq-Q.first-Q.second.size();
						std::reverse (Q.second.begin(), Q.second.end());
						temp_MaskDeq.push_front (Q);		}
						);
		MaskDeq.swap ( temp_MaskDeq );
		return *this;
	}

	TwoBitSequence& complement () /// complement sequence
	{
		sequence.flip ();   /// in-place flip
		return *this;
	}
	
	TwoBitSequence& reverse_complement () /// reverse complement sequence
	{
		this->reverse ();
		this->complement ();
		return *this;
	}

	TwoBitSequence reverse_copy ()	/// reverse sequence
	{
		uint64_t len_twb  { sequence.size () };
		uint32_t len_seq = len_twb>>1;
		decltype(sequence) temp {len_twb}; /// to store the result
		decltype(sequence) mask {len_twb}; /// to mask certain bits
		mask[0] = mask[1] = 1;		  /// initialized
		int64_t len = reinterpret_cast<int64_t&> (len_twb);	/// will need the minus part
		/// human genome is 3095693983nt (10111000100001001000101010011111) 32bits,
		/// thus 64 bits is more than enough and the sign bit will not be used (but uint32_t won't be enough)
		len-=2;
		while (mask.any())  /// mask 11000...000 will be the last mask
		{
			if (len >= 0)
				temp |= ( (mask&sequence) << len);
			else
				temp |= ( (mask&sequence) >> -len);
			mask<<=2;   /// left shift to retrieve the next nt
			len-=4;	 /// two bits shifted from each side
		}

		std::deque < std::pair<uint32_t, std::pair <char, uint32_t> > > temp_ExceptDeq;
		std::for_each ( ExceptDeq.begin(), ExceptDeq.end(),
						[&] (std::pair< uint32_t, std::pair<char, uint32_t> >& Q)
						{Q.first = len_seq - Q.first - Q.second.second;
						 temp_ExceptDeq.push_front (Q);}
						);

		std::deque < std::pair<uint32_t, std::string> > temp_MaskDeq;
		std::for_each ( MaskDeq.begin(), MaskDeq.end(),
						[&] (std::pair< uint32_t, std::string>& Q)
						{Q.first = len_seq-Q.first-Q.second.size();
						std::reverse (Q.second.begin(), Q.second.end());
						temp_MaskDeq.push_front (Q);		}
						);

		return TwoBitSequence (std::move(temp), std::move(temp_ExceptDeq), std::move(temp_MaskDeq) );
	}

	TwoBitSequence complement_copy () /// complement sequence
	{
		TwoBitSequence temp (sequence.flip(), ExceptDeq, MaskDeq);//{sequence.flip()};
		sequence.flip ();   /// restore the state of sequence, things would be much easier if m_bits is not private...
		return temp;
	}

	TwoBitSequence reverse_complement_copy () /// reverse complement sequence
	{
		TwoBitSequence temp {*this};
		temp.reverse ();
		temp.complement ();
		return temp;
	}
};
#endif
