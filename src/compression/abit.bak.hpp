#ifndef ABIT_HPP_
#define ABIT_HPP_
#include <map>
#include <locale>
#include <algorithm>
#include <bitset>
#include <iostream>
#include <deque>
#include <string>
#include "boost/dynamic_bitset.hpp"
#include "boost/utility/binary.hpp"
#include "../constant_def.hpp"


template<typename IndexType, typename ContainerType >
class QSeqTable
{
public:
	
	IndexType Qseq_s, Qseq_e, max;
	typedef typename std::deque < std::pair < IndexType, ContainerType > >::iterator DequeIterType;
	DequeIterType Qseq_i;
	
	QSeqTable()
		:Container(), Qseq_s(0), Qseq_e(0), Qseq_i(Container.begin()), max(0)
	{}
	QSeqTable(IndexType m)
		:Container(), Qseq_s(0), Qseq_e(0), Qseq_i(Container.begin()), max(m)
	{}
	
	void insert( const IndexType &index, const ContainerType &value)
	{
		auto data = std::make_pair ( index, value );
		Container.push_back (data);
		if( Qseq_i != Container.begin() )  Qseq_i = Container.begin();
	}
	inline bool next(){
		if(Qseq_e != max){
			++Qseq_i;
			resetValue();
			return true;
		}
		return false;
	}
	
	inline void bsearch(IndexType pos){
		Qseq_i = std::lower_bound(Container.begin(), Container.end(), pos, 
			[](const std::pair< IndexType, ContainerType >&pr, uint32_t oprnd)
			{ return pr.first <= oprnd;}
		);
		if( !isBegin() ) --Qseq_i;
		resetValue();
	}
	inline void fsearch(IndexType pos){
		if(pos >= Qseq_e && Qseq_e != max){
			next();
			if( pos < Qseq_e ) return true;
		}
		bsearch(pos);
	}
	
	bool isBegin() {return (Container.begin() == Qseq_i);}
	bool isEnd() {return (Container.end() == Qseq_i);}
	void setMax(IndexType m){max = m;}
	DequeIterType begin(){return Container.begin();}
	DequeIterType end(){return Container.end();}
	bool empty(){return Container.empty();}
	
private:
	std::deque < std::pair < IndexType, ContainerType > > Container;
	
	inline void resetValue()
	{
		Qseq_s = (*Qseq_i).first;
		if( (Qseq_i+1) == Container.end() )
			Qseq_e = max;
		else
			Qseq_e = (*(Qseq_i+1)).first;
	}
	
};




class ABmpl
{
protected:
	std::string CharExcept;
	
	uint32_t len_seq;
	boost::dynamic_bitset<> sequence;
	
	uint32_t Qseq_s, Qseq_e;
	std::deque < std::pair < uint32_t, std::tuple <uint64_t, char, uint32_t> > >::iterator Qseq_i;
	std::deque < std::pair < uint32_t, std::tuple <uint64_t, char, uint32_t> > > ExceptDeq;
	
	QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > EscapeChar;
	
public:
	ABmpl () 	/// default constructor
		: CharExcept (""), sequence (0), ExceptDeq (), Qseq_s(0), Qseq_e(0), Qseq_i(ExceptDeq.begin()), len_seq(0)
	{}//std::cerr<<"\ndefault constructing"<<std::endl;}
	
	std::string MakeSeqString () const
	{
		std::string buffer;
		buffer.reserve (sequence.size()); 
		boost::to_string (sequence, buffer);
		return buffer;
	}
	char index_impl (uint32_t p) const  //Implement detail for operator[]
	{
		if (p >= len_seq)
			throw std::out_of_range ("[Error] : ABmpl : index longer than size");
		uint64_t pos {p<<1};
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
	
	void Printer (void)
	{
		std::cerr<<"ExceptDeq info"<<std::endl;
		std::for_each ( EscapeChar.begin(), 
			EscapeChar.end(),
			[] (const std::pair< uint32_t, std::tuple <uint64_t, char, uint32_t> > Q)
			{
				std::cerr<<Q.first<<'\n'
				<< std::get<0>(Q.second) << '\t'
				<< std::get<1>(Q.second) << '\t'
				<< std::get<2>(Q.second)
				<<std::endl;
			}
		);
	}
	uint32_t size()
	{
		return len_seq;
	}
	uint32_t sizeSeq()
	{
		return sequence.size();
	} 
};

class ABSequence 
	: public ABmpl
{
private:
	bool FillSeq ( char ChrIn, uint64_t pos, uint32_t length = 1)
	{
		bool first, second;
		switch (ChrIn)
		{
			case 'A':
			case 'a':
				first  = 0;
				second = 0;
				break;

			case 'C':
			case 'c':
				first  = 1;
				second = 0;
				break;

			case 'G':
			case 'g':
				first  = 0;
				second = 1;
				break;

			case 'T':
			case 't':
			case 'U':
			case 'u':
				first  = 1;
				second = 1;
				break;
			case 'N':
				return false;
				break;
		}
		for(uint32_t i = 0; i < length; ++i)
		{
			sequence[pos + (i << 1)] = first;
			sequence[pos + (i << 1) + 1] = second;
		}
		return true;
	}
public:
	uint32_t bt;
	uint32_t ft;
	ABSequence ()
	{}
	ABSequence (const std::string& str_seq) /// default constructor
		: ABmpl (), bt(0), ft(0)
	{
		
		uint32_t len_record = 1;
		len_seq = str_seq.size();
		uint64_t pos(0);
		sequence.resize (len_seq << 1);
		uint32_t i = 1;
		
		char tmp_repeat_char = str_seq[0];
		uint32_t tmp_repeat_length = 1;
		uint32_t tmp_repeat_limit = 100; // > 5 will be record
		
		EscapeChar.setMax(len_seq);

		for (; i <= len_seq ; ++ i)
		{
			if(str_seq[i] == tmp_repeat_char && i != (len_seq) )
			{
				++tmp_repeat_length;
			}
			else
			{
				if(tmp_repeat_length < tmp_repeat_limit)
				{
					//too short to save to table
					if ( FillSeq ( str_seq[ i - tmp_repeat_length ], pos, tmp_repeat_length) )
					{
						// success
						pos += (tmp_repeat_length << 1);
					}
					else
					{
						// this is N, save to table
						EscapeChar.insert((i-tmp_repeat_length), std::make_tuple ( pos, tmp_repeat_char, tmp_repeat_length ) );
					}// if ( FillSeq ( str_seq[ i - tmp_repeat_length ], pos, tmp_repeat_length) )
				}//if (tmp_repeat_times <= tmp_repeat_limit)
				else
				{
					// enough length, save to table
					EscapeChar.insert((i-tmp_repeat_length), std::make_tuple ( pos, tmp_repeat_char, tmp_repeat_length ) );
				}//else (tmp_repeat_times <= tmp_repeat_limit)
				
				tmp_repeat_char = str_seq[i];
				tmp_repeat_length = 1;
				
			}// else (str_seq[i] == tmp_repeat_char)
		}// for (uint32_t i = 1; i < len_seq ; ++ i)

		sequence.resize (pos);
	}// ABSequence (const std::string& str_seq)
	
  inline char operator[] (uint32_t pos)
  {
  	if(EscapeChar.empty())  return index_impl(pos);
  	
  	uint32_t &start = EscapeChar.Qseq_s;
  	uint32_t &end = EscapeChar.Qseq_e;
  	
  	if(EscapeChar.isBegin() && pos < start)  return index_impl(pos);
  	
  	if( !(pos >= start && pos < end) )  EscapeChar.bsearch(pos);
  	
  	if(EscapeChar.isBegin() && pos < start)  return index_impl(pos);

  	uint32_t repeat_end_pos(  (*EscapeChar.Qseq_i).first + std::get<2>( (*EscapeChar.Qseq_i).second ) - 1  );
  	
	  if(pos <= repeat_end_pos)
	  	return std::get<1>( (*EscapeChar.Qseq_i).second );
	  else
	  	return index_impl(  pos - repeat_end_pos + ((std::get<0>( (*EscapeChar.Qseq_i).second )) >> 1) -1  );
  	
  }
	
};









#endif