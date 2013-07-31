#ifndef ABIT_HPP_
#define ABIT_HPP_
#include <map>
#include <locale>
#include <algorithm>
#include <bitset>
#include <iostream>
#include <deque>
#include <string>
#include <set>
#include "boost/dynamic_bitset.hpp"
#include "boost/utility/binary.hpp"
#include "../constant_def.hpp"
#include <cctype>

template<typename IndexType, typename ContainerType >
class QSeqTable
{
public:
	
	mutable IndexType Qseq_s, Qseq_e, max;
	//typedef typename std::deque < std::pair < IndexType, ContainerType > >::iterator DequeIterType;
	typedef typename std::map < IndexType, ContainerType	>::iterator DequeIterType;
	mutable DequeIterType Qseq_i;
	uint32_t f,b;
	QSeqTable()
		:Container(), Qseq_s(0), Qseq_e(0), Qseq_i(Container.begin()), max(0),f(0),b(0)
	{}
	QSeqTable(IndexType m)
		:Container(), Qseq_s(0), Qseq_e(0), Qseq_i(Container.begin()), max(m),f(0),b(0)
	{}
	QSeqTable(QSeqTable &QseqTableIn)
		:Container(QseqTableIn.Container), Qseq_s(0), Qseq_e(0), Qseq_i(Container.begin()), max(QseqTableIn.max),f(0),b(0)
	{}
	////
	inline void insert( const IndexType &index, const ContainerType &value)
	{
		// use make_pair is slower
		//auto data = std::make_pair ( index, value );
		//Container.push_back (data);
		//Container.insert (data);
		Container.insert ( {{index,value}});
		if( Qseq_i != Container.begin() )	Qseq_i = Container.begin();
	}
	inline void next()
	{
		if(Qseq_e != max){
			++Qseq_i;
			resetValue();
		}
	}
	
	inline void lower_bound_search(const IndexType &pos)
	{
		if( !(pos >= Qseq_s && pos < Qseq_e))
		{
			//++b;
			Qseq_i = Container.lower_bound(pos);
			if( !isBegin() && Qseq_i->first != pos)
				--Qseq_i;
			
			/*
			Qseq_i = std::lower_bound(Container.begin(), Container.end(), pos, 
				[](const std::pair< IndexType, ContainerType >&pr, uint32_t oprnd)
				{ return pr.first < oprnd;}
			);
			*/
			
			resetValue();
		}
	}
	inline void jump_lower_bound(const IndexType &pos)
	{
		/*
		++Qseq_i;
		if((*Qseq_i).first == pos)
		{
			resetValue();
		}
		else if( !(pos >= Qseq_s && pos < Qseq_e))
		{
			Qseq_i = Container.lower_bound(pos);
			resetValue();
		}
		else
		{
			--Qseq_i;
		}
		*/
		Qseq_i = Container.lower_bound(pos);
		if( !isBegin() && Qseq_i->first != pos)
				--Qseq_i;
		resetValue();
		//
	}
	inline void fsearch(const IndexType &pos)
	{
		//std::cout << Qseq_s << ":" << Qseq_e << std::endl;
		if( !(pos >= Qseq_s && pos < Qseq_e) )
		{
			//++f;
			//if(pos >= Qseq_e && Qseq_e != max && Qseq_e != 0){
			if(pos >= Qseq_e && Qseq_e != max && Qseq_e != 0){
				++Qseq_i;
				resetValue();
				if( pos >= Qseq_e ) lower_bound_search(pos);
			}
			else
			{
				lower_bound_search(pos);
			}
			
		}
	}
	
	inline void Qfsearch(const IndexType &pos)
	{
		if( !(pos >= Qseq_s && pos < Qseq_e) )
		{
			if(pos >= Qseq_e && Qseq_e != max && Qseq_e != 0){
				++Qseq_i;
				resetValue();
			}
		}
	}
	
	inline	bool isBegin() const {return (Container.begin() == Qseq_i);}
	inline bool isEnd() const {return (Container.end() == Qseq_i);}
	inline void setMax(IndexType m){max = m;}
	inline void toBegin(){Qseq_i = Container.begin();resetValue();}
	inline DequeIterType begin(){return Container.begin();}
	inline DequeIterType end(){return Container.end();}
	inline DequeIterType find(IndexType v){return Container.find(v);}
	inline bool empty(){return Container.empty();}
	inline IndexType size(){return Container.size();}
	
	
private:
	//std::deque < std::pair < IndexType, ContainerType > > Container;
	std::map < IndexType, ContainerType > Container;
	
	inline void resetValue()
	{
		if(Container.empty())
		{
			Qseq_s = 0;
			Qseq_e = max;
		}
		else
		{	
			Qseq_s = (*Qseq_i).first;
			++Qseq_i;
			if( (Qseq_i) == Container.end() )
			{
				Qseq_e = max;
			}
			else
			{
				Qseq_e = (*(Qseq_i)).first;
			}
			--Qseq_i;
		}
	}
	
};


class ABmpl
{
protected:
	std::string CharExcept;
	
	uint32_t len_seq;
	boost::dynamic_bitset<> sequence;
	
	//uint32_t Qseq_s, Qseq_e;
	//std::deque < std::pair < uint32_t, std::tuple <uint64_t, char, uint32_t> > >::iterator Qseq_i;
	//std::deque < std::pair < uint32_t, std::tuple <uint64_t, char, uint32_t> > > ExceptDeq;
	
	QSeqTable < uint32_t, uint32_t > LowerChar;
	//QSeqTable < uint32_t, uint32_t > LowerPos;
	bool ignore_lower;
	
public:
	//std::vector<uint32_t> LowerPos;
	QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > EscapeChar;
	ABmpl () 	/// default constructor
		: CharExcept (""), sequence (0), len_seq(0), ignore_lower(false), EscapeChar(), LowerChar()
	{}//std::cerr<<"\ndefault constructing"<<std::endl;}
	
	
	std::string MakeSeqString () const
	{
		std::string buffer;
		buffer.reserve (sequence.size()); 
		boost::to_string (sequence, buffer);
		return buffer;
	}
	inline char index_impl (uint32_t p, bool islower = false)	//Implement detail for operator[]
	{
		if (p >= len_seq)
			throw std::out_of_range ("[Error] : ABmpl : index longer than size");
		uint64_t pos {p<<1};
		bool first = sequence [pos];
		bool second = sequence [pos+1];
		if (first && second)
			return islower?'t':'T';
		else if (first && !second)
			return islower?'c':'C';
		else if (!first && second)
			return islower?'g':'G';
		else if (!first && !second)
			return islower?'a':'A';
	}

	void Printer (void)
	{
		
		std::cout<<"EscapeChar info"<<std::endl;
		std::cout << EscapeChar.size() << std::endl;
		std::for_each ( EscapeChar.begin(), 
			EscapeChar.end(),
			[] (const std::pair< uint32_t, std::tuple <uint64_t, char, uint32_t> > Q)
			{
				std::cerr<<Q.first<<'\t'
				<< std::get<0>(Q.second) << '\t'
				<< std::get<1>(Q.second) << '\t'
				<< std::get<2>(Q.second)
				<<std::endl;
			}
		);
		
		/*
		std::cout<<"LowerChar info"<<std::endl;
		std::cout << LowerChar.size() << std::endl;
		std::for_each ( LowerChar.begin(), 
			LowerChar.end(),
			[this] (const std::pair< uint32_t, uint32_t > Q)
			{
				if(Q.second > 1000)
				{
					for(int32_t i = Q.first ;i < Q.first+100; ++i)
					{
						std::cout << index_impl(i);
					}
					std::cout << std::endl;
				}
			}
		);
		*/
		
	}
	inline uint32_t size()
	{
		return len_seq;
	}
	inline uint32_t sizeSeq()
	{
		return sequence.size();
	}
	
	boost::dynamic_bitset<> &bitSequence()
	{
		return sequence;
	}
	inline bool getBitSequence(uint64_t p)
	{
		return sequence[p];
	}
	
};




template <class T = std::string>
class ABSequence
	: public T
{
public:
	ABSequence(T &sequence)
		:T (sequence)
	{
		std::cerr << "Seqence Size : " << this->size() << std::endl;
	}
	ABSequence()
		:T ()
	{
		//std::cerr << this->size() << std::endl;
	}
	std::string& getContent()
	{
		return *this;
	}
	/*
	inline char compare( uint32_t a, uint32_t b, uint32_t len)
	{
		//if(len == 0) len = this->size();
		
		for(auto i(0); i< len; ++i)
		{
			if( a+i >= this->size() )
				return '<';
			else if( b+i >= this->size() )
				return '>';
			//std::cerr << len << std::endl;
			//std::cerr << (a+i) << ":" << (b+i) << ":" << (*this)[b+i] << std::endl;
			if( (*this)[(a+i)] > (*this)[(b+i)] )
				return '>';
			else if( (*this)[(a+i)] < (*this)[(b+i)] )
				return '<';
			//std::cerr << len << std::endl;
		}
		return '=';
	}
	inline bool compare( uint32_t a, uint32_t b)
	{
		for(auto i(0); i< this->size(); ++i)
		{
			if( a+i >= this->size() )
				return true;
			else if( b+i >= this->size() )
				return false;
			
			if( (*this)[(a+i)] > (*this)[(b+i)] )
				return false;
			else if( (*this)[(a+i)] < (*this)[(b+i)] )
				return true;
		}
	}
	*/
};






template<>
class ABSequence<ABmpl>
	: public ABmpl
{
private:
	inline bool FillSeq ( char ChrIn, uint64_t pos, uint32_t length = 1)
	{
		bool first, second;
		switch (ChrIn)
		{
			case 'A':
			case 'a':
				first	= 0;
				second = 0;
				break;

			case 'C':
			case 'c':
				first	= 1;
				second = 0;
				break;

			case 'G':
			case 'g':
				first	= 0;
				second = 1;
				break;

			case 'T':
			case 't':
			case 'U':
			case 'u':
				first	= 1;
				second = 1;
				break;
			case 'N':
				return false;
				break;
		}
		for(uint32_t i = 0; i < length; ++i)
		{
			//if(first) sequence.set( pos + (i << 1) );
			//if(second) sequence.set( pos + (i << 1) +1 );
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
	ABSequence (const std::string& str_seq, bool ignlow = false) /// default constructor
		: ABmpl (), bt(0), ft(0)
	{
		
		uint32_t len_record = 1;
		len_seq = str_seq.size();
		uint64_t pos(0);
		sequence.resize (len_seq << 1);
		uint32_t i (1);
		
		char tmp_repeat_char (str_seq[0]);
		uint32_t tmp_repeat_length (1);
		uint32_t tmp_repeat_limit (100); // > 20 will be record
		
		std::string lower_pattern ("atcgu");
		uint32_t tmp_lower_length (0);
		ignore_lower = ignlow;
		
		EscapeChar.setMax(len_seq);
		LowerChar.setMax(len_seq);
		
		if( lower_pattern.find (str_seq[0]) != std::string::npos )
			++tmp_lower_length;

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
					//if(ignore_lower) tmp_repeat_char = toupper(tmp_repeat_char);
					EscapeChar.insert((i-tmp_repeat_length), std::make_tuple ( pos, tmp_repeat_char, tmp_repeat_length ) );
					//effect lower
					if( !ignore_lower && islower(tmp_repeat_char))	tmp_lower_length -= tmp_repeat_length;
				}//else (tmp_repeat_times <= tmp_repeat_limit)
				
				tmp_repeat_char = str_seq[i];
				tmp_repeat_length = 1;
				
			}// else (str_seq[i] == tmp_repeat_char)
			
			if( ! ignore_lower)
			{
				//if( lower_pattern.find ( str_seq[i] ) != std::string::npos)
				//faster, 3s // 24s -> 21s
				//if( str_seq[i]=='a' || str_seq[i]=='t' || str_seq[i]=='g' || str_seq[i]=='c' || str_seq[i]=='u')
				//faster, 2s // 21s -> 19s
				if(islower(str_seq[i]))
				{
					//lower case
					++tmp_lower_length;
				}
				else
				{
					//uper case
					if(tmp_lower_length != 0)
					{
						// record to table
						LowerChar.insert(( (pos >> 1) - tmp_lower_length), tmp_lower_length );
						tmp_lower_length = 0;
					}
				}
			}
		}// for (uint32_t i = 1; i < len_seq ; ++ i)

		sequence.resize (pos);
	}// ABSequence (const std::string& str_seq)
	
	inline char getSeqChar (const uint32_t &pos)	//Implement detail for operator[]
	{	
		bool islower(false);
		
		if( ! ignore_lower) 
		{
			const uint32_t &start = LowerChar.Qseq_s;
			const uint32_t &end = LowerChar.Qseq_e;
			
			LowerChar.lower_bound_search(pos);
			
			if(pos >= start && pos < (start + (*LowerChar.Qseq_i).second ) )
				islower = true;
		}	
		return index_impl(pos, islower);
	}
	inline char mt_at (const uint32_t &pos)
	{
		QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > EscapeCharTmp(EscapeChar);
		if(EscapeCharTmp.empty())	return getSeqChar(pos);
		
		uint32_t &start = EscapeCharTmp.Qseq_s;
		uint32_t &end = EscapeCharTmp.Qseq_e;
		
		EscapeCharTmp.lower_bound_search(pos);
		if(EscapeCharTmp.isBegin() && pos < start)	return getSeqChar(pos);
		
		uint32_t repeat_end_pos(	(*EscapeCharTmp.Qseq_i).first + std::get<2>( (*EscapeCharTmp.Qseq_i).second ) - 1	);
		if(pos <= repeat_end_pos )
		{
			return std::get<1>( (*EscapeCharTmp.Qseq_i).second );
		}
		else
		{
			return getSeqChar(	pos - repeat_end_pos + ((std::get<0>( (*EscapeCharTmp.Qseq_i).second )) >> 1) -1	);
		}
	}
	inline char operator[] (const uint32_t &pos)
	{
		if(EscapeChar.empty())	return getSeqChar(pos);
		
		uint32_t &start = EscapeChar.Qseq_s;
		uint32_t &end = EscapeChar.Qseq_e;
		
		EscapeChar.lower_bound_search(pos);
		if(EscapeChar.isBegin() && pos < start)	return getSeqChar(pos);
		
		uint32_t repeat_end_pos(	(*EscapeChar.Qseq_i).first + std::get<2>( (*EscapeChar.Qseq_i).second ) - 1	);
		if(pos <= repeat_end_pos )
		{
			return std::get<1>( (*EscapeChar.Qseq_i).second );
		}
		else
		{
			return getSeqChar(	pos - repeat_end_pos + ((std::get<0>( (*EscapeChar.Qseq_i).second )) >> 1) -1	);
		}
	}
	inline void toString(std::string &buffer)
	{
		buffer.clear();
		buffer.reserve(len_seq);
		
		uint32_t pos(0);
		uint64_t p(0);
		EscapeChar.toBegin();
		LowerChar.toBegin();
		
		
		bool islower(false);
		
		while(pos < len_seq)
		{
			//Escape
			if(pos == EscapeChar.Qseq_s)
			{
				buffer += std::string( std::get<2>( (*EscapeChar.Qseq_i).second ), std::get<1>( (*EscapeChar.Qseq_i).second ) );
				pos += std::get<2>( (*EscapeChar.Qseq_i).second );
				EscapeChar.next();
			}
			else
			{
				if( p == LowerChar.Qseq_s + (*LowerChar.Qseq_i).second )	
					LowerChar.next();
				//if(p >= LowerChar.Qseq_s )
				if(p >= LowerChar.Qseq_s && p < LowerChar.Qseq_s + (*LowerChar.Qseq_i).second )
					islower = true;
				else
					islower = false;

				buffer.push_back( index_impl( p, islower) );
				++pos;
				++p;
			}
		}
	}
	
	
	/* comapre const length, return '>'(a>b), '<'(a<b), '='(==) */
	inline char compare(uint32_t a, uint32_t b, uint32_t len,
		QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > &quick_table_A,
		QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > &quick_table_B
	)
	{
		// two EscapeChar iterator
		//static QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > quick_table_A(EscapeChar);
		//static QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > quick_table_B(EscapeChar);
		
		// repeat_end_pos is the last position of escape chars
		uint32_t repeat_end_pos_A(0), repeat_end_pos_B(0), min(0), ori_a(a), ori_b(b);
		char returnValue('=');
		// the value of the position
		char value_A('\0'), value_B('\0');
		
		//move iterator to right position
		quick_table_A.lower_bound_search(a);
		quick_table_B.lower_bound_search(b);
		

		while( a < (ori_a+len) && b < (ori_b+len) && a < len_seq && b < len_seq)
		{
			//std::cout << a << ":" << ":"<< ":" << b << std::endl;	
			// If A_position is less then the last position of escape chars,
			// it can escape to read a char until a position is bigger then the pos of repeat end.
			// And A special method is '=', because it needs to read a char to initial the system.
			// To begin, repeat_end_pos_A is 0, so any A_position is bigger then 0, or equal to 0.
			if( a >= repeat_end_pos_A)
			{
				// if A_position >= Qseq_e, move to next iterator.
				if( a >= quick_table_A.Qseq_e)
					quick_table_A.next();
					
				// the least Qseq_s of escape table might is not 0, so it need to exam the A_position
				// smaller than quick_table_A.Qseq_s. A_position maybe smaller than Qseq start.
				if(a < quick_table_A.Qseq_s || quick_table_A.size() ==0)
				{	
					value_A = index_impl(a);
					
				}
				else
				// In normal escaper table block, between Qseq start and Qseq end.
				// In detail, you can see function operator[]
				{
					repeat_end_pos_A = (*quick_table_A.Qseq_i).first + std::get<2>( (*quick_table_A.Qseq_i).second ) - 1;
					
					// In escaper chars block.
					if(a <= repeat_end_pos_A)
						value_A = toupper(std::get<1>( (*quick_table_A.Qseq_i).second ));
					else
						value_A = index_impl(	a - repeat_end_pos_A + ((std::get<0>( (*quick_table_A.Qseq_i).second )) >> 1) -1	);					
				}
			}
			// The same with A_position
			if( b >= repeat_end_pos_B )
			{
				if( b >= quick_table_B.Qseq_e)
					quick_table_B.next();
				
				if(b < quick_table_B.Qseq_s || quick_table_B.size() ==0)
				{
					value_B = index_impl(b);
				}
				else
				{
					repeat_end_pos_B = (*quick_table_B.Qseq_i).first + std::get<2>( (*quick_table_B.Qseq_i).second ) - 1;

					if(b <= repeat_end_pos_B)
						value_B = toupper(std::get<1>( (*quick_table_B.Qseq_i).second ));
					else
						value_B = index_impl(	b - repeat_end_pos_B + ((std::get<0>( (*quick_table_B.Qseq_i).second )) >> 1) -1	);
				}
			}
			// compare finish
			//std::cout << value_A << ":" << value_B << "\t" << a << ":" << b << std::endl;
			if( value_B > value_A )
			{
				returnValue = '<';
				break;
			}
			else if ( value_A > value_B)
			{
				returnValue = '>';
				break;
			}
			// From now, it means two value is euqal. value_A == value_B
			// If two value is euqal and both in escape chars block, Jump!
			if(a < repeat_end_pos_A && b < repeat_end_pos_B)
			{
				min = std::min( repeat_end_pos_A - a, repeat_end_pos_B - b ) ;
				a += min;
				b += min;
			}
			++a;
			++b;
		}
	
		if(returnValue != '=') return returnValue;
		
		if( ( a >= (ori_a+len) || b >=(ori_b+len) ) && (ori_a+len) != len_seq && (ori_b+len) != len_seq )
		{
			returnValue = '=';
		}
		else
		{
			if(a > b)
				returnValue = '<';
			else if(a < b)
				returnValue = '>';
		}
		
		//std::cout << a << ":" << ori_a << std::endl;
		
		
		return returnValue;
	}
	
	
	
	/* input 2 sequence location to compare string value */
	inline bool compare(uint32_t a, uint32_t b,
		QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > &quick_table_A,
		QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > &quick_table_B
	)
	{
		// two EscapeChar iterator
		//static QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > quick_table_A(EscapeChar);
		//static QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > quick_table_B(EscapeChar);
		
		// repeat_end_pos is the last position of escape chars
		uint32_t repeat_end_pos_A(0), repeat_end_pos_B(0), min(0);
		bool returnValue(false);
		// the value of the position
		char value_A('\0'), value_B('\0');
		
		//move iterator to right position
		quick_table_A.lower_bound_search(a);
		quick_table_B.lower_bound_search(b);

		while(a < len_seq && b < len_seq)
		{
			// If A_position is less then the last position of escape chars,
			// it can escape to read a char until a position is bigger then the pos of repeat end.
			// And A special method is '=', because it needs to read a char to initial the system.
			// To begin, repeat_end_pos_A is 0, so any A_position is bigger then 0, or equal to 0.
			if( a >= repeat_end_pos_A)
			{
				// if A_position >= Qseq_e, move to next iterator.
				if( a >= quick_table_A.Qseq_e)
					quick_table_A.next();
					
				// the least Qseq_s of escape table might is not 0, so it need to exam the A_position
				// smaller than quick_table_A.Qseq_s. A_position maybe smaller than Qseq start.
				if(a < quick_table_A.Qseq_s )
				{	
					value_A = index_impl(a);
				}
				else
				// In normal escaper table block, between Qseq start and Qseq end.
				// In detail, you can see function operator[]
				{
					repeat_end_pos_A = (*quick_table_A.Qseq_i).first + std::get<2>( (*quick_table_A.Qseq_i).second ) - 1;
					// In escaper chars block.
					if(a <= repeat_end_pos_A)
						value_A = toupper(std::get<1>( (*quick_table_A.Qseq_i).second ));
					else
						value_A = index_impl(	a - repeat_end_pos_A + ((std::get<0>( (*quick_table_A.Qseq_i).second )) >> 1) -1	);					
				}
			}
			// The same with A_position
			if( b >= repeat_end_pos_B )
			{
				if( b >= quick_table_B.Qseq_e)
					quick_table_B.next();
				
				if(b < quick_table_B.Qseq_s )
				{
					value_B = index_impl(b);
				}
				else
				{
					repeat_end_pos_B = (*quick_table_B.Qseq_i).first + std::get<2>( (*quick_table_B.Qseq_i).second ) - 1;
					if(b <= repeat_end_pos_B)
						value_B = toupper(std::get<1>( (*quick_table_B.Qseq_i).second ));
					else
						value_B = index_impl(	b - repeat_end_pos_B + ((std::get<0>( (*quick_table_B.Qseq_i).second )) >> 1) -1	);
				}
			}
			// compare finish
			if( value_B > value_A )
			{
				returnValue = true;
				break;
			}
			else if ( value_A > value_B)
			{
				returnValue = false;
				break;
			}
			// From now, it means two value is euqal. value_A == value_B
			// If two value is euqal and both in escape chars block, Jump!
			if(a < repeat_end_pos_A && b < repeat_end_pos_B)
			{
				min = std::min( repeat_end_pos_A - a, repeat_end_pos_B - b ) ;
				a += min;
				b += min;
			}
			++a;
			++b;
		}
		
		if(a == len_seq || b == len_seq)
		{
			if( b > a )
				returnValue = false;
			else if ( b < a)
				returnValue = true;
		}
		return returnValue;
	}
	
	// the same as before, but using static QseqTables
	inline bool compare(uint32_t a, uint32_t b)
	{
		static QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > quick_table_A(EscapeChar);
		static QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > quick_table_B(EscapeChar);
		return compare(a, b, quick_table_A, quick_table_B);
	}
	
	
	// for run time jump and fix length
	inline char compare(uint32_t a, uint32_t b, uint32_t len,
		QSeqTable < uint32_t, std::map<uint32_t, uint32_t> > &RunTimeJumpTable,
		QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > &quick_table_A,
		QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > &quick_table_B
	)
	{
		// two EscapeChar iterator
		//static QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > quick_table_A(EscapeChar);
		//static QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > quick_table_B(EscapeChar);
		
		// repeat_end_pos is the last position of escape chars
		uint32_t repeat_end_pos_A(0), repeat_end_pos_B(0), min(0), jumpCount(0), jumpCountWithoutN(0), ori_a(a), ori_b(b);
		
		// the value of the position
		char value_A('\0'), value_B('\0'), returnValue('=');
		
		//move iterator to right position
		quick_table_A.lower_bound_search(a);
		quick_table_B.lower_bound_search(b);
		
		std::map<uint32_t, uint32_t>::iterator sub_jumpIt(0);
		bool flag(false), jump(false), find(false);
		uint32_t jumpInfo_start(0), jumpInfo_distance(0), jumpInfo_length(0),min_ab(0),distance_ab(0);
		
		min_ab = std::min(a,b);
		distance_ab = abs(a-b);
		//distance_ab = b-a;
		
		//search for jump table to jump
		//RunTimeJumpTable.lower_bound(std::min(a,b));
		//pos, distance, length

		//jumpInfo_distance = ((*RunTimeJumpTable.Qseq_i).second).first;
		//auto tmpIt = (RunTimeJumpTable.Qseq_i);
		//if(a==214530711 && b==155157145)	jump=true;
		
		
		if(RunTimeJumpTable.size() != 0)
		{
			RunTimeJumpTable.jump_lower_bound(min_ab);
			if(min_ab >= RunTimeJumpTable.Qseq_s )
			{
				jumpInfo_start = RunTimeJumpTable.Qseq_s;
				sub_jumpIt = ((*RunTimeJumpTable.Qseq_i).second).find(distance_ab);
				if(sub_jumpIt != ((*RunTimeJumpTable.Qseq_i).second).end() )
				{
					find = true;
					jumpInfo_distance = distance_ab;
					jumpInfo_length = sub_jumpIt->second;
					if( min_ab < ( jumpInfo_start + jumpInfo_length -1 ) )
					{
						a += jumpInfo_length - (min_ab - jumpInfo_start);
						b += jumpInfo_length - (min_ab - jumpInfo_start);
					}
				}
			}
		}
		
		
		int32_t run_min_ab(0), run_jump_distance(0);
		run_min_ab = std::min(a,b);
		
		while( a < (ori_a+len) && b < (ori_b+len) && a < len_seq && b < len_seq)
		//while(a < len_seq && b < len_seq)
		{
			// Run time jump
			if(RunTimeJumpTable.size() != 0)
			{
				RunTimeJumpTable.lower_bound_search(run_min_ab);
				
				if(run_min_ab >= RunTimeJumpTable.Qseq_s )
				{
					jumpInfo_start = RunTimeJumpTable.Qseq_s;
					sub_jumpIt = ((*RunTimeJumpTable.Qseq_i).second).find(distance_ab);
					if(sub_jumpIt != ((*RunTimeJumpTable.Qseq_i).second).end() )
					{
						//find = true;
						jumpInfo_distance = distance_ab;
						jumpInfo_length = sub_jumpIt->second;
						if( run_min_ab < ( jumpInfo_start + jumpInfo_length -1 ) )
						{
							//std::cout << "jump:" << jumpInfo_length - (a - jumpInfo_start) << "\n";
							run_jump_distance = jumpInfo_length - (run_min_ab - jumpInfo_start)-1;
							a += run_jump_distance;
							b += run_jump_distance;
							run_min_ab += run_jump_distance;
							jumpCount += run_jump_distance;
							jumpCountWithoutN += run_jump_distance;
						}
					}
				}
			}
			
		
			// If A_position is less then the last position of escape chars,
			// it can escape to read a char until a position is bigger then the pos of repeat end.
			// And A special method is '=', because it needs to read a char to initial the system.
			// To begin, repeat_end_pos_A is 0, so any A_position is bigger then 0, or equal to 0.
			if( a >= repeat_end_pos_A)
			{
				// if A_position >= Qseq_e, move to next iterator.
				quick_table_A.lower_bound_search(a);
				//if( a >= quick_table_A.Qseq_e)
				//	quick_table_A.next();
					
				// the least Qseq_s of escape table might is not 0, so it need to exam the A_position
				// smaller than quick_table_A.Qseq_s. A_position maybe smaller than Qseq start.
				if(a < quick_table_A.Qseq_s || quick_table_A.size()==0)
				{	
					value_A = index_impl(a);
				}
				else
				// In normal escaper table block, between Qseq start and Qseq end.
				// In detail, you can see function operator[]
				{
					repeat_end_pos_A = (*quick_table_A.Qseq_i).first + std::get<2>( (*quick_table_A.Qseq_i).second ) - 1;
					// In escaper chars block.
					if(a <= repeat_end_pos_A)
						value_A = toupper(std::get<1>( (*quick_table_A.Qseq_i).second ));
					else
						value_A = index_impl(	a - repeat_end_pos_A + ((std::get<0>( (*quick_table_A.Qseq_i).second )) >> 1) -1	);					
				}
			}
			// The same with A_position
			if( b >= repeat_end_pos_B )
			{
				quick_table_B.lower_bound_search(b);
				//if( b >= quick_table_B.Qseq_e)
				//	quick_table_B.next();
				
				if(b < quick_table_B.Qseq_s || quick_table_B.size()==0)
				{
					value_B = index_impl(b);
				}
				else
				{
					repeat_end_pos_B = (*quick_table_B.Qseq_i).first + std::get<2>( (*quick_table_B.Qseq_i).second ) - 1;
					if(b <= repeat_end_pos_B)
						value_B = toupper(std::get<1>( (*quick_table_B.Qseq_i).second ));
					else
						value_B = index_impl(	b - repeat_end_pos_B + ((std::get<0>( (*quick_table_B.Qseq_i).second )) >> 1) -1	);
				}
			}
			// compare finish
			if( value_B > value_A )
			{
				returnValue = '<';
				break;
			}
			else if ( value_A > value_B)
			{
				returnValue = '>';
				break;
			}
			
			// From now, it means two value is euqal. value_A == value_B
			// If two value is euqal and both in escape chars block, Jump!
			if(a < repeat_end_pos_A && b < repeat_end_pos_B)
			{
				min = std::min( repeat_end_pos_A - a, repeat_end_pos_B - b ) ;
				a += min;
				b += min;
				jumpCount	+= min;
				run_min_ab += min;
			}
			++run_min_ab;
			++a;
			++b;
			++jumpCount;
			++jumpCountWithoutN;
			if ( !flag && !find && jumpCountWithoutN >= 5000 )
			{
				flag = true;
			}
				
		}
		
		
		if(returnValue != '=') return returnValue;
		
		if( ( a >= (ori_a+len) || b >=(ori_b+len) ) && (ori_a+len) != len_seq && (ori_b+len) != len_seq )
		{
			returnValue = '=';
		}
		else
		{
			if(a > b)
				returnValue = '<';
			else if(a < b)
				returnValue = '>';
		}
		
		//pos, distance, length
		
		if(flag)
		{
			//if(RunTimeJumpTable.Qseq_s != min_ab	|| (RunTimeJumpTable.Qseq_s==0 && min_ab ==0))
			//if(RunTimeJumpTable.find(min_ab) == RunTimeJumpTable.end() )
			if(RunTimeJumpTable.Qseq_s != min_ab)
			{
				
				//new
				//std::cerr << min_ab << ":" << jumpInfo_start << std::endl;
				RunTimeJumpTable.insert( min_ab, {{distance_ab, jumpCount}} );
				//if(ori_a == 156337867)	std::cout << "Y5:" << RunTimeJumpTable.Qseq_s << ":" << min_ab << std::endl;
			}
			else
			{
				
				//std::cerr << min_ab << ":" << jumpInfo_start << std::endl;
				//have been insert
				((*RunTimeJumpTable.Qseq_i).second).insert( {distance_ab, jumpCount} );
				//if(ori_a == 156337867)	std::cout << "Y6" << std::endl;
			}
			//if(ori_a == 156337867)	
			//std::cout << (*this)[min_ab] << "->" << (*this)[a] << ":" << ori_a << ":" << ori_b << ":" <<	min_ab << ":" << distance_ab << ":" << jumpCount << std::endl;
			//RunTimeJumpTable.insert( std::min(ori_a,ori_b), std::make_pair( abs(ori_b - ori_a), jumpCount+1 ));
			//std::cout << "insert:"<< jumpCount << "\n";
		}
		
		return returnValue;
	}
	
	
	
	// for run time jump
	inline bool compare(uint32_t a, uint32_t b,
		//QSeqTable < uint32_t, std::pair<uint32_t, uint32_t> > &RunTimeJumpTable,
		QSeqTable < uint32_t, std::map<uint32_t, uint32_t> > &RunTimeJumpTable,
		QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > &quick_table_A,
		QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > &quick_table_B
	)
	{
		// two EscapeChar iterator
		//static QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > quick_table_A(EscapeChar);
		//static QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > quick_table_B(EscapeChar);
		
		// repeat_end_pos is the last position of escape chars
		uint32_t repeat_end_pos_A(0), repeat_end_pos_B(0), min(0), jumpCount(0), jumpCountWithoutN(0), ori_a(a), ori_b(b);
		
		// the value of the position
		char value_A('\0'), value_B('\0');
		
		//move iterator to right position
		quick_table_A.lower_bound_search(a);
		quick_table_B.lower_bound_search(b);
		
		std::map<uint32_t, uint32_t>::iterator sub_jumpIt(0);
		bool returnValue(false), flag(false), jump(false), find(false);
		uint32_t jumpInfo_start(0), jumpInfo_distance(0), jumpInfo_length(0),min_ab(0),distance_ab(0);
		
		min_ab = std::min(a,b);
		distance_ab = abs(a-b);
		//distance_ab = b-a;
		
		//search for jump table to jump
		//RunTimeJumpTable.lower_bound(std::min(a,b));
		//pos, distance, length

		//jumpInfo_distance = ((*RunTimeJumpTable.Qseq_i).second).first;
		//auto tmpIt = (RunTimeJumpTable.Qseq_i);
		//if(a==214530711 && b==155157145)	jump=true;
		
		if(RunTimeJumpTable.size() != 0)
		{
			RunTimeJumpTable.jump_lower_bound(min_ab);
			if(min_ab >= RunTimeJumpTable.Qseq_s )
			{
				jumpInfo_start = RunTimeJumpTable.Qseq_s;
				sub_jumpIt = ((*RunTimeJumpTable.Qseq_i).second).find(distance_ab);
				if(sub_jumpIt != ((*RunTimeJumpTable.Qseq_i).second).end() )
				{
					find = true;
					jumpInfo_distance = distance_ab;
					jumpInfo_length = sub_jumpIt->second;
					if( min_ab < ( jumpInfo_start + jumpInfo_length -1 ) )
					{
						a += jumpInfo_length - (min_ab - jumpInfo_start);
						b += jumpInfo_length - (min_ab - jumpInfo_start);
					}
				}
			}
		}
		
		int32_t run_min_ab(0), run_jump_distance(0);
		run_min_ab = std::min(a,b);
		while(a < len_seq && b < len_seq)
		{
			// Run time jump
			if(RunTimeJumpTable.size() != 0)
			{
				RunTimeJumpTable.lower_bound_search(run_min_ab);
				
				if(run_min_ab >= RunTimeJumpTable.Qseq_s )
				{
					jumpInfo_start = RunTimeJumpTable.Qseq_s;
					sub_jumpIt = ((*RunTimeJumpTable.Qseq_i).second).find(distance_ab);
					if(sub_jumpIt != ((*RunTimeJumpTable.Qseq_i).second).end() )
					{
						//find = true;
						jumpInfo_distance = distance_ab;
						jumpInfo_length = sub_jumpIt->second;
						if( run_min_ab < ( jumpInfo_start + jumpInfo_length -1 ) )
						{
							//std::cout << "jump:" << jumpInfo_length - (a - jumpInfo_start) << "\n";
							run_jump_distance = jumpInfo_length - (run_min_ab - jumpInfo_start)-1;
							a += run_jump_distance;
							b += run_jump_distance;
							run_min_ab += run_jump_distance;
							jumpCount += run_jump_distance;
							jumpCountWithoutN += run_jump_distance;
						}
					}
				}
			}
			
		
			// If A_position is less then the last position of escape chars,
			// it can escape to read a char until a position is bigger then the pos of repeat end.
			// And A special method is '=', because it needs to read a char to initial the system.
			// To begin, repeat_end_pos_A is 0, so any A_position is bigger then 0, or equal to 0.
			if( a >= repeat_end_pos_A)
			{
				// if A_position >= Qseq_e, move to next iterator.
				quick_table_A.lower_bound_search(a);
				//if( a >= quick_table_A.Qseq_e)
				//	quick_table_A.next();
					
				// the least Qseq_s of escape table might is not 0, so it need to exam the A_position
				// smaller than quick_table_A.Qseq_s. A_position maybe smaller than Qseq start.
				if(a < quick_table_A.Qseq_s || quick_table_A.size()==0)
				{	
					value_A = index_impl(a);
				}
				else
				// In normal escaper table block, between Qseq start and Qseq end.
				// In detail, you can see function operator[]
				{
					repeat_end_pos_A = (*quick_table_A.Qseq_i).first + std::get<2>( (*quick_table_A.Qseq_i).second ) - 1;
					// In escaper chars block.
					if(a <= repeat_end_pos_A)
						value_A = toupper(std::get<1>( (*quick_table_A.Qseq_i).second ));
					else
						value_A = index_impl(	a - repeat_end_pos_A + ((std::get<0>( (*quick_table_A.Qseq_i).second )) >> 1) -1	);					
				}
			}
			// The same with A_position
			if( b >= repeat_end_pos_B )
			{
				quick_table_B.lower_bound_search(b);
				//if( b >= quick_table_B.Qseq_e)
				//	quick_table_B.next();
				
				if(b < quick_table_B.Qseq_s || quick_table_B.size()==0)
				{
					value_B = index_impl(b);
				}
				else
				{
					repeat_end_pos_B = (*quick_table_B.Qseq_i).first + std::get<2>( (*quick_table_B.Qseq_i).second ) - 1;
					if(b <= repeat_end_pos_B)
						value_B = toupper(std::get<1>( (*quick_table_B.Qseq_i).second ));
					else
						value_B = index_impl(	b - repeat_end_pos_B + ((std::get<0>( (*quick_table_B.Qseq_i).second )) >> 1) -1	);
				}
			}
			// compare finish
			if( value_B > value_A )
			{
				returnValue = true;
				break;
			}
			else if ( value_A > value_B)
			{
				returnValue = false;
				break;
			}
			
			// From now, it means two value is euqal. value_A == value_B
			// If two value is euqal and both in escape chars block, Jump!
			if(a < repeat_end_pos_A && b < repeat_end_pos_B)
			{
				min = std::min( repeat_end_pos_A - a, repeat_end_pos_B - b ) ;
				a += min;
				b += min;
				jumpCount	+= min;
				run_min_ab += min;
			}
			++run_min_ab;
			++a;
			++b;
			++jumpCount;
			++jumpCountWithoutN;
			if ( !flag && !find && jumpCountWithoutN >= 1000 )
			{
				flag = true;
			}
				
		}
		
		if(a == len_seq || b == len_seq)
		{
			if( b > a )
				returnValue = false;
			else if ( b < a)
				returnValue = true;
		}
		
		//pos, distance, length
		
		if(flag)
		{
			//if(RunTimeJumpTable.Qseq_s != min_ab	|| (RunTimeJumpTable.Qseq_s==0 && min_ab ==0))
			//if(RunTimeJumpTable.find(min_ab) == RunTimeJumpTable.end() )
			if(RunTimeJumpTable.Qseq_s != min_ab)
			{
				
				//new
				//std::cerr << min_ab << ":" << jumpInfo_start << std::endl;
				RunTimeJumpTable.insert( min_ab, {{distance_ab, jumpCount}} );
				//if(ori_a == 156337867)	std::cout << "Y5:" << RunTimeJumpTable.Qseq_s << ":" << min_ab << std::endl;
			}
			else
			{
				
				//std::cerr << min_ab << ":" << jumpInfo_start << std::endl;
				//have been insert
				((*RunTimeJumpTable.Qseq_i).second).insert( {distance_ab, jumpCount} );
				//if(ori_a == 156337867)	std::cout << "Y6" << std::endl;
			}
			//if(ori_a == 156337867)	
			//std::cout << (*this)[min_ab] << "->" << (*this)[a] << ":" << ori_a << ":" << ori_b << ":" <<	min_ab << ":" << distance_ab << ":" << jumpCount << std::endl;
			//RunTimeJumpTable.insert( std::min(ori_a,ori_b), std::make_pair( abs(ori_b - ori_a), jumpCount+1 ));
			//std::cout << "insert:"<< jumpCount << "\n";
		}
		if(jump) std::cout << returnValue << std::endl;
		return returnValue;
	}
	// A stupid method (very slow) to compare string, but it is the true answer.
	
	
	
	//for limit length
	inline bool compare(uint32_t a, std::string &comapre_string, bool sameValue,
		QSeqTable < uint32_t, std::tuple<uint64_t, char, uint32_t> > &quick_table_A
	)
	{
		uint32_t b(0),repeat_end_pos_A(0);

		char value_B('\0'), value_A('\0');
		quick_table_A.lower_bound_search(a);

		while(a < len_seq)
		{
			if( a >= repeat_end_pos_A)
			{
				if( a >= quick_table_A.Qseq_e)
				{
					quick_table_A.next();
				}
				if(a < quick_table_A.Qseq_s )
				{	
					value_A = index_impl(a);
				}
				else
				{
					repeat_end_pos_A = (*quick_table_A.Qseq_i).first + std::get<2>( (*quick_table_A.Qseq_i).second ) - 1;
					// In escaper chars block.
					if(a <= repeat_end_pos_A)
					{
						value_A = toupper(std::get<1>( (*quick_table_A.Qseq_i).second ));
					}
					else
					{
						value_A = index_impl(	a - repeat_end_pos_A + ((std::get<0>( (*quick_table_A.Qseq_i).second )) >> 1) -1	);
					}		
				}
			}
				
			value_B = comapre_string[b];
			// compare finish
			if( value_B > value_A )
			{
				return true;
			}
				
			else if ( value_A > value_B)
			{
				return false;
			}
				
			if(b == comapre_string.length()-1)
			{
				return sameValue;
			}
				
			++a;
			++b;
		}
		//the length of a is smaller than b
		//if(b <= comapre_string.length())
		//	return true;
		return true;
	}
	
	
	inline bool compareRT2(uint32_t a, uint32_t b)
	{
		char A,B;
		while(a != len_seq && b != len_seq)
		{
			A = toupper((*this)[a]);
			B = toupper((*this)[b]);
			if(B > A )
			{
				return true;
			}
			else if ( A > B)
			{
				return false;
			}
			++a;
			++b;
		}
		if(b > a)	return false;
		return true;
	}
	
	void test()
	{
		std::cout << "EscapeChar: b:" << EscapeChar.b << " f: " << EscapeChar.f << std::endl; 
		std::cout << "LowerChar: b:" << LowerChar.b << " f: " << LowerChar.f << std::endl; 
	}
	
};









#endif
