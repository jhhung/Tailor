#ifndef ABWT_TABLE_HPP_
#define ABWT_TABLE_HPP_
#include <algorithm>
#include <bitset>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <locale>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <utility>
//#include <unordered_map>

#include <boost/ref.hpp>
#include <ctime>
#include <cmath>
#include <memory>
#include "compression/abit.hpp"
#include "compression/jbit.hpp"
//#include "thread_pool.hpp"
//#include "difference_cover.hpp"
#include "constant_def.hpp"
#include "mkq_sort.hpp"
#include "difference_cover.hpp"
#include "split_sort.hpp"


#include "boost/serialization/utility.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/unordered_map.hpp"
//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/archive/binary_oarchive.hpp>

//#define INTTYPE uint64_t

#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/device/file.hpp"
#include "boost/iostreams/filter/zlib.hpp"
#include "boost/serialization/map.hpp"


class ABWT_table
{
public:


	INTTYPE interval, char_size, powerV, first_location;
	std::string bwt;
	std::vector<INTTYPE> c_function;
	std::vector<INTTYPE> c_functions;
	//std::vector<INTTYPE> c_function_sp;
	
	std::vector< std::vector<INTTYPE> > occ_function;
	std::vector< std::pair<INTTYPE, INTTYPE> > location_table;
	
	std::vector<int> occ_char;
	
	//boost::unordered_map<char, INTTYPE> mtable;
	std::vector<INTTYPE> mtable;
	std::vector< std::vector< std::vector<INTTYPE> > > occ_jbwt;
	
	std::string jbwt_idx_char;
	
	std::vector<uint8_t> fbwt;
	std::vector<uint8_t> jbwt_seq;
	
	std::shared_ptr<JBit> jbwt;
	
	// FIXME: this is just temperarily storing the real size of the genome
	uint64_t _realSize = 0;
	// starting site -> chromosome
	std::map <uint64_t, std::string> chr_start_pos {};
	// chr -> chr size
	std::map <std::string, uint64_t> chr_length {};
	// unambiguous segment sequence starting position of each chr
	std::map <uint64_t, uint64_t > chr_umbiguous_starting_length {};

	double get_c_time;
	double get_occ_time;
	
	//
	ABWT_table(INTTYPE iv)
		: interval(iv)
		, char_size(256) 
		, c_function(char_size, 0)//, occ_function(bwt_size/interval+1, std::vector<INTTYPE>(char_size, 0))//, location_table(bwt_size/interval+1,0)
		, c_functions()
		, get_c_time()
		, get_occ_time()
		//, mtable {{'A',0},{'C',1},{'G',2},{'T',3}}
		//, occ_jbwt(256, std::vector<INTTYPE>(256,0))
		, mtable(256)
		, occ_jbwt(256, std::vector< std::vector<INTTYPE> >(5,std::vector<INTTYPE>(256,0) ) )
		, occ_char({'A','C','G','T'})
		
	{
		interval_set(interval);
		mtable['A']=0;
		mtable['C']=1;
		mtable['G']=2;
		mtable['T']=3;
	}
	ABWT_table()
		: interval(0)
		, char_size(256) 
		, c_function(char_size, 0)//, occ_function(bwt_size/interval+1, std::vector<INTTYPE>(char_size, 0))//, location_table(bwt_size/interval+1,0)
		, c_functions()
		, get_c_time()
		, get_occ_time()
		//, mtable {{'A',0},{'C',1},{'G',2},{'T',3}}
		//, occ_jbwt(256, std::vector<INTTYPE>(256,0))
		, mtable(256)
		, occ_jbwt(256, std::vector< std::vector<INTTYPE> >(5,std::vector<INTTYPE>(256,0) ) )
		, occ_char({'A','C','G','T'})
	{
		mtable['A']=0;
		mtable['C']=1;
		mtable['G']=2;
		mtable['T']=3;
	}
	


	void readChrStartPos (const std::string& chrStartPosFile) {
		std::ifstream in {chrStartPosFile};
		std::string line {}, chr {};
		uint64_t startPos {0};
		while (getline (in, line)) {
			std::stringstream ss {line};
			ss >> chr >> startPos;
			chr_start_pos.insert (std::make_pair (startPos, chr));
		}
	}

	void readChrLen (const std::string& chrLenFile) {
		std::ifstream in {chrLenFile};
		std::string line {}, chr {};
		uint64_t length {0};
		while (getline (in, line)) {
			std::stringstream ss {line};
			ss >> chr >> length;
			// FIXME: replace _realSize...
			_realSize += length;
			if (chr_length.find (chr) == chr_length.end())
				chr_length.insert (std::make_pair (chr, length));
			else {
				std::cerr << "Error: duplicated chromosome name" << std::endl;
				exit (1);
			}
		}
	}
	void readNPosLen (const std::string& fileName) {
		boost::iostreams::filtering_istream fis;
		fis.push (boost::iostreams::zlib_decompressor());
		fis.push (boost::iostreams::file_source (fileName));
		boost::archive::binary_iarchive iar (fis);
		iar >> chr_umbiguous_starting_length;
//		std::ifstream in {fileName};
//		uint64_t a, b;
//		while (in.good ()) {
//			in >> a >> b;
//			chr_umbiguous_starting_length.insert (std::make_pair (a,b));
//		}
	}
	void interval_set(INTTYPE iv)
	{
		interval = iv;
		
		if(interval == 2)
			powerV = 1;
		else if(interval == 4)
			powerV = 2;
		else if(interval == 8)
			powerV = 3;
		else if(interval == 16)
			powerV = 4;
		else if(interval == 32)
			powerV = 5;
		else if(interval == 64)
			powerV = 6;
		else if(interval == 128)
			powerV = 7;
		else if(interval == 256)
			powerV = 8;
		else if(interval == 512)
			powerV = 9;
		else if(interval == 1024)
			powerV = 10;
		
	}
	
	void using_jbwt()
	{
		//std::shared_ptr<JBit> kk(new JBit(bwt));
		//jbwt = std::shared_ptr<JBit>( new JBit(bwt) );
		
		
		std::string all_char("ACGT");
		for( char c : all_char)
		{
			for( std::pair< const uint8_t, std::array<char, 4> > &byte : jbwt->byte2word )
			{
				std::vector<INTTYPE> char_number(5,0);
				//for(char c_in_array : byte.second)
				for(int i(0); i<4; ++i)
				{
					char c_in_array(byte.second[i]);	
					if(c_in_array == c)
					{
						for(int j(0); j<5; ++j)
						{
							if( i < j)
								++char_number[j];
						}
					}
				}
				for(int j(0); j<5; ++j)
				{
					occ_jbwt[c][j][byte.first] = char_number[j];
				}
				//std::cerr << "{{" << (int)c << "," << (int)byte.first << "}, " << char_number << "}" << std::endl;
				//occ_jbwt.insert(	std::make_pair(c, byte.first), char_number	);
			}
		}
		
		jbwt_idx_char = "";
		for(int i(0); i<256; i++)
		{
			char w1 = all_char[ (i&255) >> 6 ];
			char w2 = all_char[ (i&63) >> 4 ];
			char w3 = all_char[ (i&15) >> 2 ];
			char w4 = all_char[ (i&3) >> 0 ];
			jbwt_idx_char += w1;
			jbwt_idx_char += w2;
			jbwt_idx_char += w3;
			jbwt_idx_char += w4;
		}
		
		std::cerr << "jbwt ok" << std::endl;
		
		//
	}
	
	template< class SEQTYPE>
	void saveSEQ(std::string filename, SEQTYPE &seq)
	{
		std::ofstream fp(filename, std::ios::binary);
		boost::archive::binary_oarchive archive_fp( fp );
		archive_fp & seq.getContent();
		fp.close();
	}
	template< class SEQTYPE>
	void readSEQ(std::string filename, SEQTYPE &seq)
	{
		std::ifstream fp(filename, std::ios::binary);
		boost::archive::binary_iarchive archive_fp( fp );
		archive_fp & seq.getContent();
		fp.close();
	}

	void saveBWT(std::string filename)
	{
		std::ofstream fp(filename, std::ios::binary);
		boost::archive::binary_oarchive archive_fp( fp );
		archive_fp & bwt;
		fp.close();
	}
	void readBWT(std::string filename)
	{
		std::ifstream fp(filename, std::ios::binary);
		boost::archive::binary_iarchive archive_fp( fp );
		archive_fp & bwt;
		fp.close();
	}

	void saveTable(std::string filename)
	{
		std::ofstream fp(filename, std::ios::binary);
		boost::archive::binary_oarchive archive_fp( fp );
		archive_fp & interval;
		archive_fp & c_function;
		archive_fp & c_functions;
		//archive_fp & c_function_sp;
		archive_fp & occ_function;
		archive_fp & location_table;
		archive_fp & fbwt;
		archive_fp & first_location;
		archive_fp & occ_jbwt;
		archive_fp & jbwt->seq_;
		archive_fp & jbwt_idx_char;
		fp.close();
	}
	void readTable(std::string filename)
	{
		std::ifstream fp(filename, std::ios::binary);
		boost::archive::binary_iarchive archive_fp( fp );
		archive_fp & interval;
		archive_fp & c_function;
		archive_fp & c_functions;
		//archive_fp & c_function_sp;
		archive_fp & occ_function;
		archive_fp & location_table;
		archive_fp & fbwt;
		archive_fp & first_location;
		archive_fp & occ_jbwt;
		archive_fp & jbwt_seq;
		archive_fp & jbwt_idx_char;
		fp.close();
		
		interval_set(interval);
	}
	
	
	
	
	template<class SEQTYPE>
	void createAllTable(SEQTYPE &seq, std::vector<std::string>& filenames)
	{
		jbwt = std::shared_ptr<JBit>( new JBit( seq.size() ) );
	
		occ_function.resize(char_size, std::vector<INTTYPE>() );//256
		location_table.reserve(seq.size()/interval+1);
		
		INTTYPE idx(0);
		char tmp_char('$'), bwt_char;
		std::string tmp_str, tmp_strs("");
		INTTYPE tmp_str_i(0), tmp_strs_i(0);
		
		INTTYPE tmp_cs(0);
		std::vector<INTTYPE> tmp_occ_count(char_size, 0);
		fbwt.resize( seq.size() ,0);
		
		INTTYPE c_functions_interval(12);
		
		c_functions.resize( std::pow(4,c_functions_interval) +1,0);
		
		//c_function_sp.resize( std::pow(4,c_functions_interval) ,0);
		
		//std::vector< std::pair<char, uint8_t> > test_compress;
		char tmp_c_char='\0';
		uint8_t tmp_c_int = 0;
		INTTYPE ttttt = 0;
		
		for(std::string& filename : filenames)
		{
			std::ifstream in (filename, std::ios::binary);
			boost::archive::binary_iarchive tmp_archive(in);	
			std::vector<INTTYPE> sorted_table;
			
			tmp_archive & sorted_table;
			in.close();
			
			for(INTTYPE c_idx : sorted_table)
			//for(INTTYPE c_idx : seq_table)
			//for(INTTYPE i(1); i<seq_table.size(); ++i) // no first line ($...)
			{	
				// BWT char
				if(c_idx == 0)
				{
					bwt_char = 'A'; //$
					first_location = idx;
				}
				else
					bwt_char = seq[c_idx-1];
				
				//bwt += bwt_char;
				jbwt->push_back(bwt_char);
				
				//fbwt
				fbwt[idx] = (c_idx & (interval-1));
				
				
				// C function
				if(seq[c_idx] != tmp_char)
				{
					tmp_char = seq[c_idx];
					std::cerr << "ooo - " << tmp_char << " " << (INTTYPE)tmp_char << " " << idx << std::endl;
					c_function[ (INTTYPE)tmp_char ] = idx;
				}
				
				//C functions
				
				//tmp_strs = seq.substr(c_idx, c_functions_interval);
				tmp_strs_i = c_idx;
				//if(tmp_strs != tmp_str && tmp_strs.size() == c_functions_interval)
				if( str_idx_compare(seq, tmp_strs_i, tmp_str_i, c_functions_interval) )
				{
					
					tmp_cs = 0;
					int is_exist_$(0);
					for(INTTYPE i(0); i < c_functions_interval; ++i)
					{
						//if(tmp_strs[i] == '$' )
						if(seq[tmp_strs_i+i] == '$')
						{
							is_exist_$ = 1;
							break;
						}
						
						//tmp_cs += mtable[ tmp_strs[ i ] ] * std::pow(4,(c_functions_interval-1-i)) ;
						tmp_cs += mtable[ seq[tmp_strs_i+i] ] * std::pow(4,(c_functions_interval-1-i)) ;
					}
					if( !is_exist_$ )
					{
						//c_function_sp[tmp_cs] = 
						c_functions[ tmp_cs ] = idx;
						tmp_str = tmp_strs;
						tmp_str_i = tmp_strs_i;
					}
						
				}
				
				// OCC function
				if( (idx & (interval-1) ) == 0)
				{
					//for (int v(0); v < occ_char.size(); ++v)
					//for(INTTYPE v(0); v < char_size; ++v)
					for (int v(0); v < 4; ++v)
					{
						occ_function[ occ_char[v] ].push_back(tmp_occ_count[ occ_char[v] ]);
						//occ_function[ v ].push_back(tmp_occ_count[v]);//==0? 0: tmp_occ_count[v]-1;
					}
				}
				tmp_occ_count[ bwt_char ] ++;
				
				//location table
				if( (c_idx & (interval-1) ) == 0)
				{
					location_table.push_back( {idx, c_idx} );
				}
				
				idx++;
			}
		}
		
		jbwt->last_push_back();
		
		c_function[c_function.size()-1] = seq.size();
		
		INTTYPE tmp( 0 );
		
		for(INTTYPE i(c_function.size()-1); i > 0; --i)
		{
			//std::cerr << "c function : i :" << i << " c : " << c_function[i] << std::endl;
			if(c_function[i] != 0) // magic number
				tmp = c_function[i];
			//else
			//	occ_function[i].clear();
			c_function[i] = tmp;
			//std::cerr << "c function : i :" << i << " c : " << c_function[i] << std::endl;
		}
		
		c_functions[c_functions.size()-1] = seq.size();
		
		tmp = c_functions[c_functions.size()-1];
		for(INTTYPE i(c_functions.size()-1); i > 0; --i)
		{
			if(c_functions[i] != 0)
				tmp = c_functions[i];
			c_functions[i] = tmp;
		}
		std::cerr << "Creating Jbwt..." << std::endl;
		using_jbwt();
	}


	
	template<class SEQTYPE>
	inline bool str_idx_compare(SEQTYPE &seq, INTTYPE a, INTTYPE b, INTTYPE len)
	{
		if(a+len >= seq.size() || b+len >= seq.size())
			return false;
		if(a==b)
			return false;
		for(INTTYPE i=0; i<len; i++)
		{
			if(seq[a+i] != seq[b+i])
				return true;
		}
		return false;
	}
	inline INTTYPE get_c(INTTYPE i) const
	{
		//clock_t start = clock();
		return c_function[ (INTTYPE)bwt[i] ];
		//clock_t end = clock();
		//get_c_time += double(end -start);
	}
	inline INTTYPE get_c(char c) const 
	{
		return c_function[ c ];
	}
	//
	inline char get_jbwt_char(INTTYPE i) const 
	{
		INTTYPE jbwt_idx( (i >> 2) );
		uint8_t chars = jbwt_seq[jbwt_idx];
		return jbwt_idx_char[ (chars << 2) + (i & 3) ];
		
	}
	INTTYPE get_occ_using_jbwt(INTTYPE i, char c = '\0', int show_error=0) const
	{
		
		//if(show_error)
		//	std::cerr <<"i : "<<i<< " , (i & (interval-1) ) : "<< (i & (interval-1) ) << " c " << c << std::endl;
		//clock_t start = clock();
		
		//if (c=='\0')
		//	c=bwt[i];
			
			
		INTTYPE pre_interval( (i >> powerV) );
		INTTYPE j( (pre_interval << powerV) );
		INTTYPE count( occ_function[ c ][ pre_interval ] );
		if(i > first_location &&	c =='A')
			--count;
		if(i == j)
			return count;
		
		INTTYPE jbwt_idx (j>>2), jbwt_idx_end(i>>2);
		
		INTTYPE tmp_count = count; 
		
		for( ; jbwt_idx != jbwt_idx_end; j+=4, ++jbwt_idx)
		{
			//count += occ_jbwt [ c ][4] [ jbwt->seq_[ jbwt_idx ] ];
			count += occ_jbwt [ c ][4] [ jbwt_seq[ jbwt_idx ] ];
			
			/*
			if(show_error)
			{
				for(auto kk=0;kk<4;kk++)
				{
					if (bwt[j+kk] == c )
						tmp_count++;
				}
				std::cerr << "jbwt_idx" << jbwt_idx << "j" << j << std::endl;
				if(tmp_count != count)
					std::cerr << "OH NO" << std::endl;
			}
			*/
		}
		/*
		if(show_error)
		{
			std::cerr << "O jbwt_idx" << jbwt_idx << "j" << j << std::endl;
			
			
			auto kk=j;
			tmp_count=0;
			for(;kk != i;kk++)
			{
				std::cerr << c <<std::endl;
				if (bwt[kk] == c )
					tmp_count++;
			}
			std::cerr << "j" << j << "i" << i << " seq " << "jbwt_idx" << jbwt_idx << " " <<  (int)jbwt->seq_[ jbwt_idx ]<<  std::endl;
			std::cerr << "ccccc" << occ_jbwt [ c ][i-j] [ jbwt->seq_[ jbwt_idx ] ] << " tmp_count " << tmp_count << std::endl;
			
		}
		*/
		
		//count += occ_jbwt [ c ][i-j] [ jbwt->seq_[ jbwt_idx ] ];
		count += occ_jbwt [ c ][i-j] [ jbwt_seq[ jbwt_idx ] ];
		
		return count;
			
		
	}
	
	inline INTTYPE get_occ(INTTYPE i, char c = '\0') const
	{
		//std::cerr <<"i : "<<i<< " , (i & (interval-1) ) : "<< (i & (interval-1) ) << std::endl;
		//clock_t start = clock();
		if (c=='\0')
			c=bwt[i];
			
			INTTYPE pre_interval( (i >> powerV) );
			INTTYPE j( (pre_interval << powerV) );
			INTTYPE count( occ_function[ (INTTYPE) c ][ pre_interval ] );
			if(i > first_location &&	c =='A')
				--count;
			if(i == j)
				return count;
			//std::cerr << "OOori count " << count << std::endl;
			
			for(; j != i; ++j)
			{
				if (bwt[j] == c )
					++count;
			}
			return count;	
	}
	
	INTTYPE back_tracking_using_jbwt(INTTYPE i) const
	{
		char c = get_jbwt_char(i);
		return get_c(c) + get_occ_using_jbwt(i, c);
	}
	inline INTTYPE back_tracking(INTTYPE i) const
	{
		//std::cerr << "i: " << i << "c: " << get_c(i) << " occ: " << get_occ(i) << std::endl;;
		return get_c(i) + get_occ(i);
	}
};

#endif
