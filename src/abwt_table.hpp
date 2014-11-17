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
#include "mkq_sort.hpp"
#include "split_sort.hpp"

#include "boost/serialization/utility.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/unordered_map.hpp"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>


#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/device/file.hpp"
#include "boost/iostreams/filter/zlib.hpp"
#include "boost/serialization/map.hpp"


class ABWT_table
{
public:

	typedef ABSequence<std::string> SEQTYPE;
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
	INTTYPE _realSize = 0;
	// starting site -> chromosome
	std::map <INTTYPE, std::string> chr_start_pos {};
	// chr -> chr size
	std::map <std::string, INTTYPE> chr_length {};
	// unambiguous segment sequence starting position of each chr
	std::map <INTTYPE, INTTYPE > chr_umbiguous_starting_length {};

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

	void readChrStartPos (const std::string& chrStartPosFile);

	void readChrLen (const std::string& chrLenFile);
	void readNPosLen (const std::string& fileName);
	void interval_set(INTTYPE iv);
	void using_jbwt();

//	template< class SEQTYPE>
	void saveSEQ(std::string filename, SEQTYPE &seq);

//	template< class SEQTYPE>
	void readSEQ(std::string filename, SEQTYPE &seq);
	void saveBWT(std::string filename);
	void readBWT(std::string filename);

	void saveTable(std::string filename);
	void readTable(std::string filename);


//	template<class SEQTYPE>
	void createAllTable(SEQTYPE &seq, std::vector<std::string>& filenames);


//	template<class SEQTYPE>
	void createAllTable(SEQTYPE &seq, std::vector<INTTYPE> &seq_table);


//	template<class SEQTYPE>
	inline bool str_idx_compare(SEQTYPE &seq, INTTYPE a, INTTYPE b, INTTYPE len);

	inline INTTYPE get_c(INTTYPE i) const;

	inline INTTYPE get_c(char c) const;

	inline char get_jbwt_char(INTTYPE i) const;

	INTTYPE get_occ_using_jbwt(INTTYPE i, char c='\0', int show_error=0) const;

	inline INTTYPE get_occ(INTTYPE i, char c='\0') const;

	INTTYPE back_tracking_using_jbwt(INTTYPE i) const;
	inline INTTYPE back_tracking(INTTYPE i) const;
};

#endif
