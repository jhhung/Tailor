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

#ifndef ABWT_SEARCH_HPP_
#define ABWT_SEARCH_HPP_

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
#include <boost/ref.hpp>
#include <ctime>
#include <cmath>
#include <memory>
#include <type_traits>
#include "compression/abit.hpp"
#include "compression/jbit.hpp"
#include "mkq_sort.hpp"
#include "difference_cover.hpp"
#include "split_sort.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/utility.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/unordered_map.hpp"
#include "abwt_format.hpp"
#include "constant_def.hpp"


template <class TABLE>
class ABWT_search
{
public:
	const TABLE& abwt_table_;
	mutable std::pair<INTTYPE, INTTYPE> start_end_pos_;
	std::array<INTTYPE,256> mtable_;
	INTTYPE c_functions_interval;
	std::array<char, 4> all_char;
	INTTYPE count;
	double binary_search_time;
	ABWT_search(TABLE& table)
		: abwt_table_ (table)
		, start_end_pos_()
		, binary_search_time()
		, c_functions_interval(12)
		, all_char{ {'A', 'C', 'G', 'T'} }
		,count(0)
	{
		mtable_['A'] = 0;
		mtable_['C'] = 1;
		mtable_['G'] = 2;
		mtable_['T'] = 3;
	}

	inline void init_exact_match( char c ) const
	{
		start_end_pos_.first = abwt_table_.c_function[ c ];
		start_end_pos_.second = abwt_table_.c_function[ c+1 ];
	}

	inline void init_exact_match( std::string tmp_str8 ) const
	{
		INTTYPE tmp_cs(0);
		for(INTTYPE i(0); i < c_functions_interval; ++i)
		{
			tmp_cs += mtable_[ tmp_str8[ i ] ] << ((c_functions_interval-1-i)<<1);
		}
		start_end_pos_.first = abwt_table_.c_functions[ tmp_cs ];
		start_end_pos_.second = abwt_table_.c_functions[ tmp_cs+1 ];
	}
	inline void init_exact_match( INTTYPE tmp_cs ) const
	{
		start_end_pos_.first = abwt_table_.c_functions[ tmp_cs ];
		start_end_pos_.second = abwt_table_.c_functions[ tmp_cs+1 ];
	}
	inline void exec_exact_match( char c ) const
	{
		start_end_pos_.first = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( start_end_pos_.first, c );
		start_end_pos_.second = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( start_end_pos_.second, c );
	}
	inline long find_nearest_mark_sbwt(INTTYPE trace_point) const
	{
		INTTYPE traceback_count(0);
		INTTYPE flocation = abwt_table_.fbwt[trace_point];
		for(uint8_t i(0); i<flocation;++i)
		{
			trace_point = abwt_table_.back_tracking_using_jbwt(trace_point);
			++traceback_count;
		}
		auto pp = std::lower_bound(
			abwt_table_.location_table.begin(),
			abwt_table_.location_table.end(),
			std::pair<INTTYPE, INTTYPE>{trace_point, 0},
				[]( const std::pair<INTTYPE, INTTYPE>& a, const std::pair<INTTYPE, INTTYPE>& b)
				{
					return a.first < b.first;
				}
		);
		if ( pp->first == trace_point )
		{
			return traceback_count + pp->second;
		}
		return -1;
	}
	inline INTTYPE find_nearest_mark (INTTYPE trace_point) const
	{
		INTTYPE traceback_count(0);
		INTTYPE flocation = abwt_table_.fbwt[trace_point];

		for(uint8_t i(0); i<flocation;++i)
		{
			trace_point = abwt_table_.back_tracking_using_jbwt(trace_point);
			++traceback_count;
		}

		auto pp = std::lower_bound(
			abwt_table_.location_table.begin(),
			abwt_table_.location_table.end(),
			std::pair<INTTYPE, INTTYPE>{trace_point, 0},
				[]( const std::pair<INTTYPE, INTTYPE>& a, const std::pair<INTTYPE, INTTYPE>& b)
				{
					return a.first < b.first;
				}
		);
		traceback_count += pp->second;
		return traceback_count;

		while(1)//trace to the nearest upstream location_table
		{
			auto pp = std::lower_bound(
				abwt_table_.location_table.begin(),
				abwt_table_.location_table.end(),
				std::pair<INTTYPE, INTTYPE>{trace_point, 0},
					[]( const std::pair<INTTYPE, INTTYPE>& a, const std::pair<INTTYPE, INTTYPE>& b)
					{
						return a.first < b.first;
					}
			);

			if ( pp->first == trace_point )
			{
				traceback_count += pp->second;
				break;
			}
			else//no hit yet
			{
				trace_point = abwt_table_.back_tracking_using_jbwt(trace_point);
				traceback_count++;
			}
		}
		return traceback_count;
	}

	inline void start_all_possible(std::string& query, std::vector<INTTYPE>& result)
	{

		init_exact_match( query.substr(query.size()-c_functions_interval,c_functions_interval) );

		INTTYPE first_pos_f(start_end_pos_.first);
		INTTYPE first_pos_e(start_end_pos_.second);
		all_possible(result, query, query.length()-1-c_functions_interval, first_pos_f, first_pos_e, 0);
	}

	inline void all_possible(std::vector<INTTYPE>& result, std::string &query, int i, INTTYPE pos_f, INTTYPE pos_e, INTTYPE mn)
	{
		char c(query[i]);

		INTTYPE n_pos_f = abwt_table_.c_function[ c ] + abwt_table_.get_occ( pos_f, c );
		INTTYPE n_pos_e = abwt_table_.c_function[ c ] + abwt_table_.get_occ( pos_e, c );

		if(i == 0)
		{
			push_result(result, n_pos_f, n_pos_e);
			return;
		}
		if(pos_f < pos_e)
		{
			all_possible(result, query, i-1, n_pos_f, n_pos_e, mn);
		}
		if(mn == 7 )
			return;

		for(int cn(0); cn < all_char.size(); ++cn)
		{
			if( c == all_char[cn])
				continue;

			query[i] = all_char[cn];
			all_possible(result, query, i, pos_f, pos_e, mn+1);
			query[i] = c;
		}
	}

	inline void start_one_mismatch(std::string& query, std::vector<INTTYPE>& result)
	{
		init_exact_match( query.substr(query.size()-c_functions_interval,c_functions_interval) );
		INTTYPE first_pos_f(start_end_pos_.first);
		INTTYPE first_pos_e(start_end_pos_.second);

		INTTYPE kk(0);
		for (int i=query.length()-1-c_functions_interval; i>=0 && start_end_pos_.first < start_end_pos_.second; i--)

		{
			exec_exact_match(query[i]);
		}
		//std::cout << "Aresult.size()" << result.size() << " " << start_end_pos_.second - start_end_pos_.first << std::endl;}

		push_result(result, start_end_pos_.first, start_end_pos_.second);
		//std::cout << "Bresult.size()" << result.size() << std::endl;
		if(result.size() != 0)
			return;
		//std::cout << "exact_match_rec" << std::endl;
		exact_match_rec(result, query, query.length()-1-c_functions_interval, first_pos_f, first_pos_e, 0);

	}

	inline void start_one_mismatch2(std::string& query, std::vector<INTTYPE>& result)
	{

		init_exact_match( query.substr(query.size()-c_functions_interval,c_functions_interval) );

		INTTYPE first_pos_f(start_end_pos_.first);
		INTTYPE first_pos_e(start_end_pos_.second);

		INTTYPE kk(0);
		for (int i=query.length()-1-c_functions_interval; i>=0 && start_end_pos_.first < start_end_pos_.second; i--)
		{
			exec_exact_match(query[i]);

		}
		push_result(result, start_end_pos_.first, start_end_pos_.second);
		if(result.size() != 0)
			return;
		exact_match_rec(result, query, query.length()-1-c_functions_interval, first_pos_f, first_pos_e, 0);
	}
	inline void exact_match_rec(std::vector<INTTYPE>& result, std::string &query, int i, INTTYPE pos_f, INTTYPE pos_e, INTTYPE mn)
	{
		// now i
		char c(query[i]);
		INTTYPE n_pos_f = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_f, c );
		INTTYPE n_pos_e = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_e, c );

		if(i == 0)//tail 也沒中，表示seed後都沒有中，可以略過 seed 之後
		{
			push_result(result, n_pos_f, n_pos_e);
			return;
		}
		if(pos_f < pos_e) //y, yes match
		{
			//tail 也沒中，表示seed後都沒有中，可以略過 seed 之後
			exact_match_rec(result, query, i-1, n_pos_f, n_pos_e, mn);
		}
		//having matches, not to mismatch, break
		if(result.size() != 0)
		{
			return;
		}
		//misamtch number limit
		if(mn == 1)
		{
			return;
		}
		for(int cn(0); cn < all_char.size(); ++cn)
		{
			if( query[i+1] == all_char[cn])
				continue;
			query[i] = all_char[cn];
			exact_match_rec(result, query, i, pos_f, pos_e, mn+1);
			query[i] = c;
		}
	}

	//(results, fq, query_, query.length()-1-1, pos_f, pos_e, 0, prefixMatchLen);
	//result:: pos, mismatch_pos_i, char, tail_pos
	inline void tailor_mismatch_rec
		( std::vector< std::tuple<INTTYPE, int, char, int> >& result
		, std::string &query
		, INTTYPE i
		, INTTYPE pos_f
		, INTTYPE pos_e
		, int mn
		, int minimalPrefixLength
		, int mismatch_pos_i
		, char mismatch_char
		)
	{
		// now i
		char c(query[i]);
		INTTYPE n_pos_f = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_f, c );
		INTTYPE n_pos_e = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_e, c );

		//之前已經做過perfact match，所以可以省略沒有 mismatch還是會全部做完，只要超過prefixMatchLen 就略過
		//tail 也沒中，表示seed後都沒有中，可以略過 seed 之後
		if(mn == 0 && i < query.size()-minimalPrefixLength)
			return;

		//不管有沒有mismatch，i==0就是已經結束，0 mismatch or 1 mismatch
		if(i == 0)
		{
			for (INTTYPE pos_idx = n_pos_f; pos_idx < n_pos_e; pos_idx++)
			{
				auto real_pos = find_nearest_mark(pos_idx);
				result.push_back(std::make_tuple(real_pos, mismatch_pos_i, mismatch_char, 0));
			}
			return;
		}
		if(pos_f < pos_e) //y, yes old have range, do next pos
		{
			//tail 也沒中，表示seed後都沒有中，可以略過 seed 之後
			tailor_mismatch_rec(result, query, i-1, n_pos_f, n_pos_e, mn, minimalPrefixLength, mismatch_pos_i, mismatch_char);
		}
		//there are matches, not to do more mismatch, break
		if(result.size() != 0)
		{
			return;
		}
		//misamtch number limit, maybe a tail ( 2 mismatch)
		if(mn == 1)
		{
			// next pos 都沒有 range了，i 位置已經是錯的
			// 超過 min prefix match len 外，可以是 tail
			if(i < query.size()-minimalPrefixLength)
			{
				for (INTTYPE pos_idx = n_pos_f; pos_idx < n_pos_e; pos_idx++)
				{
					auto real_pos = find_nearest_mark(pos_idx);
					result.push_back(std::make_tuple(real_pos, mismatch_pos_i, mismatch_char, i));
				}

			}
			return;
		}
		for(int cn(0); cn < all_char.size(); ++cn)
		{
			if( query[i] == all_char[cn])
				continue;
			query[i] = all_char[cn];
			tailor_mismatch_rec(result, query, i, pos_f, pos_e, mn+1, minimalPrefixLength, i, all_char[cn]);
			query[i] = c;
		}
	}
	inline void position_to_sam(std::stringstream* out, std::vector< std::tuple<INTTYPE, int, char, int> >& results, const Fastq& fq, std::string &query, int prefixMatch)
	{
		//pos, mismatch_pos_i, char, tail_pos
		for(auto &result : results)
		{
			auto position = std::get <0>(result);
			auto mismatch_pos_i = std::get <1>(result);
			auto mismatch_char = std::get <2>(result);
			auto tail_pos = std::get <3>(result);

			int queryPosition = tail_pos;
			int prefixMatchLen = query.size() - queryPosition;
			INTTYPE NH_tag = results.size();
			std::string _query = fq.getSeq();

			bool isRC = false;
			/// the second comparsion is to suppress weird bug of TTTTTTTTTTTT mapping to position == 2*abwt_table_._realSize
			if (position >= this->abwt_table_._realSize && position < (abwt_table_._realSize<<1))
			{
				isRC = true;
				position = this->abwt_table_._realSize*2 - position - prefixMatchLen;
			} else if (position < this->abwt_table_._realSize) {
				isRC = false;
			} else {
				continue;
			}
			auto tailSeq = _query.substr(prefixMatchLen);
			auto lowerIter = this->abwt_table_.chr_start_pos.upper_bound (position);
			std::advance (lowerIter, -1);
			auto chr = lowerIter->second;
			auto lowerIter3 = this->abwt_table_.chr_start_pos.upper_bound (position + prefixMatchLen -1);
			std::advance (lowerIter3, -1);
			auto chr3 = lowerIter3->second;
			if (chr != chr3) continue;
			auto NLowerIter = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position);
			std::advance (NLowerIter, -1);
			auto NLowerIter3 = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position + prefixMatchLen -1);
			std::advance (NLowerIter3, -1);
			if (NLowerIter != NLowerIter3) continue;
			position = position - lowerIter->first + NLowerIter->second;

			std::string MD_tag("");
			std::string CIGAR_tag("");

			//反
			if (!isRC)
			{ /// same as start_tailing_match_AS
				//TODO: redefine MAPQ
				if(queryPosition == 0)
					CIGAR_tag = std::to_string (prefixMatchLen) + 'M';
				else
					CIGAR_tag = std::to_string (queryPosition) + 'S' + std::to_string (prefixMatchLen) + 'M';

				// 有 mismatch 也有 tail
				int pre = query.size() - mismatch_pos_i - 1;
				int pro = mismatch_pos_i - tail_pos;

				//if(pro!=0)
				MD_tag += std::to_string(pro);
				MD_tag += mismatch_char;
				//if(pre!=0)
				MD_tag += std::to_string(pre);
				//MD_tag += tailSeq;

				//if(pre == 0)
				//{
					//mismatch 在最後一個字，屬於 tail
					//MD_tag = "";
					//tailSeq = mismatch_char;
					//CIGAR_tag = std::string("1S") + std::to_string (prefixMatchLen-1) + 'M';
				//}
				//std::cout << " A position " << position+1 << " MD_tag " << MD_tag << std::endl;
				*out << Sam { fq.getName (),
					Sam::SAM_FLAG::REVERSE_COMPLEMENTED,
					std::move (chr),
					position + 1,
					255 - (queryPosition),
					std::move(CIGAR_tag),
					"*",
					0,
					0,
					query,
					fq.getRevQuality (),
					NH_tag,
					std::move(tailSeq),
					std::move(MD_tag)};
			}
			else
			{ /// same as start_tailing_match_S
				//TODO: redefine MAPQ
				switch (mismatch_char)
				{
					case 'A': mismatch_char = 'T'; break;
					case 'T': mismatch_char = 'A'; break;
					case 'C': mismatch_char = 'G'; break;
					case 'G': mismatch_char = 'C'; break;
				}
				if(queryPosition == 0)
					CIGAR_tag = std::to_string (prefixMatchLen) + 'M';
				else
					CIGAR_tag = std::to_string (prefixMatchLen) + 'M' + std::to_string (queryPosition) + 'S';

				int pre = query.size() - mismatch_pos_i - 1;
				int pro = mismatch_pos_i - tail_pos;
				//if(pre!=0)
				MD_tag += std::to_string(pre);
				MD_tag += mismatch_char;
				//if(pro!=0)
				MD_tag += std::to_string(pro);
				//MD_tag += tailSeq;

				*out << Sam {fq.getName (),
					Sam::SAM_FLAG::MAPPED,
					std::move (chr),
					position+1,
					255 - (queryPosition) ,
					std::move(CIGAR_tag),
					"*",
					0,
					0,
					_query,
					fq.getQuality (),
					NH_tag,
					std::move(tailSeq),
					std::move(MD_tag)};
			}
		}

	}

	inline void exact_match_rec(std::vector<INTTYPE>& result, std::string &query, int i, INTTYPE pos_f, INTTYPE pos_e, INTTYPE b_pos_f, INTTYPE b_pos_e, INTTYPE mn)
	{
		char c(query[i]);

		INTTYPE n_pos_f = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_f, c );
		INTTYPE n_pos_e = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_e, c );

		if(i == 0)
		{
			push_result(result, n_pos_f, n_pos_e);
			return;
		}
		if(pos_f < pos_e)
		{
			exact_match_rec(result, query, i-1, n_pos_f, n_pos_e, pos_f, pos_e, mn);
		}
	}

	inline void push_result(std::vector<INTTYPE>& result, INTTYPE pos_f, INTTYPE pos_e)
	{
		for (long i = pos_f; i < pos_e; i++)
		{
			result.push_back( find_nearest_mark(i) );
		}
	}

	inline void start_sbwt_match( std::string& query, std::vector<INTTYPE>& result )
	{
		init_exact_match( query.substr(query.size()-c_functions_interval,c_functions_interval) );

		for (int i=query.length()-1-c_functions_interval; i>=0 && start_end_pos_.first < start_end_pos_.second; i--)
		{
			exec_exact_match(query[i]);
		}
		if (start_end_pos_.second - start_end_pos_.first > 100)
			return;

		for(INTTYPE i = start_end_pos_.first; i < start_end_pos_.second; ++i)
		{
			if(abwt_table_.fbwt[i] != 0)
				continue;
			auto pp = std::lower_bound(
				abwt_table_.location_table.begin(),
				abwt_table_.location_table.end(),
				std::pair<INTTYPE, INTTYPE>{ i, 0},
					[]( const std::pair<INTTYPE, INTTYPE>& a, const std::pair<INTTYPE, INTTYPE>& b)
					{
						return a.first < b.first;
					}
			);
			if(pp->first == i)
			{
				result.push_back(pp->second);
				if(start_end_pos_.second - start_end_pos_.first == 1)
					return;
			}

		}

		for(int cn(0); cn < all_char.size(); ++cn)
		{
			find_possible(result, 1, all_char[cn], start_end_pos_.first, start_end_pos_.second);
		}

	}

	inline void find_possible(std::vector<INTTYPE>& result, int len, char c, INTTYPE pos_f, INTTYPE pos_e)
	{
		if(len >= abwt_table_.interval)
			return;
		if (result.size() == start_end_pos_.second - start_end_pos_.first)
			return;

		INTTYPE n_pos_f = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_f, c );
		INTTYPE n_pos_e = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_e, c );
		if(n_pos_f >= n_pos_e) // no result
		{
			// no result
			return;
		}

		for(INTTYPE i = n_pos_f; i < n_pos_e; ++i)
		{
			if(abwt_table_.fbwt[i] != 0)
				continue;
			auto pp = std::lower_bound(
				abwt_table_.location_table.begin(),
				abwt_table_.location_table.end(),
				std::pair<INTTYPE, INTTYPE>{ i, 0},
					[]( const std::pair<INTTYPE, INTTYPE>& a, const std::pair<INTTYPE, INTTYPE>& b)
					{
						return a.first < b.first;
					}
			);
			if(pp->first == i)
			{
				result.push_back(len + pp->second);
				if(n_pos_e - n_pos_f == 1 )
					return;
			}

			if(len == abwt_table_.interval-1)
				return;

		}
		if(len == abwt_table_.interval-1)
			return;
		for(int cn(0); cn < all_char.size(); ++cn)
		{
			find_possible(result, len+1, all_char[cn], n_pos_f, n_pos_e);
		}
	}

	inline void start_exact_match( std::string& query, std::vector<INTTYPE>& result )
	{

		init_exact_match( query.substr(query.size()-c_functions_interval,c_functions_interval) );


		for (int i=query.length()-1-c_functions_interval; i>=0 && start_end_pos_.first < start_end_pos_.second; i--)
		{
			exec_exact_match(query[i]);
		}
		if (start_end_pos_.second - start_end_pos_.first > 100)
			return;
		for (INTTYPE i = start_end_pos_.first; i < start_end_pos_.second; i++)
		{
			result.push_back( find_nearest_mark(i) );
		}
	}




	/*@ a version of exec_exact_match that return the previous value of start_end_pos_ */
	std::pair<INTTYPE, INTTYPE> exec_exact_match2 ( char c ) const {
		auto old_start_end_pos_ = start_end_pos_;
		start_end_pos_.first = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( start_end_pos_.first, c );
		start_end_pos_.second = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( start_end_pos_.second, c );
		return old_start_end_pos_;
	}

	// tailing searching version for dual BWT
    template <class allowMismatch>
	void start_tailing_match_Dual (const Fastq& fq, std::stringstream* out, int minimalPrefixLen)//, int maximumTailLen=3)
	{
		std::string _query = fq.getSeq ();

		int minimalPrefixLength = minimalPrefixLen;
		//int maximumTailLength = maximumTailLen;
		//minimalPrefixLength = std::max( (int)(_query.size() - maximumTailLength), minimalPrefixLength);

		/* reverse complement query string */
		std::string query {_query.crbegin(), _query.crend()};
		for (auto & c : query) {
			switch (c) {
				case 'A': c = 'T'; break;
				case 'T': c = 'A'; break;
				case 'C': c = 'G'; break;
				case 'G': c = 'C'; break;
				default : return;
			}
		} /* end of RC */
		bool isRC = false;
		/* the front c_functions_interval has to be exact match */
		//init_exact_match( query.substr (query.size() - c_functions_interval, c_functions_interval) );
		init_exact_match(query.back());
		
		//std::cout << " GG " << start_end_pos_.second - start_end_pos_.first << " " << fq.getName() << std::endl;
		/* if not found, exit */
		if ( start_end_pos_.first >= start_end_pos_.second) {
			return;
		}
		/* queryPosition need to be tested again */
		//int queryPosition = query.length() -1 - c_functions_interval;
		int queryPosition = query.length() -1 - 1;
		/* last step need to be recorded */
		std::vector< std::pair<INTTYPE, INTTYPE> > last_start_end_pos_record;
		
		auto last_start_end_pos_ = start_end_pos_;
		
		last_start_end_pos_record.push_back(last_start_end_pos_);
		for (; queryPosition >= 0 && start_end_pos_.first < start_end_pos_.second; --queryPosition) {
			last_start_end_pos_ = this->exec_exact_match2(query[queryPosition]);
			last_start_end_pos_record.push_back(last_start_end_pos_);
		}
		/// substract an extra one when exiting the loop, so add it back
		auto prefixMatchLen = _query.size() - 1 - (queryPosition + 1);
		//std::cout << "mismatch " << allowMismatch::value << " " << prefixMatchLen+1 << std::endl;
		if (prefixMatchLen+1 < minimalPrefixLength && allowMismatch::value )
		{
			//std::cout << "do mismatch" << std::endl;
			//return;
// too many mismatch
// DO mismatch
			//no tail and do mismatch
			std::vector< std::tuple<INTTYPE, int, char, int> > results;
			// 取的最大的範圍 從 0 到 A
			init_exact_match( 'A' );
			INTTYPE pos_f(start_end_pos_.first);
			init_exact_match( 'T' );
			INTTYPE pos_e(start_end_pos_.second);

			//std::vector<INTTYPE>& result, std::string &query, INTTYPE i, INTTYPE pos_limit, INTTYPE pos_f, INTTYPE pos_e, INTTYPE mn
			this->tailor_mismatch_rec(results, query, query.length()-1, pos_f, pos_e, 0, minimalPrefixLength, query.length(), 'A');
			//this->tailor_mismatch_rec(results, query, query.length()-1-1, pos_f, pos_e, 0, minimalPrefixLength, query.length(), 'A');
			//position_to_sam(std::vector< std::tuple<INTTYPE, int, char, int> >& results, const Fastq& fq, std::string &query, int prefixMatch)
			this->position_to_sam(out, results, fq, query, prefixMatchLen);

			return;
		}
		//if (prefixMatchLen+1 < minimalPrefixLength)
		//	return;
		
		
		//std::cout << start_end_pos_.second - start_end_pos_.first << " " << fq.getName() << " queryPosition "<< queryPosition << " prefixMatchLen "<< prefixMatchLen << std::endl;

		bool is_any_perfect = false;
		/// found perfect match
		//if (queryPosition == -1) {
		if(start_end_pos_.first < start_end_pos_.second) {
			auto NH_tag = start_end_pos_.second - start_end_pos_.first;
			for (INTTYPE i = start_end_pos_.first; i < start_end_pos_.second; i++) {
				auto position = this->find_nearest_mark(i);
				if (position >= this->abwt_table_._realSize) {
					isRC = true;
					position = this->abwt_table_._realSize*2 - position - _query.size ();
				} else {
					isRC = false;
				}
				
				if (!isRC) {
					auto lowerIter = this->abwt_table_.chr_start_pos.upper_bound (position);
					std::advance (lowerIter, -1);
					auto chr = lowerIter->second;

					auto lowerIter3 = this->abwt_table_.chr_start_pos.upper_bound (position + _query.size () - 1);
					std::advance (lowerIter3, -1);
					auto chr3 = lowerIter3->second;
					if (chr != chr3) continue;

					auto NLowerIter = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position);
					std::advance (NLowerIter, -1);

					auto NLowerIter3 = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position + _query.size () - 1);
					std::advance (NLowerIter3, -1);
					if (NLowerIter != NLowerIter3) continue;

					position = position - lowerIter->first + NLowerIter->second;
					*out << Sam {fq.getName (),
						Sam::SAM_FLAG::REVERSE_COMPLEMENTED,
						std::move (chr),
						position+1,
						255,
						std::to_string (_query.size ()) + 'M',
						"*",
						0,
						0,
						query,
						fq.getRevQuality (),
						NH_tag
					};
					is_any_perfect = true;
				} else {
					auto lowerIter = this->abwt_table_.chr_start_pos.upper_bound (position);
					std::advance (lowerIter, -1);
					auto chr = lowerIter->second;
					auto lowerIter3 = this->abwt_table_.chr_start_pos.upper_bound (position + _query.size () - 1);
					std::advance (lowerIter3, -1);
					auto chr3 = lowerIter3->second;
					if (chr != chr3) continue;
					auto NLowerIter = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position);
					std::advance (NLowerIter, -1);
					auto NLowerIter3 = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position + _query.size () - 1);
					std::advance (NLowerIter3, -1);
					
					if (NLowerIter != NLowerIter3) continue;
					position = position - lowerIter->first + NLowerIter->second;

					*out << Sam {fq.getName (),
						Sam::SAM_FLAG::MAPPED,
						std::move (chr),
						position+1,
						255,
						std::to_string (_query.size ()) + 'M',
						"*",
						0,
						0,
						_query,
						fq.getQuality (),
						NH_tag
					};
					is_any_perfect = true;
				}
			}
			//return ;
		}

/// begin recording tailing
		if (start_end_pos_.first >= start_end_pos_.second || !is_any_perfect) {

			++queryPosition; /// substract an extra one when exiting the loop, so add it back

			//測試是不是真的 tail，還是只是 mismatch
			if(queryPosition != 0  && allowMismatch::value)
			{
				int result_size = 0;
				for(int cn(0); cn < all_char.size(); ++cn)
				{
					auto tmp_i = queryPosition;
					//
					if( query[tmp_i] == all_char[cn])
						continue;
					start_end_pos_ = last_start_end_pos_;
					auto test_start_end_pos_ = this->exec_exact_match2(all_char[cn]);
					--tmp_i;
					for (; tmp_i >= 0 && test_start_end_pos_.first < test_start_end_pos_.second; --tmp_i)
					{
						test_start_end_pos_ = this->exec_exact_match2(query[tmp_i]);
					}
					if(tmp_i == -1)
					{
						std::vector< std::tuple<INTTYPE, int, char, int> > results2;
						// 只是 mismatch
						for (INTTYPE i = start_end_pos_.first; i < start_end_pos_.second; i++)
						{
							auto position = this->find_nearest_mark(i);
							results2.push_back( std::make_tuple(position, queryPosition, all_char[cn], 0) );
						}
						this->position_to_sam(out, results2, fq, query, prefixMatchLen);
						result_size += results2.size();
					}
				}
				if(result_size != 0)
					return;
			}

			//std::cout << " last_start_end_pos_record.size() " <<  last_start_end_pos_record.size() << std::endl;
			bool is_any_result = false;
			prefixMatchLen++;
			queryPosition--;
			
			/// some read with tail will align to like: 
			// ATCGGTAA (with N: ATCGNNNNNGTAC)
			// ATCGGCCC (Bug, ATCGG "last_start_end_pos_record"  record second G, but in genome, prefix match is ATCG, because of NNNNN...)
			// before, did not record "last_start_end_pos_record" with fitst G, so debug: make "last_start_end_pos_record", record all the last_start_end_poses.
			for(auto idx = last_start_end_pos_record.size(); idx >= (minimalPrefixLength-c_functions_interval); idx--)
			{
				prefixMatchLen--;
				queryPosition++;
				
				last_start_end_pos_ = last_start_end_pos_record[idx-1];
				
				auto NH_tag = last_start_end_pos_.second - last_start_end_pos_.first; // record the theoretically hits
				// tail後還有 mismatches，判定為真 tail
				for (INTTYPE i = last_start_end_pos_.first; i < last_start_end_pos_.second; i++) {
					auto position = this->find_nearest_mark(i);
					
					if (position >= this->abwt_table_._realSize && position < (abwt_table_._realSize<<1)) { /// the second comparsion is to suppress weird bug of TTTTTTTTTTTT mapping to position == 2*abwt_table_._realSize
						isRC = true;
						position = this->abwt_table_._realSize*2 - position - prefixMatchLen;
					} else if (position < this->abwt_table_._realSize) {
						isRC = false;
					} else {
						continue;
					}
					if (!isRC) { /// same as start_tailing_match_AS
						auto tailSeq = _query.substr(prefixMatchLen);
						
						auto lowerIter = this->abwt_table_.chr_start_pos.upper_bound (position);
						std::advance (lowerIter, -1);
						auto chr = lowerIter->second;
						auto lowerIter3 = this->abwt_table_.chr_start_pos.upper_bound (position + prefixMatchLen -1);
						std::advance (lowerIter3, -1);
						auto chr3 = lowerIter3->second;
						if (chr != chr3) continue;
	
						auto NLowerIter = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position);
						std::advance (NLowerIter, -1);
						auto NLowerIter3 = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position + prefixMatchLen -1);
						std::advance (NLowerIter3, -1);
						if (NLowerIter != NLowerIter3) continue;
						
						position = position - lowerIter->first + NLowerIter->second;
						
						//TODO: redefine MAPQ
						if( minimalPrefixLen > prefixMatchLen)
						{	
							is_any_result = true;
							continue;
						}
						else;
						*out << Sam { fq.getName (),
							Sam::SAM_FLAG::REVERSE_COMPLEMENTED,
							std::move (chr),
							position + 1,
							255 - queryPosition - 1,
							std::to_string (queryPosition+1) + 'S' + std::to_string (prefixMatchLen) + 'M',
							"*",
							0,
							0,
							query,
							fq.getRevQuality (),
							NH_tag,
							std::move(tailSeq)};
						is_any_result = true;
					} else { /// same as start_tailing_match_S
						auto tailSeq = _query.substr(prefixMatchLen);
						
						auto lowerIter = this->abwt_table_.chr_start_pos.upper_bound (position);
						std::advance (lowerIter, -1);
						auto chr = lowerIter->second;
						auto lowerIter3 = this->abwt_table_.chr_start_pos.upper_bound (position + prefixMatchLen -1);
						std::advance (lowerIter3, -1);
						auto chr3 = lowerIter3->second;
						if (chr != chr3) continue;
	
						auto NLowerIter = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position);
						std::advance (NLowerIter, -1);
						auto NLowerIter3 = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position + prefixMatchLen -1);
						std::advance (NLowerIter3, -1);
						
						if (NLowerIter != NLowerIter3)
							continue;
						
						position = position - lowerIter->first + NLowerIter->second;
						
						if( minimalPrefixLen > prefixMatchLen)
						{
							is_any_result = true;
							continue;
						}
						else;
						//TODO: redefine MAPQ
						*out << Sam {fq.getName (),
							Sam::SAM_FLAG::MAPPED,
							std::move (chr),
							position+1,
							255 - queryPosition - 1,
							std::to_string (prefixMatchLen) + 'M' + std::to_string (queryPosition+1) + 'S',
							"*",
							0,
							0,
							_query,
							fq.getQuality (),
							NH_tag,
							std::move(tailSeq)};
						is_any_result = true;
					}//if
					
				}//for
				if(is_any_result)
					break;
			}//for
			
			return;
		}

	}
};


#endif
