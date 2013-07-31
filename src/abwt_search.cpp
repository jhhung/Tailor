/// @file abwt_search.hpp
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
#include "difference_cover.hpp"
#include "split_sort.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/utility.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/unordered_map.hpp"
//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/archive/binary_oarchive.hpp>

#define INTTYPE uint64_t
#include "abwt_format.hpp"

template <class TABLE>
class ABWT_search
{
public:
	const TABLE& abwt_table_;
	mutable std::pair<INTTYPE, INTTYPE> start_end_pos_;
	//boost::unordered_map<char, INTTYPE> mtable_;
	std::array<INTTYPE,256> mtable_;
	INTTYPE c_functions_interval;
	std::array<char, 4> all_char;
	
	//std::vector<int> keys;
	
	INTTYPE count;
	
	double binary_search_time;
	ABWT_search(TABLE& table)
		: abwt_table_ (table)
		, start_end_pos_() //initialize as 0,0 by default
		//, mtable_({{'A',0},{'C',1},{'G',2},{'T',3}})
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
		//
		INTTYPE tmp_cs(0);
		for(INTTYPE i(0); i < c_functions_interval; ++i)
		{
			tmp_cs += mtable_[ tmp_str8[ i ] ] << ((c_functions_interval-1-i)<<1); //slower.....@@
			//tmp_c8 += mtable_[ tmp_str8[ i ] ] * std::pow(mtable_.size(),(7-i)) ;
			//std::cerr << tmp_str8[ i ] << " * " << std::pow(mtable_.size(),(7-i)) << " = " << tmp_c8 << std::endl;
		}
		//std::cerr << "tmp_cs " << tmp_cs << " size " << abwt_table_.c_functions.size() << std::endl;
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
		//INTTYPE tmp_a = start_end_pos_.first;
		//INTTYPE tmp_b = start_end_pos_.second;
		
		//INTTYPE a = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( start_end_pos_.first, c );
		//INTTYPE b = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( start_end_pos_.second, c );
		
		//start_end_pos_.first = abwt_table_.c_function[ c ] + abwt_table_.get_occ( start_end_pos_.first, c );
		//start_end_pos_.second = abwt_table_.c_function[ c ] + abwt_table_.get_occ( start_end_pos_.second, c );

		
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
		
		//INTTYPE tmp = trace_point;
		INTTYPE traceback_count(0);
		INTTYPE flocation = abwt_table_.fbwt[trace_point];
		
		for(uint8_t i(0); i<flocation;++i)
		{
			//trace_point = abwt_table_.back_tracking(trace_point);
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
			
			//std::cerr << "last trace_point : " << trace_point << " flocation : " << flocation << std::endl;
			
			auto pp = std::lower_bound( 
				abwt_table_.location_table.begin(), 
				abwt_table_.location_table.end(), 
				std::pair<INTTYPE, INTTYPE>{trace_point, 0}, 
					[]( const std::pair<INTTYPE, INTTYPE>& a, const std::pair<INTTYPE, INTTYPE>& b)
					{
						return a.first < b.first;
					} 
			);
			
			///std::cerr << "trace_point : " << trace_point << " f " << pp->first << " s " << pp->second << std::endl;




			if ( pp->first == trace_point )
			{
				//std::cerr << "traceback_count : " << traceback_count << std::endl;

				//traceback_count += hit.first->second;
				traceback_count += pp->second;
				break;
			}
			else//no hit yet
			{
				//trace_point = abwt_table_.back_tracking(trace_point);
				trace_point = abwt_table_.back_tracking_using_jbwt(trace_point);
				traceback_count++;

			}

		}
		return traceback_count;
	}
	
	inline void start_all_possible(std::string& query, std::vector<INTTYPE>& result)
	{	
			
		init_exact_match( query.substr(query.size()-c_functions_interval,c_functions_interval) );
		//init_exact_match( query.back() );
		
		INTTYPE first_pos_f(start_end_pos_.first);
		INTTYPE first_pos_e(start_end_pos_.second);

		
		all_possible(result, query, query.length()-1-c_functions_interval, first_pos_f, first_pos_e, 0);
		
		//push_result();
	}
	
	inline void all_possible(std::vector<INTTYPE>& result, std::string &query, int i, INTTYPE pos_f, INTTYPE pos_e, INTTYPE mn)
	{
		//std::cerr << "i" << i << std::endl;
		char c(query[i]);
		
		//INTTYPE n_pos_f = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_f, c );
		//INTTYPE n_pos_e = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_e, c );
		INTTYPE n_pos_f = abwt_table_.c_function[ c ] + abwt_table_.get_occ( pos_f, c );
		INTTYPE n_pos_e = abwt_table_.c_function[ c ] + abwt_table_.get_occ( pos_e, c );
		
		if(i == 0)
		{
			push_result(result, n_pos_f, n_pos_e);
			return;
		}
		if(pos_f < pos_e) //y
		{
			all_possible(result, query, i-1, n_pos_f, n_pos_e, mn);
		}
		if(mn == 7 )//1)
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
		//init_exact_match( query.back() );

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

		//push_result();
	}
		
	inline void start_one_mismatch2(std::string& query, std::vector<INTTYPE>& result)
	{	
			
		init_exact_match( query.substr(query.size()-c_functions_interval,c_functions_interval) );
		//init_exact_match( query.back() );
		
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
		
		//push_result();
	}
	inline void exact_match_rec(std::vector<INTTYPE>& result, std::string &query, int i, INTTYPE pos_f, INTTYPE pos_e, INTTYPE mn)
	{
		char c(query[i]);
		
		INTTYPE n_pos_f = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_f, c );
		INTTYPE n_pos_e = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_e, c );
		
		if(i == 0)
		{
			push_result(result, n_pos_f, n_pos_e);
			return;
		}
		if(pos_f < pos_e) //y
		{
			exact_match_rec(result, query, i-1, n_pos_f, n_pos_e, mn);
		}

		if(mn == 1)
			return;
		for(int cn(0); cn < all_char.size(); ++cn)
		{
			if( query[i+1] == all_char[cn])
				continue;
				
			query[i] = all_char[cn];
			exact_match_rec(result, query, i, pos_f, pos_e, mn+1);
			query[i] = c;
		}	
	}
	
	inline void exact_match_rec(std::vector<INTTYPE>& result, std::string &query, int i, INTTYPE pos_f, INTTYPE pos_e, INTTYPE b_pos_f, INTTYPE b_pos_e,	INTTYPE mn)
	{
		char c(query[i]);
		
		INTTYPE n_pos_f = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_f, c );
		INTTYPE n_pos_e = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( pos_e, c );
		
		if(i == 0)
		{
			//std::cerr << "ok print" << std::endl;
			push_result(result, n_pos_f, n_pos_e);
			return;
		}		
		if(pos_f < pos_e) //y
		{
			exact_match_rec(result, query, i-1, n_pos_f, n_pos_e, pos_f, pos_e, mn);

		}

	}
		
		
	inline void push_result(std::vector<INTTYPE>& result, INTTYPE pos_f, INTTYPE pos_e)
	{
		for (int i = pos_f; i < pos_e; i++)
		{
			result.push_back( find_nearest_mark(i) );
		}
	}
	

		
	inline void start_sbwt_match( std::string& query, std::vector<INTTYPE>& result )
	{
		//init_exact_match( query.back() );
		init_exact_match( query.substr(query.size()-c_functions_interval,c_functions_interval) );
		
		for (int i=query.length()-1-c_functions_interval; i>=0 && start_end_pos_.first < start_end_pos_.second; i--)
		//for (int i=query.length()-1-1; i>=0 && start_end_pos_.first < start_end_pos_.second; i--)
		{
			exec_exact_match(query[i]);	
		}
		//once there is any hit, this loop will run and print every hits
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
		
		//INTTYPE n_pos_f = abwt_table_.c_function[ c ] + abwt_table_.get_occ( pos_f, c );
		//INTTYPE n_pos_e = abwt_table_.c_function[ c ] + abwt_table_.get_occ( pos_e, c );
		
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

		//std::cerr << "len" << len << " size " << n_pos_e-n_pos_f << " c " << c << std::endl;



		for(int cn(0); cn < all_char.size(); ++cn)
		{
			find_possible(result, len+1, all_char[cn], n_pos_f, n_pos_e);
		}
	}


	inline void start_exact_match( std::string& query, std::vector<INTTYPE>& result )
	{
		//std::cerr << "start_exact_match" << std::endl;
		//initialize exact matching process
		init_exact_match( query.substr(query.size()-c_functions_interval,c_functions_interval) );
		//init_exact_match( query.back() );
		//std::cerr << "start_pos : " << start_end_pos_.first << " end_pos : " << start_end_pos_.second << std::endl;

		//extending matching query
		//for (int i=query.length()-1-1; i>=0 && start_end_pos_.first < start_end_pos_.second; i--)
		for (int i=query.length()-1-c_functions_interval; i>=0 && start_end_pos_.first < start_end_pos_.second; i--)
		{
			//c = query[i];
			//std::cerr << "start: "<< start_end_pos_.first << ":" << start_end_pos_.second << std::endl;
			exec_exact_match(query[i]);
		}
		//once there is any hit, this loop will run and print every hits
		if (start_end_pos_.second - start_end_pos_.first > 100)
			return;

		//std::cerr << "start: "<< start_end_pos_.first << ":" << start_end_pos_.second << std::endl;
		for (INTTYPE i = start_end_pos_.first; i < start_end_pos_.second; i++)
		{
			result.push_back( find_nearest_mark(i) );
		}
	}



	/////////////////////////////////////// Bo Han 	///////////////////////////////////////

	/*@ a version of exec_exact_match that return the previous value of start_end_pos_ */
	std::pair<INTTYPE, INTTYPE> exec_exact_match2 ( char c ) const {
		auto old_start_end_pos_ = start_end_pos_;
		start_end_pos_.first = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( start_end_pos_.first, c );
		start_end_pos_.second = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( start_end_pos_.second, c );
		return old_start_end_pos_;
	}

	/*@ using sense index, but since we RC the query sequence first, this is actually searching - strand */
	//TODO: result is a vector of Sam class
	void start_tailing_match_S (const Fastq& fq, std::ostream* out) const {
		std::string _query = fq.getSeq ();
		/* reverse complement query string */
		std::string query {_query.crbegin(), _query.crend()};
		for (auto & c : query) {
			switch (c) {
			case 'A': c = 'T'; break;
			case 'T': c = 'A'; break;
			case 'C': c = 'G'; break;
			case 'G': c = 'C'; break;
			default : throw "illegal char";
			}
		} /* end of RC */

		/* the front c_functions_interval has to be exact match */
		init_exact_match( query.substr (query.size() - c_functions_interval, c_functions_interval) );
		/* if not found, exit */
		if ( start_end_pos_.first >= start_end_pos_.second) {
			return;
		}
		/* queryPosition need to be tested again */
		int queryPosition = query.length() -1 - c_functions_interval;
		/* last step need to be recorded */
		auto last_start_end_pos_ = start_end_pos_;
		for (; queryPosition >= 0 && start_end_pos_.first < start_end_pos_.second; --queryPosition) {
			last_start_end_pos_ = this->exec_exact_match2(query[queryPosition]);
		}
/* recording tailing */
		if (start_end_pos_.first >= start_end_pos_.second) {
			++queryPosition; /// substract an extra one when exiting the loop, so add it back
			for (int i = last_start_end_pos_.first; i < last_start_end_pos_.second; i++) {
				auto position = this->find_nearest_mark(i);
				auto tailSeq = _query.substr(_query.size() - 1 - queryPosition);
				auto lowerIter = this->abwt_table_.chr_start_pos.upper_bound (position);
				std::advance (lowerIter, -1);
				auto chr = lowerIter->second;
				auto lowerIter3 = this->abwt_table_.chr_start_pos.upper_bound (position + _query.size () - 2 - queryPosition);
				std::advance (lowerIter3, -1);
				auto chr3 = lowerIter3->second;
				if (chr != chr3) continue;

				auto NLowerIter = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position);
				std::advance (NLowerIter, -1);
				auto NLowerIter3 = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position + _query.size () - 2 - queryPosition);
				std::advance (NLowerIter3, -1);
				if (NLowerIter != NLowerIter3) continue;
				position = position - lowerIter->first + NLowerIter->second;
				//TODO: redefine MAPQ
				*out << Sam {fq.getName (),
					Sam::SAM_FLAG::REVERSE_COMPLEMENTED,
					chr,
					position + queryPosition,
						255,
						std::to_string (_query.size() - 1 - queryPosition) + 'M' + std::to_string (queryPosition+1) + 'S',
						"*",
						0,
						0,
						_query,
						fq.getQuality (),
						std::unordered_map<std::string, std::string> {{"TL", tailSeq}}
				};
			}
			return;
		}
/* recording perfect match */
		if (queryPosition < 0) {
			for (int i = start_end_pos_.first; i < start_end_pos_.second; i++) {
				auto position = this->find_nearest_mark(i);

				auto lowerIter = this->abwt_table_.chr_start_pos.upper_bound (position);
				std::advance (lowerIter, -1);
				auto chr = lowerIter->second;
				auto lowerIter3 = this->abwt_table_.chr_start_pos.upper_bound (position + _query.size () -1 );
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
						chr,
						position+1,
						255,
						std::to_string (_query.size ()) + 'M',
						"*",
						0,
						0,
						_query,
						fq.getQuality ()
				};
			}
			return ;
		}
	}



	/*@ using antisense index, since we RC the query first, so this is actually + strand mapping */
	//TODO: result is a vector of Sam class
	void start_tailing_match_AS (const Fastq& fq, std::ostream* out) const {
		std::string _query = fq.getSeq ();
		/* reverse complement query string */
		std::string query {_query.crbegin(), _query.crend()};
		for (auto & c : query) {
			switch (c) {
			case 'A': c = 'T'; break;
			case 'T': c = 'A'; break;
			case 'C': c = 'G'; break;
			case 'G': c = 'C'; break;
			default : throw "illegal char";
			}
		} /* end of RC */
		/* the front c_functions_interval has to be exact match */
		init_exact_match( query.substr (query.size() - c_functions_interval, c_functions_interval) );
		/* if not found, exit */
		if ( start_end_pos_.first >= start_end_pos_.second) {
			return;
		}
		/* queryPosition need to be tested again */
		int queryPosition = query.length() -1 - c_functions_interval;
		/* last step need to be recorded */
		auto last_start_end_pos_ = start_end_pos_;
		for (; queryPosition >= 0 && start_end_pos_.first < start_end_pos_.second; --queryPosition) {
			last_start_end_pos_ = this->exec_exact_match2(query[queryPosition]);
		}
		if (start_end_pos_.first >= start_end_pos_.second) { /// begin recording tailing
			++queryPosition; /// substract an extra one when exiting the loop, so add it back
			for (int i = last_start_end_pos_.first; i < last_start_end_pos_.second; i++) {
				auto position = this->abwt_table_._realSize - this->find_nearest_mark(i) - _query.size () + queryPosition + 1;
				auto tailSeq = _query.substr(_query.size() - 1 - queryPosition);

				auto lowerIter = this->abwt_table_.chr_start_pos.upper_bound (position);
				std::advance (lowerIter, -1);
				auto chr = lowerIter->second;
				auto lowerIter3 = this->abwt_table_.chr_start_pos.upper_bound (position + _query.size () - 2 - queryPosition );
				std::advance (lowerIter3, -1);
				auto chr3 = lowerIter3->second;
				if (chr != chr3) continue;
				auto NLowerIter = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position);
				std::advance (NLowerIter, -1);
				auto NLowerIter3 = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position + _query.size () - 2 - queryPosition);
				std::advance (NLowerIter3, -1);
				if (NLowerIter != NLowerIter3) continue;
				position = position - lowerIter->first + NLowerIter->second;
				//TODO: redefine MAPQ
				*out << Sam {fq.getName (),
						Sam::SAM_FLAG::MAPPED,
						chr,
						position+1,
						255,
						std::to_string (_query.size() - 1 - queryPosition) + 'M' + std::to_string (queryPosition+1) + 'S',
						"*",
						0,
						0,
						_query,
						fq.getQuality (),
						std::unordered_map<std::string, std::string> {{"TL", tailSeq}}
				};
			}
			return;
		}
		if (queryPosition < 0) { /// found perfect match
			for (int i = start_end_pos_.first; i < start_end_pos_.second; i++) {

				auto position = this->abwt_table_._realSize - this->find_nearest_mark(i) - _query.size () + queryPosition + 1;
				auto lowerIter = this->abwt_table_.chr_start_pos.upper_bound (position);
				std::advance (lowerIter, -1);
				auto chr = lowerIter->second;
				auto lowerIter3 = this->abwt_table_.chr_start_pos.upper_bound (position + _query.size () - 1 );
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
					chr,
						position+1,
						255,
						std::to_string (_query.size ()) + 'M',
						"*",
						0,
						0,
						_query,
						fq.getQuality ()
				};
			}
			return ;
		}
	}

	// tailing searching version for dual BWT
	void start_tailing_match_Dual (const Fastq& fq, std::ostream* out, int minimalPrefixLength) const {
		std::string _query = fq.getSeq ();
		/* reverse complement query string */
		std::string query {_query.crbegin(), _query.crend()};
		for (auto & c : query) {
			switch (c) {
			case 'A': c = 'T'; break;
			case 'T': c = 'A'; break;
			case 'C': c = 'G'; break;
			case 'G': c = 'C'; break;
			default : throw "illegal char";
			}
		} /* end of RC */
		bool isRC = false;
		/* the front c_functions_interval has to be exact match */
		init_exact_match( query.substr (query.size() - c_functions_interval, c_functions_interval) );
		/* if not found, exit */
		if ( start_end_pos_.first >= start_end_pos_.second) {
			return;
		}
		/* queryPosition need to be tested again */
		int queryPosition = query.length() -1 - c_functions_interval;
		/* last step need to be recorded */
		auto last_start_end_pos_ = start_end_pos_;
		for (; queryPosition >= 0 && start_end_pos_.first < start_end_pos_.second; --queryPosition) {
			last_start_end_pos_ = this->exec_exact_match2(query[queryPosition]);
		}
/// begin recording tailing
		if (start_end_pos_.first >= start_end_pos_.second) {
			++queryPosition; /// substract an extra one when exiting the loop, so add it back
			for (int i = last_start_end_pos_.first; i < last_start_end_pos_.second; i++) {
				auto position = this->find_nearest_mark(i);
				auto prefixMatchLen = _query.size() - 1 - queryPosition;

				if (prefixMatchLen < minimalPrefixLength)
					continue;

				if (position > this->abwt_table_._realSize && position < (abwt_table_._realSize<<1)) { /// the second comparsion is to suppress weird bug of TTTTTTTTTTTT mapping to position == 2*abwt_table_._realSize
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
					*out << Sam { fq.getName (),
						Sam::SAM_FLAG::REVERSE_COMPLEMENTED,
						chr,
						position + queryPosition,
						255,
						std::to_string (prefixMatchLen) + 'M' + std::to_string (queryPosition+1) + 'S',
						"*",
						0,
						0,
						_query,
						fq.getQuality (),
						std::unordered_map<std::string, std::string> {{"TL", tailSeq}}
					};
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
					if (NLowerIter != NLowerIter3) continue;

					position = position - lowerIter->first + NLowerIter->second;
					//TODO: redefine MAPQ
					*out << Sam {fq.getName (),
						Sam::SAM_FLAG::MAPPED,
						chr,
						position+1,
						255,
						std::to_string (prefixMatchLen) + 'M' + std::to_string (queryPosition+1) + 'S',
						"*",
						0,
						0,
						_query,
						fq.getQuality (),
						std::unordered_map<std::string, std::string> {{"TL", tailSeq}}
					};
				}
			}
			return;
		}
/// found perfect match
		if (queryPosition == -1) {
			for (int i = start_end_pos_.first; i < start_end_pos_.second; i++) {
				auto position = this->find_nearest_mark(i);
				if (position > this->abwt_table_._realSize) {
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
						chr,
						position+1,
						255,
						std::to_string (_query.size ()) + 'M',
						"*",
						0,
						0,
						_query,
						fq.getQuality ()
					};
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
						chr,
						position+1,
						255,
						std::to_string (_query.size ()) + 'M',
						"*",
						0,
						0,
						_query,
						fq.getQuality ()
					};
				}
			}
			return ;
		}
	}
};


#endif
