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
#ifndef ABWT_HPP_
#define ABWT_HPP_
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
#include "mkq_sort.hpp"
#include "abwt_table.hpp"
#include "abwt_search.hpp"

#include "boost/serialization/vector.hpp"
#include "boost/serialization/utility.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/unordered_map.hpp"
//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/archive/binary_oarchive.hpp>

/**
 * @struct ABWT<>
 * @brief 全特化樣板，無功能
 * @tparam SeqType 通用版本目前無實作，只用來做全特化
 *
 */

//#define INTTYPE uint64_t


template< class SeqType>
class ABWT
{};


/**
 * @brief 正常版本BWT演算法實作，使用 ABSequence<std::string>
 * @tparam SeqType=ABSequence<std::string> > 正常版本為std::string模式，不壓縮
 *
 */
template<>
class ABWT< ABSequence<> >
{

public:

	/// @brief 定義 SeqType 為 ABSequence
	typedef ABSequence<> SeqType;

	/// @brief 定義 TableType 為儲存bwt location (int)的容器
	typedef std::vector< INTTYPE > TableType;

	/// @brief 定義在做dcs表時，sort前dcs length個字所用的 sort type \n
	/// 這邊使用 multikey quick sort ，並且要開啓紀錄相等模式
	typedef Multikey_quicksort
					<	SeqType,
						TableType,
						Record_rank_enable,
						Compare_default,
						Sort_small_n_disable,
						Bucket_sort
					>	PreSortType;

	/// @brief 定義主體sort suffix 的 sort type \n
	/// 這邊使用 multikey quick sort ，並且使用 dcs compare
	typedef Multikey_quicksort
					<	SeqType,
						TableType,
						Record_rank_disable,
						Compare_DCS,
						Sort_small_n_enable,
						Bucket_sort
					>	MainSortType;

	/// @brief 原始序列儲存
	SeqType &seq;

	/// @brief 排序 All Suffix後，儲存每個Suffix的第一個字在原始sequence的位置(uint32_t) \n
	/// 此時 int 對應到 seq 的字元應該是： $AAAA...CC...G...TT... \n
	/// 而前一個字元(int -1)對應到的字元則為 BWT \n
	TableType seq_table;

	/// @brief ABWT_table 為處理seq_table的class，主要負責產生所有在 search所需要的表 \n
	/// @brief 像是： bwt(char), o table, occ table, location table 等等，和一些加速的 table
	ABWT_table abwtt;

	/// @brief 紀錄原始sequence的總長度
	INTTYPE len_seq;

	/// @brief difference cover 的參數，為2的倍數，關係到dcs的大小，速度
	INTTYPE dcs_v;

	/// @brief 排序時，分割每群的大小 \n
	/// @brief 在排序時，是分成好幾群，依序排序，可以有效減少記憶體用量，也可以多cpu執行
	INTTYPE split_size;

	/// @brief 記錄時間
	clock_t start, stop;

	/// @brief prefix name for bwt file
	std::string _prefixName;

//	TODO: should I put these information in ABWT or in a separate table?
//	/// @brief record the offsets of Ns
//	std::set <uint32_t> _validRegions;
//	/// @brief record the region of chromosomes
//	std::map <uint32_t, std::string> _chrSizes;

	/**
	 * @memberof class ABWT< ABSequence<std::string> >
	 * @brief 建構式與程式主體
	 * @param sequence ABSequence<std::string> 原始待排序序列
	 * @param dcs_size uint32_t dcs 參數，攸關dcs大小速度
	 * @param interval_size uint32_t 在建abwt表時所使用的間隔參數，越小越快，記憶體使用越大
	 */
	ABWT (SeqType &sequence, uint32_t dcs_size=512, uint32_t interval_size=64, std::string prefixName="")
		:seq(sequence), len_seq(seq.size()) , dcs_v(dcs_size), abwtt(interval_size), split_size(30000000), _prefixName {prefixName}

	{

		//part 1
		/// @brief 測試
		start = clock();
		// make dcs table
		DIFFERENCE_COVER < SeqType, PreSortType > dcs(seq, dcs_v);

		std::function<bool(INTTYPE, INTTYPE)> dcs_compare_func = std::bind( &DIFFERENCE_COVER < SeqType, PreSortType >::compare, dcs, std::placeholders::_1, std::placeholders::_2);

		stop = clock();
		std::clog << "===================================" << std::endl;
		std::clog << "Create dcs table, time:" << double(stop - start) / CLOCKS_PER_SEC << std::endl;
		std::clog << "===================================" << std::endl;


		//part 2

		start = clock();

		// sort suffix using split sort
		std::cout << "split_size : " << split_size << std::endl;
		//new split sort obj


		Split_sort<	SeqType, std::vector< INTTYPE >, MainSortType > split_sort (seq,dcs_v);
		//setting
		split_sort.getTableV = [](INTTYPE& a){
			return a;
		};
		split_sort.dcs_compare = dcs_compare_func;
		//start split sort
		split_sort.split_by_size_init(split_size);

		// copy sorted result to seq_table(bwt)
		//seq_table.reserve( seq.size() );
		//seq_table.swap(split_sort.Results);
		std::vector<std::string> filenames = split_sort.archive_name;

		//release dcs memory
		std::cout << "release dcs..." << std::endl;
		dcs.release();
		std::cout << "release split sort..." << std::endl;
		split_sort.release();


		stop = clock();
		std::clog << "===================================" << std::endl;
		std::clog << "Create BWT, time:" << double(stop - start) / CLOCKS_PER_SEC << std::endl;
		std::clog << "===================================" << std::endl;


		//part 3

		start = clock();
		// save abwt table

		abwtt.createAllTable(seq, filenames);
		//abwtt.saveBWT("t_bwt.bwt");
		//abwtt.saveTable("t_table.bwt");
		//abwtt.saveSEQ("t_seq.bwt", seq);
		abwtt.saveBWT( _prefixName + "t_bwt.bwt");
		abwtt.saveTable( _prefixName + "t_table.bwt");
		abwtt.saveSEQ( _prefixName + "t_seq.bwt", seq);


		stop = clock();
		std::clog << "===================================" << std::endl;
		std::clog << "Create abwt Table, time:" << double(stop - start) / CLOCKS_PER_SEC << std::endl;
		std::clog << "===================================" << std::endl;

		split_sort.clean_up ();

		return ;

	}


};

#endif
