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
#include <boost/ref.hpp>
#include <ctime>
#include <cmath>
#include <memory>
#include "compression/abit.hpp"
#include "compression/jbit.hpp"
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

template< class SeqType>
class ABWT
{};

template<>
class ABWT< ABSequence<> >
{

public:
	typedef ABSequence<> SeqType;
	typedef std::vector< INTTYPE > TableType;
	typedef Multikey_quicksort
					<	SeqType,
						TableType,
						Record_rank_enable,
						Compare_default,
						Sort_small_n_disable,
						Bucket_sort
					>	PreSortType;
	typedef Multikey_quicksort
					<	SeqType,
						TableType,
						Record_rank_disable,
						Compare_DCS,
						Sort_small_n_enable,
						Bucket_sort
					>	MainSortType;
	SeqType &seq;
	TableType seq_table;
	ABWT_table abwtt;
	INTTYPE len_seq; 
	INTTYPE dcs_v; 
	INTTYPE split_size;
	clock_t start, stop;
	std::string _prefixName;
	ABWT (SeqType &sequence, uint32_t dcs_size=512, uint32_t interval_size=64, std::string prefixName="")
		:seq(sequence), len_seq(seq.size()) , dcs_v(dcs_size), abwtt(interval_size), split_size(30000000), _prefixName {prefixName}
	{		
		start = clock();
		// make dcs table
		DIFFERENCE_COVER < SeqType, PreSortType > dcs(seq, dcs_v);
		std::function<bool(INTTYPE, INTTYPE)> dcs_compare_func = std::bind( &DIFFERENCE_COVER < SeqType, PreSortType >::compare, dcs, std::placeholders::_1, std::placeholders::_2);	
		stop = clock();
		std::clog << "===================================" << std::endl;
		std::clog << "Create dcs table, time:" << double(stop - start) / CLOCKS_PER_SEC << std::endl;
		std::clog << "===================================" << std::endl;
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
		stop = clock();
		std::clog << "===================================" << std::endl;
		std::clog << "Create BWT, time:" << double(stop - start) / CLOCKS_PER_SEC << std::endl;
		std::clog << "===================================" << std::endl;
		start = clock();
		abwtt.createAllTable(seq, filenames);
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
