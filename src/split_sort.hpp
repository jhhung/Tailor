#ifndef SPLIT_SORT_HPP_
#define SPLIT_SORT_HPP_
#include <fstream>
#include <map>
#include <list>
#include <random>
#include <algorithm>
#include <functional>
#include <boost/ref.hpp>
#include <utility>
#include <cstdlib>
#include <queue>
#include <tuple>
#include <ctime>
#include "mkq_sort.hpp"
#include "compression/abit.hpp"

#include "boost/serialization/vector.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/filesystem.hpp"
#include <memory>

//#define INTTYPE uint64_t

template <class SeqType, class TableType, class SortType>
class Split_sort
{
private:

	SeqType &seq;
	TableType &vec, self_vec;
	
	INTTYPE average_size, split_num, random_num, compare_length, ref_size, ori_table_size;
	std::vector<INTTYPE> random_table, split_table;
	
	std::vector< TableType > SeqTables;
	
	
	INTTYPE limit_depth;

public:
	TableType Results;
	std::function<INTTYPE(typename TableType::value_type&)> getTableV;
	std::function<bool(INTTYPE, INTTYPE)> dcs_compare;
	
	
	std::vector<std::string> archive_name;
	std::vector<INTTYPE> object_in_each_archive;
	std::vector<INTTYPE> size_in_each_archive;
	std::vector< std::shared_ptr<std::ofstream> > ofile_list;
	std::vector< std::shared_ptr<std::ifstream> > ifile_list;
	std::vector< std::shared_ptr< boost::archive::binary_oarchive> > oarchive_list;
	std::vector< std::shared_ptr< boost::archive::binary_iarchive> > iarchive_list;
	
	

	Split_sort(SeqType& sequence, TableType& table, INTTYPE ld)
		:seq(sequence), vec(table), limit_depth(ld)
	{
		
	}
	Split_sort(SeqType& sequence, INTTYPE ld)
		:seq(sequence), vec(self_vec), limit_depth(ld)
	{
		
	}
	
	void split_by_size_init(const INTTYPE s)
	{
		clock_t start, stop;
		average_size = s;
		split_num = seq.size() / average_size;
		if(split_num == 0) split_num = 1;
		random_num = split_num * 4;
		SeqTables.resize(split_num, std::vector<typename TableType::value_type> () );
		
		make_random_table();
		sort_random_table();
		make_split_table();
		
		start = clock();
		std::clog << "classify start:" <<std::endl;
		classify_seq_tables();
		stop = clock();
		std::clog << "classify end: "<< double(stop - start)/CLOCKS_PER_SEC << std::endl;
		
		INTTYPE i(0);
		for(auto &K : SeqTables)
		{
			std::clog << "group: "<< i << ":" << K.size() << std::endl;
			++i;
		}
		
		start = clock();
		std::clog << "sort final start:" <<std::endl;
		
		//34310
		//archive_name.push_back("split_0");
		//archive_name.push_back("split_1");
		//object_in_each_archive.push_back(34310);
		//object_in_each_archive.push_back(34000);
		
		mkq_sort();
		
		stop = clock();
		std::clog << "sort final end: "<< double(stop - start)/CLOCKS_PER_SEC << std::endl;

		merge();
	}
	
	inline void make_random_table()
	{
		for(INTTYPE i(0); i < random_num; ++i)
		{
			INTTYPE val = std::rand() % seq.size();
			random_table.push_back(val);
		}
		std::sort(random_table.begin(),random_table.end());
		std::unique(random_table.begin(), random_table.end());
	}
	
	inline void sort_random_table()
	{
		
		std::sort(random_table.begin(), random_table.end(),
			[&]( typename TableType::value_type A, typename TableType::value_type B){
				INTTYPE a(getTableV(A)), b(getTableV(B));
				for(INTTYPE i(0); i <= limit_depth; ++i)
				{	
					if( seq[(a+i)] > seq[(b+i)] )
						return false;
					else if( seq[(a+i)] < seq[(b+i)] )
						return true;
				}
				return dcs_compare(a,b);
			}
		);
	}
	
	
	
	inline void make_split_table()
	{
		
		for(INTTYPE i(1); i < (split_num); ++i)
		{
			INTTYPE idx = random_table.size() / (split_num) * i;
			split_table.push_back(random_table[idx]);
		
		}
		for (INTTYPE i(0); i < (split_num); ++i)
		{
				archive_name.push_back( std::string("split_")+boost::lexical_cast<std::string>(i) );
				ofile_list.push_back( std::shared_ptr<std::ofstream> ( new std::ofstream(archive_name[i], std::ios::binary) ) );
				oarchive_list.push_back( std::shared_ptr<boost::archive::binary_oarchive> ( new boost::archive::binary_oarchive( *ofile_list[i] ) ) );
				object_in_each_archive.push_back(0);
				size_in_each_archive.push_back(0);
		}
		std::cerr << "table size : " <<	SeqTables.size() << std::endl;
	}

	void clean_up () const {
		for(INTTYPE i(0); i < (split_num); ++i) {
			boost::filesystem::remove_all (std::string("split_") + std::to_string(i));
		}
	}

	inline void classify_seq_tables()
	{	
		clock_t start, stop;
		std::vector<INTTYPE>::iterator split_table_it(split_table.begin());
		INTTYPE idx(0);
		
		for(INTTYPE K(0); K<seq.size(); ++K)
		{
			split_table_it = std::lower_bound(split_table.begin(), split_table.end(), K ,
				[&](INTTYPE a, INTTYPE b){
					for(INTTYPE i(0); i <= limit_depth; ++i)
					{	
						if( seq[(a+i)] > seq[(b+i)] )
							return false;
						else if( seq[(a+i)] < seq[(b+i)] )
							return true;
					}
					return dcs_compare(a,b);
				}
			);
			idx = split_table_it - split_table.begin() ;
			SeqTables[idx].push_back(K);
			
			if (SeqTables[idx].size()>1000)
			{
					size_in_each_archive[idx] += SeqTables[idx].size();
					*oarchive_list[idx] & SeqTables[idx];
					object_in_each_archive[idx]++;
					SeqTables[idx].clear();
			};
			
		}
		
		idx=0;
		for(auto &KK : SeqTables)
		{
			
			size_in_each_archive[idx] += KK.size();
			*oarchive_list[idx] & KK;
			ofile_list[idx]->close();
			std::cerr<<"done outputing: "<<archive_name[idx]<< "size: "<<size_in_each_archive[idx]<<std::endl;
			object_in_each_archive[idx]++;
			KK.clear();
			idx++;
		}
		
		
	}
	
	inline void mkq_sort()
	{
		clock_t start, stop;
		INTTYPE i(0);
		for(TableType &K : SeqTables)
		{
			start = clock();
			std::ifstream in( archive_name[i], std::ios::binary );
			std::cerr<<"reading: "<<archive_name[i] << " " << object_in_each_archive[i] <<std::endl;
			boost::archive::binary_iarchive temp_archive(in);
			for (INTTYPE j=0; j<object_in_each_archive[i]; ++j)
			{
				TableType temp_vector;
				temp_archive & temp_vector;
				K.insert( K.end(), temp_vector.begin(), temp_vector.end()	); 
			}
			in.close();
			
			std::clog << "Sort start:" << i << " size : " << K.size() <<std::endl;
			
			SortType sorter(seq, K, limit_depth, dcs_compare);
																							
			sorter.getTableV = [](INTTYPE& a){
					return a;
			};
			sorter.sort(0,K.size(),0);
			
			stop = clock();

			std::cerr << "Sort end: "<< double(stop - start)/CLOCKS_PER_SEC << std::endl;
			
			std::ofstream out( archive_name[i], std::ios::binary );
			boost::archive::binary_oarchive temp_archive_o(out);
			temp_archive_o & K;
			out.close();
			
			
			TableType temp;
			K.swap(temp);
			
			std::clog << "Sort end: "<< double(stop - start)/CLOCKS_PER_SEC << std::endl;
			++i;
		}
	}
	inline void merge()
	{
		Results.resize(seq.size());
		INTTYPE start(0), end(0);
		for(TableType &K : SeqTables)
		{
			end += K.size();
			std::swap_ranges(Results.begin()+start, Results.begin()+end, K.begin());
			start = end;
		}
	}
};


#endif
