#ifndef DIFFERENCE_COVER_HPP_
#define DIFFERENCE_COVER_HPP_
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

#include "compression/abit.hpp"
#include "mkq_sort.hpp"
#include "split_sort.hpp"
#include "lssort.hpp"
//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/archive/binary_oarchive.hpp>
//
//


//#define INTTYPE uint64_t

template <typename T>
void FreeAll( T & t ) {
		T tmp;
		t.swap( tmp );
}

class DIFFERENCE_COVER_IMPL
{
public:

	INTTYPE dcs_v, lv_size;
	std::vector<int> dcs_D;
	std::vector<int> dcs_Di;
	std::vector<int> dcs_dh;
	
	std::vector< INTTYPE > iD;
	std::vector< INTTYPE > iDs;
	std::vector< int > sp;
	std::vector< INTTYPE > ISAp;
	
	clock_t Tstart,start, stop;
	
	INTTYPE m_k, m_q;
	
	
	std::set<std::pair<INTTYPE,INTTYPE>> iD_same_rank;
	
	
	INTTYPE const_V, const_powerV;
	
	INTTYPE	limit_repeat_length;
	

	DIFFERENCE_COVER_IMPL ()
		:Tstart(clock())
	{}
	void def_dcs (int v)
	{
		if(v == 8)
		{
			const_V = 7;
			const_powerV = 3;
		}
		if(v == 16)
		{
			const_V = 15;
			const_powerV = 4;
		}
		if(v == 32)
		{
			const_V = 31;
			const_powerV = 5;
		}
		if(v == 64)
		{
			const_V = 63;
			const_powerV = 6;
		}
		if(v == 128)
		{
			const_V = 127;
			const_powerV = 7;
		}
		if(v == 256)
		{
			const_V = 255;
			const_powerV = 8;
		}
		if(v == 512)
		{
			const_V = 511;
			const_powerV = 9;
		}
		if(v == 1024)
		{
			const_V = 1023;
			const_powerV = 10;
		}
		
		std::vector< std::vector<int> > dcs(1025);
		
		dcs[7] = {1, 2, 4};
		dcs[8] = {1, 2, 3, 5};
		dcs[16] = {1, 2, 3, 6, 9};
		dcs[32] = {1, 2, 3, 4, 8, 12, 20};
		dcs[64] = {1, 2, 3, 6, 15, 17, 35, 43, 60};
		dcs[128] = {1,2,4,8,18,41,56,65,76,86,105,110,118};
		dcs[256] = {1,2,4,8,13,21,31,45,66,81,90,97,115,123,129,151,197,198,202,220};
		dcs[512] = {1,2,3,4,5,10,19,28,37,46,65,84,103,122,141,160,179,198,217,227,237,247,257,267,268,269,270,271};
		dcs[1024] = {1,2,3,4,5,6,7,14,27,40,53,66,79,92,119,146,173,200,227,254,281,308,335,362,389,416,443,457,471,485,499,513,527,541,542,543,544,545,546,547};		
		
		dcs_v = v;
		dcs_D = dcs[v];
		make_dcs_Di_table();
		make_dcs_dh_table();
		
	}
	inline void make_dcs_dh_table()
	{
		//#include "Dh_1048576.hpp";
		//#include "Dh_1024.hpp";
		bool ck=false;
		for(INTTYPE i(0); i<dcs_v; ++i)
		{
			ck=false;
			for(INTTYPE j(0); j<dcs_D.size();++j)
			{
				if(ck==true) break;
				for(INTTYPE k(0); k<dcs_D.size();++k)
				{
					if(mod_v( dcs_D[k] - dcs_D[j] ) == i)
					{
						dcs_dh.push_back(dcs_D[j]);
						ck=true;
						break;
					}
				}
			}
		}
	}
	inline void make_dcs_Di_table()
	{
		dcs_Di.resize( dcs_D[dcs_D.size()-1]+1, 0);
		for(int i(0),j(0); i<=dcs_D[dcs_D.size()-1]; ++i )
		{
			if(dcs_D[j]==i)
			{
				dcs_Di[i] = j;
				++j;
			}		
		}
	}
	inline bool compare(INTTYPE i, INTTYPE j)	{return compare_dcs_level_from_location(i,j);}
	inline INTTYPE i2ui(INTTYPE a)
	{
		INTTYPE mod(dcs_Di[ (a & const_V) ]);
		return	(a >> const_powerV) + mod*m_k + (m_q<mod?m_q:mod); //std::min(m_q, mod);
	}
	inline INTTYPE mod_v ( INTTYPE value )
	{
			return value & const_V;
	}
	inline INTTYPE mod_v ( INTTYPE value, INTTYPE value2 ) //value - value2
	{
		if (value < value2)
		{
			INTTYPE temp = ( -( (value2-value) & const_V) + const_V+1 );
			if (temp == const_V+1)
				return 0;
			else
				return temp;		
		}
		else
			return (value - value2) & const_V;
	}
	inline INTTYPE count_lv_size(INTTYPE seq_size)
	{
		INTTYPE lv_size = (seq_size / dcs_v)*dcs_D.size();
		for( INTTYPE i : dcs_D)
		{
			if(seq_size > lv_size/dcs_D.size()*dcs_v+i) 
				++lv_size;
		}
		return lv_size;	
	}
	
	inline bool compare_dcs_level_from_location	(INTTYPE i, INTTYPE j)
	{
		INTTYPE delta =	mod_v( dcs_dh[ mod_v( j, i ) ] , i ) ;
		return iDs[changeIndex2iD( i + delta )]	<	iDs[changeIndex2iD( j + delta )] ;
	}
	
	inline INTTYPE changeIndex2iD(INTTYPE i)
	{
		return dcs_Di[ (i & const_V) ]	+ (i >> const_powerV) * dcs_D.size();
	}
	
	inline bool compare_iD_rank_from_sp(INTTYPE a, INTTYPE b)
	{
		INTTYPE aa, bb;
		while(a < sp.size() && b < sp.size())
		{
			aa=sp[a];
			bb=sp[b];
			
			if(aa > bb )
				return false;
			else if(aa < bb)
				return true;
			++a;
			++b;
		}
		
		if(a>b) return true;
		else return false;
	}	


//old 1
	inline void make_dcs_iD_table(INTTYPE seq_size)
	{
		lv_size = count_lv_size(seq_size);
		
		start = clock();
		
		iD.reserve( lv_size );
		
		INTTYPE ui(0);
		for(size_t i=0; i< lv_size; ++i)
		{
			size_t mod = i % dcs_D.size();
			size_t except = i / dcs_D.size();
			//ui = except + mod*(lv_size/ dcs_D.size()) + std::min(lv_size%dcs_D.size(), mod) ;
			iD.push_back( dcs_D[mod]+except*dcs_v );
		}
		
		m_k = lv_size / dcs_D.size();
		m_q = lv_size % dcs_D.size();
	}


};





template<typename SEQTYPE, typename SORTTYPE>
class DIFFERENCE_COVER
	: public DIFFERENCE_COVER_IMPL
{
public:
	SEQTYPE &seq;
	INTTYPE len_seq, dcs_v;
	INTTYPE tmp_v, tmp_min, tmp_max;
	std::vector<std::pair<INTTYPE,INTTYPE>> iD_sames;
	
	DIFFERENCE_COVER (SEQTYPE &sequence, INTTYPE dcs_size=512)
		:seq(sequence), len_seq(sequence.size()), tmp_v(0), dcs_v(dcs_size)
	{
		std::cerr << "dcs : " << dcs_v << std::endl;
		def_dcs(dcs_v);
		make_dcs_table_mkq();
	}
	
	void make_dcs_table_mkq()
	{
//old 1
		make_dcs_iD_table(len_seq);
//new 2
		sort_pre_dcs_size();
//new 3
		INTTYPE rank = make_sp_table();
//new 4
		ls_sort(rank);
//new 5
		prepare_iDs();		
	}
	
//new2
	void sort_pre_dcs_size()
	{
		start = clock();
		
		INTTYPE limit_sort_length(dcs_v);
		SORTTYPE pre_sorter(seq, iD, limit_sort_length, iD_sames);
		
		// sort(start, length, depth)
		pre_sorter.sort(0,iD.size(),0);
		stop = clock();
		std::clog << "Sort iD table End, time:" << double(stop - start) / CLOCKS_PER_SEC << "\n" << std::endl;
	}

//new3	
	INTTYPE make_sp_table()
	{
		//make sp table
		INTTYPE i(0), j(0), begin(0), end(0), state(0) ,rank(0);
		sp.resize(iD.size(),0);
		if(iD_sames.size() > 0)
		{
			begin = iD_sames[0].first;
			end = iD_sames[0].second;
		
			for( INTTYPE q = 0; q<iD.size(); ++q )
			{
				sp[ i2ui( iD[q] ) ] = rank;
				
				if(i == begin)
					state = 1;
				
				if(i == end)
				{
					++j;
					begin = iD_sames[j].first;
					end = iD_sames[j].second;
					state = 0;
				}
				if(state !=1 )
					++rank;
				++i;
			}
		}else
		{
			for( INTTYPE q = 0; q<iD.size(); ++q )
			{
				sp[ i2ui( iD[q] ) ] = rank;
			}
			
		}
		//release iD
		FreeAll(iD);	
		
		return rank;
	}

//new 4	
	void ls_sort(INTTYPE rank)
	{
		start = clock();
		std::vector<int> sorted_sp;
		sp.push_back(-1);
		sorted_sp.reserve(sp.size());
		for (int i=0;i<sp.size(); ++i)
			sorted_sp.push_back(i);
		
		suffixsort( sp.data(), sorted_sp.data(), sp.size()-1, rank+1, 0 );
		stop = clock();
		
		//release sorted_sp
		FreeAll(sorted_sp);
		std::clog << "reSort iD table, fix sp rank, using sp table, End, time:" << double(stop - start) / CLOCKS_PER_SEC << "\n" << std::endl;
	}

//new 5	
	void prepare_iDs()
	{
		size_t mk(lv_size/ dcs_D.size());
		size_t mg(lv_size%dcs_D.size());
		iDs.resize(sp.size(),0);
		for(INTTYPE i(0); i < sp.size(); ++i)
		{
			size_t mod = i % dcs_D.size();
			size_t except = i / dcs_D.size();
			INTTYPE ui = except + mod * mk + std::min(mg, mod);
			iDs[i] = sp[ui];	
		}
		//release sp
		FreeAll(sp);
	}
	
};






/*
template<class SORTTYPE>
class DIFFERENCE_COVER< ABSequence<ABmpl>, SORTTYPE>
	: public DIFFERENCE_COVER_IMPL
{
public:
	ABSequence<ABmpl> &seq;
	INTTYPE len_seq;
	
	DIFFERENCE_COVER (ABSequence<ABmpl> &sequence, INTTYPE dcs_v=512)
		:seq(sequence), len_seq(sequence.size())
	{
		std::cerr << "dcs : " << dcs_v << std::endl;
		def_dcs(dcs_v);
		make_dcs_table_default(len_seq,
			[&](){
				QSeqTable < INTTYPE, std::tuple<uint64_t, char, INTTYPE> > quick_table_A(seq.EscapeChar);
				QSeqTable < INTTYPE, std::tuple<uint64_t, char, INTTYPE> > quick_table_B(seq.EscapeChar);
				std::sort (iD.begin(), iD.end(), 
					[&] (const INTTYPE &A, const INTTYPE &B)
					{
						INTTYPE a(A), b(B);

						char result = seq.compare(a,b,dcs_v,quick_table_A,quick_table_B);
						if(result=='>') return false;
						else if(result=='<') return true;
						iD_same_rank.insert( typename decltype(iD_same_rank)::value_type (std::minmax(a,b) ));
						//iD_same_rank.insert({std::min(a,b), std::max(a,b) });
						return false;
					}
				);
			}
		);
		
		
	}
	
	
	template<typename SORTFUNTYPE>
	void make_dcs_table_default(INTTYPE seq_size, SORTFUNTYPE sort_fun_with_length)
	{
//old 1
		make_dcs_iD_table(seq_size);
//old 2
		sort_dcs_iD_table_using_iD(sort_fun_with_length);
//old 3
		insert_sp2iD_table();
//old 4
		make_sp_table();
//old 5
		reSort_iD_table_with_sp();
//old 6
		fix_sp_rank();
	}
	
	
//old 2	
	template<typename SORTFUNTYPE>
	inline void sort_dcs_iD_table_using_iD(SORTFUNTYPE sort_fun_with_length)
	{
		std::clog << "Sort dcs iD table, using iD as key, Start" << std::endl;
		start = clock();
		
		sort_fun_with_length();
		
		std::clog << std::endl;
		
		
		stop = clock();
		std::clog << "Sort dcs iD table, using iD as key, End, time:" << double(stop - start) / CLOCKS_PER_SEC << "\n" << std::endl;
		
	}
	
//old 3
	inline void insert_sp2iD_table()
	{
		std::clog << "Insert sp to dcs iD table, sp, Start" << std::endl;
		start = clock();
		
		//insert sp
		sp.resize(iD.size(),0);
		
		INTTYPE rank(0), pre_value(0), now_value(0), i(0);
		for( INTTYPE &K : iD)
		{
			now_value = K;
			
			if( iD_same_rank.find( {std::min(pre_value, now_value), std::max(pre_value, now_value)} ) == iD_same_rank.end() )
			{
				++rank;
			}
			if( i2ui(K) >= sp.size() )
				std::cerr << K << ":" << i2ui(K) << ":" << changeIndex2iD(K) << std::endl;
			sp[ i2ui(K) ] = rank -1;
			pre_value = now_value;
			++i;
		}
		std::clog << std::endl;
		stop = clock();
		std::clog << "Insert sp to dcs iD table, sp, End, time:" << double(stop - start) / CLOCKS_PER_SEC << "\n" << std::endl;
		
	}

//old 4
	void make_sp_table()
	{}

//old 5
	void reSort_iD_table_with_sp()
	{
		std::clog << "reSort iD table, fix sp rank, using sp table, Start" << std::endl;
		start = clock();
		
		//re sort iD-> sp
		std::sort (iD.begin(), iD.end(), 
			[&] (const INTTYPE &a, const INTTYPE &b)
			{
				return compare_iD_rank_from_sp( a , b );
			});
		stop = clock();
		std::clog << "reSort iD table, fix sp rank, using sp table, End, time:" << double(stop - start) / CLOCKS_PER_SEC << "\n" << std::endl;
	}

//old 6
	void fix_sp_rank()
	{
		std::clog << "fix sp rank in iD table, Start" << std::endl;
		start = clock();
		
		iDs.resize(iD.size(),0);
		//fix sp
		INTTYPE index(0), pre_rnak(0), now_rank(0);
		for( INTTYPE &K : iD)
		{
			iDs[ changeIndex2iD(K) ] = index;
			++index;
		}
		stop = clock();
		std::clog << "fix sp rank in iD table, End, time:" << double(stop - start) / CLOCKS_PER_SEC << "\n" << std::endl;
	}
	

};

*/


#endif
