/// @file n_bit_compress.hpp
/// @brief define a data structure for achieving N bit encode/compression scheme and accordingly saving the encoded result.  Overall, two kinds of classes, i.e. CompressRuleTrait and NBitCompress are provided. \n The CompressRuleTrait is served as a template parameter for the encode/compression core class of NBitCompress, wherein the CompressRuleTrait ACTUALLY conveys the information of how many bits and the actual encode/compression rule that the encode/compression sceme takes. \n The NBitCompress takes the responsibility to carry out the encode/compression operation. 
/// @author JHH Corp.
#ifndef COMPRESSRULE_TRAIT_HPP_
#define COMPRESSRULE_TRAIT_HPP_
#include <iostream>
#include "boost/dynamic_bitset.hpp"
#include <unordered_map>
#include <string>
#include <cstdint>
#include <deque>
#include "../constant_def.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"

/// @struct CompressRuleTrait
/// @brief served as a encode/compression scheme selector by conveying the information of how many bits and the actual encode/compression rule that the encode/compression sceme takes.
/// @tparam N how many bits that the encode/compression sceme employes.  A default value of 2 is designed for the non-type parameter of N.
/// @tparam CompressType a enumerated index indicating what encode/compress scheme is applied.  \n For N_BITS situation, the NBitCompress class takes some characters as target characters, such as characters 'A', 'C', 'G', 'T', for which encode/conpression rules are provided in the TraiClass, and the rest of the characters as except characters, for which no encode/compression rules are provided. 
template <int N, CompressType C>
struct CompressRuleTrait
{
	template <int, CompressType>
	friend bool operator== (const CompressRuleTrait<N, C>& q, const CompressRuleTrait<N, C>& qq);
};

/// @brief a specialization version of the CompressRuleTrait with its template parameter CompressType specialized to be CompressType::N_BITS.  \n For N_BITS situation with N=2, input characters are segmented into two different groups, target and except. Characters of the target group is encoded with two-bit scheme, while the except group is handled as except.  Thus, two different data processes are employed to respectively taking care of the target and the except groups.
template <>
struct CompressRuleTrait <2, CompressType::N_BITS >
{
	static const int type = CompressType::N_BITS;
	static const int bit_count = 2;
	std::string alphabet;
	std::unordered_map < char, std::pair< bool, boost::dynamic_bitset<> > > table;
	CompressRuleTrait ()
		: table ({
	   		{'A', std::make_pair (true, boost::dynamic_bitset<> (2, 0ul))},
	   		{'a', std::make_pair (false,boost::dynamic_bitset<> (2, 0ul))},
 			{'C', std::make_pair (true, boost::dynamic_bitset<> (2, 1ul))},
			{'c', std::make_pair (false,boost::dynamic_bitset<> (2, 1ul))},
		   	{'G', std::make_pair (true, boost::dynamic_bitset<> (2, 2ul))},
		   	{'g', std::make_pair (false,boost::dynamic_bitset<> (2, 2ul))},
			{'T', std::make_pair (true, boost::dynamic_bitset<> (2, 3ul))},
			{'t', std::make_pair (false,boost::dynamic_bitset<> (2, 3ul))},
			{'U', std::make_pair (false,boost::dynamic_bitset<> (2, 3ul))},
			{'u', std::make_pair (false,boost::dynamic_bitset<> (2, 3ul))}
				})
	{
		std::for_each ( table.begin(), table.end(),
						[&] (const std::pair < char, std::pair<bool,boost::dynamic_bitset<> > > Q)
						{ alphabet+= Q.first;}
					);
	}

	friend class boost::serialization::access;
	template <typename Archive>
	void serialize (Archive& ar, const unsigned int)
	{
		ar & alphabet;
		std::for_each (table.begin(), table.end(),
						[&] (std::pair <char, std::pair<bool, boost::dynamic_bitset<> > > Q)
						{ 
							std::vector<boost::dynamic_bitset<>::block_type, boost::dynamic_bitset<>::allocator_type>* q;
							to_block_range(Q.second.second, &q);
							size_t* qq; 
							to_block_range(Q.second.second, &qq);
							ar & Q.first & Q.second.first & (*q) & (*qq);
						}
						);
	}
};

/// @brief a specialization version of the CompressRuleTrait with its template parameter CompressType specialized to be CompressType::N_BITS_MASK.  \n For N_BITS_MASK situation with N=2, an extra mask rule, such as lowercase character, is employed by the NBitCompress, so tht three different data processes are employed to respectively taking care of the target, the except, and the mask characters.
template <>
struct CompressRuleTrait <2, CompressType::N_BITS_MASK >
{
	static const int type = CompressType::N_BITS_MASK;
	static const int bit_count = 2;
	std::string alphabet;
	std::unordered_map < char, std::pair< bool, boost::dynamic_bitset<> > > table;
	CompressRuleTrait ()
		: table ({
	   		{'A', std::make_pair (true, boost::dynamic_bitset<> (2, 0ul))},
	   		{'a', std::make_pair (false,boost::dynamic_bitset<> (2, 0ul))},
 			{'C', std::make_pair (true, boost::dynamic_bitset<> (2, 1ul))},
			{'c', std::make_pair (false,boost::dynamic_bitset<> (2, 1ul))},
		   	{'G', std::make_pair (true, boost::dynamic_bitset<> (2, 2ul))},
		   	{'g', std::make_pair (false,boost::dynamic_bitset<> (2, 2ul))},
			{'T', std::make_pair (true, boost::dynamic_bitset<> (2, 3ul))},
			{'t', std::make_pair (false,boost::dynamic_bitset<> (2, 3ul))},
			{'U', std::make_pair (false,boost::dynamic_bitset<> (2, 3ul))},
			{'u', std::make_pair (false,boost::dynamic_bitset<> (2, 3ul))}
//	   		{'A', std::make_pair (true, boost::dynamic_bitset<> (2, 0ul))},
 //			{'C', std::make_pair (true, boost::dynamic_bitset<> (2, 1ul))},
//		   	{'G', std::make_pair (true, boost::dynamic_bitset<> (2, 2ul))},
//			{'T', std::make_pair (true, boost::dynamic_bitset<> (2, 3ul))},
//			{'U', std::make_pair (false, boost::dynamic_bitset<> (2, 3ul))}
				})
	{
		std::for_each ( table.begin(), table.end(),
						[&] (const std::pair < char, std::pair<bool, boost::dynamic_bitset<> > > Q)
						{ alphabet += Q.first;}
						);
	}
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & alphabet;
		std::for_each (table.begin(), table.end(),
						[&] (std::pair <char, std::pair<bool, boost::dynamic_bitset<> > > Q)
						{ 
							std::vector<boost::dynamic_bitset<>::block_type, boost::dynamic_bitset<>::allocator_type>* q;
							to_block_range(Q.second.second, &q);
							size_t* qq; 
							to_block_range(Q.second.second, &qq);
							ar & Q.first & Q.second.first & (*q) & (*qq);
						}
						);
	}
};


/// @brief a specialization version of the CompressRuleTrait with its template parameter CompressType specialized to be CompressType::LOWERCASE.  \n For LOWERCASE situation with N=2, lowercase characters a, c, g, and t are encoded into two-bits sized data format, while all the other characters are handled as except.
template <>
struct CompressRuleTrait <2, CompressType::LOWERCASE>
{
	static const int type = CompressType::LOWERCASE;
	static const int bit_count = 2;
	std::string alphabet;
	std::unordered_map < char, std::pair< bool, boost::dynamic_bitset<> > > table;
	CompressRuleTrait ()
		: table ({
	   		{'a', std::make_pair (true, boost::dynamic_bitset<> (2, 0ul))},
 			{'c', std::make_pair (true, boost::dynamic_bitset<> (2, 1ul))},
		   	{'g', std::make_pair (true, boost::dynamic_bitset<> (2, 2ul))},
			{'t', std::make_pair (true, boost::dynamic_bitset<> (2, 3ul))},
			{'u', std::make_pair (false,boost::dynamic_bitset<> (2, 3ul))},
				})
	{
		std::for_each ( table.begin(), table.end(),
						[&] (const std::pair < char, std::pair<bool,boost::dynamic_bitset<> > > Q)
						{ alphabet+= Q.first;}
					);
	}
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize (Archive& ar, const unsigned int)
	{
		std::for_each (table.begin(), table.end(),
						[&] (std::pair <char, std::pair<bool, boost::dynamic_bitset<> > > Q)
						{ 
							std::vector<boost::dynamic_bitset<>::block_type, boost::dynamic_bitset<>::allocator_type>* q;
							to_block_range(Q.second.second, &q);
							size_t* qq;
							to_block_range(Q.second.second, &qq);
							ar & Q.first & Q.second.first & (*q) & (*qq);
						}
						);
	}
};

/// @brief a specialization version of the CompressRuleTrait with its template parameter CompressType specialized to be CompressType::N_BITS.  \n For N_BITS situationwith N=6, input characters are segmented into two different groups, target and except. Characters of the target group is encoded with six-bit scheme, while the except group is handled as except. Thus, two different data processes are employed to respectively taking care of the target and the except groups.
template <>
struct CompressRuleTrait <6, CompressType::N_BITS >
{
	static const int type = CompressType::N_BITS;
	static const int bit_count = 6;
	std::string alphabet;
	std::unordered_map < char, std::pair< bool, boost::dynamic_bitset<> > > table;
	CompressRuleTrait ()
		: table ({
	   		{'A', std::make_pair (true, boost::dynamic_bitset<> (6, 0ul))},
	   		{'a', std::make_pair (false,boost::dynamic_bitset<> (6, 0ul))},
 			{'C', std::make_pair (true, boost::dynamic_bitset<> (6, 1ul))},
 			{'c', std::make_pair (false,boost::dynamic_bitset<> (6, 1ul))},
		   	{'G', std::make_pair (true, boost::dynamic_bitset<> (6, 2ul))},
		   	{'g', std::make_pair (false,boost::dynamic_bitset<> (6, 2ul))},
			{'T', std::make_pair (true, boost::dynamic_bitset<> (6, 3ul))},
			{'t', std::make_pair (false,boost::dynamic_bitset<> (6, 3ul))},
			{'U', std::make_pair (false,boost::dynamic_bitset<> (6, 3ul))},
			{'u', std::make_pair (false,boost::dynamic_bitset<> (6, 3ul))}
				})
	{
		std::for_each ( table.begin(), table.end(),
						[&] (const std::pair < char, std::pair<bool,boost::dynamic_bitset<> > > Q)
						{ alphabet+= Q.first;}
					);
	}

	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & alphabet;
		std::for_each (table.begin(), table.end(),
						[&] (std::pair <char, std::pair<bool, boost::dynamic_bitset<> > > Q)
						{ 
							std::vector<boost::dynamic_bitset<>::block_type, boost::dynamic_bitset<>::allocator_type>* q;
							to_block_range(Q.second.second, &q);
							size_t* qq; 
							to_block_range(Q.second.second, &qq);
							ar & Q.first & Q.second.first & (*q) & (*qq);
						}
						);
	}
};

/// @brief a specialization version of the CompressRuleTrait with its template parameter CompressType specialized to be CompressType::N_BITS_MASK.  \n For N_BITS_MASK situation with N=6, an extra mask rule, such as lowercase character, is employed by the NBitCompress, so tht three different data processes are employed to respectively taking care of the target, the except, and the mask characters.
template <>
struct CompressRuleTrait <6, CompressType::N_BITS_MASK >
{
	static const int type = CompressType::N_BITS_MASK;
	static const int bit_count = 6;
//	static const uint8_t Value = 1;
	std::string alphabet;
	std::unordered_map < char, std::pair< bool, boost::dynamic_bitset<> > > table;
	CompressRuleTrait ()
		: table ({
	   		{'A', std::make_pair (true, boost::dynamic_bitset<> (6, 0ul))},
	   		{'a', std::make_pair (false,boost::dynamic_bitset<> (6, 0ul))},
 			{'C', std::make_pair (true, boost::dynamic_bitset<> (6, 1ul))},
			{'c', std::make_pair (false,boost::dynamic_bitset<> (6, 1ul))},
		   	{'G', std::make_pair (true, boost::dynamic_bitset<> (6, 2ul))},
		   	{'g', std::make_pair (false,boost::dynamic_bitset<> (6, 2ul))},
			{'T', std::make_pair (true, boost::dynamic_bitset<> (6, 3ul))},
			{'t', std::make_pair (false,boost::dynamic_bitset<> (6, 3ul))},
			{'U', std::make_pair (false,boost::dynamic_bitset<> (6, 3ul))},
			{'u', std::make_pair (false,boost::dynamic_bitset<> (6, 3ul))}
//	   		{'A', std::make_pair (true, boost::dynamic_bitset<> (6, 0ul))},
 //			{'C', std::make_pair (true, boost::dynamic_bitset<> (6, 1ul))},
//		   	{'G', std::make_pair (true, boost::dynamic_bitset<> (6, 2ul))},
//			{'T', std::make_pair (true, boost::dynamic_bitset<> (6, 3ul))},
//			{'U', std::make_pair (false, boost::dynamic_bitset<> (6, 3ul))}
				})
	{
		std::for_each ( table.begin(), table.end(),
						[&] (const std::pair < char, std::pair<bool, boost::dynamic_bitset<> > > Q)
						{ alphabet += Q.first;}
						);
	}
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & alphabet;
		std::for_each (table.begin(), table.end(),
						[&] (std::pair <char, std::pair<bool, boost::dynamic_bitset<> > > Q)
						{ 
							std::vector<boost::dynamic_bitset<>::block_type, boost::dynamic_bitset<>::allocator_type>* q;
							to_block_range(Q.second.second, &q);
							size_t* qq; 
							to_block_range(Q.second.second, &qq);
							ar & Q.first & Q.second.first & (*q) & (*qq);
						}
						);
	}
};

/// @brief a specialization version of the CompressRuleTrait with its template parameter CompressType specialized to be CompressType::LOWERCASE.  \n For LOWERCASE situation with N=6, lowercase characters a, c, g, and t are encoded into two-bits sized data format, while all the other characters are handled as except.
template <>
struct CompressRuleTrait <6, CompressType::LOWERCASE>
{
	static const int type = CompressType::LOWERCASE;
	static const int bit_count = 6;
	std::string alphabet;
	std::unordered_map < char, std::pair< bool, boost::dynamic_bitset<> > > table;
	CompressRuleTrait ()
		: table ({
	   		{'a', std::make_pair (true, boost::dynamic_bitset<> (6, 0ul))},
 			{'c', std::make_pair (true, boost::dynamic_bitset<> (6, 1ul))},
		   	{'g', std::make_pair (true, boost::dynamic_bitset<> (6, 2ul))},
			{'t', std::make_pair (true, boost::dynamic_bitset<> (6, 3ul))},
			{'u', std::make_pair (false,boost::dynamic_bitset<> (6, 3ul))},
				})
	{
		std::for_each ( table.begin(), table.end(),
						[&] (const std::pair < char, std::pair<bool,boost::dynamic_bitset<> > > Q)
						{ alphabet+= Q.first;}
					);
	}
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize (Archive& ar, const unsigned int)
	{
		std::for_each (table.begin(), table.end(),
						[&] (std::pair <char, std::pair<bool, boost::dynamic_bitset<> > > Q)
						{ 
							std::vector<boost::dynamic_bitset<>::block_type, boost::dynamic_bitset<>::allocator_type>* q;
							to_block_range(Q.second.second, &q);
							size_t* qq;
							to_block_range(Q.second.second, &qq);
							ar & Q.first & Q.second.first & (*q) & (*qq);
						}
						);
	}
};

template <int N, CompressType C>
bool operator== (const CompressRuleTrait<N, C>& q, const CompressRuleTrait<N, C>& qq)
{
	if ( q.type == qq.type && q.bit_count == qq.bit_count && q.alphabet == qq.alphabet && q.table == qq.table )
		return true;
	else
		return false;
}
#endif
