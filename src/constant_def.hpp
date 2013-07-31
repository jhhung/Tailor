/// @file constant_def.hpp
/// @brief provide enum format_types, CompressType, and WrapperType as navigating indexes for respectively achieving tag techniques on specialized FileReader, compression, and Wrapper template classes
/// @author JHH Corp.
#ifndef CONSTANT_DEF_HPP_
#define CONSTANT_DEF_HPP_

#include <string>

#ifndef myALPHABET
#define myALPHABET
std::string ALPHABET_$ACGNT("$ACGNT");
std::string ALPHABET_3CHAR(" \"#$%&\'()*+,-./012345;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz");
#endif

typedef uint64_t INTTYPE;
#define gAlphabet &ALPHABET_$ACGNT
///@enum format_types
///@brief format_types, assigning enumerated int, e.g. 0, 1, 2, 3, ... for indicating different formats, e.g. FastA, FastQ, Bed, Wig, ... .
enum format_types
{
    FASTA,
    FASTQ,
    BED,
	WIG,
	SAM
};
///@enum CompressTypes
///@brief CompressType, assigning enumerated int, e.g. 0, 1, ... for indicating different compression schemes, e.g. Plain, TWO_BITS, ... .
enum CompressType
{
    Plain,
    N_BITS, ///two bits, A->00 C->01 G->10 T->11, the order of binary is the same as the lexigraphical order
	N_BITS_MASK,
	LOWERCASE
};
///@enum WrapperType
///@brief WrapperType, assigning enumerated int, e.g. 0, 1, ... for indicating different wrapper formats, e.g. TUPLE_WRAPPER, VECTOR_WRAPPER, ... .
enum WrapperType
{
	TUPLE_WRAPPER = 1,
	VECTOR_WRAPPER =2,
};
///@enum ParallelTypes
///@brief ParallelType, assigning enumerated int for indicating different file_reader_intf policies, e.g. normal, multi-thread, ... .
enum ParallelTypes
{
	NORMAL,
	M_T,
	M_P_I
};

enum TwoBit_types
{   
    Normal,
    Mask
};

enum CompressFormat
{
    GZ,
    PLAIN,
};

enum SOURCE_TYPE
{
    IFSTREAM_TYPE,
    CURL_TYPE_GZ,
    CURL_TYPE   //not tested yet
};

#endif
