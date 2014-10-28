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

#ifndef JBIT_HPP_
#define JBIT_HPP_
#include <map>
#include <unordered_map>
#include "boost/unordered_map.hpp"
#include <locale>
#include <algorithm>
#include <bitset>
#include <iostream>
#include <array>
#include <fstream>
#include <deque>
#include <string>
#include <set>
//#include "../constant_def.hpp"
#include <cctype>
#include <bitset>
///this class is currently used for only JHH's personal attempt of expediting BWT

class Byte2word
{
public:
	Byte2word()
	{
		std::string bits;
		std::string words;
		std::vector<std::string> all_char     {"A", "C", "G", "T", "$"};
		std::vector<std::string> all_char_bit {"00", "11", "01", "10", "11"};
		//XX$, X$$, $$$
		/*
		GGG e
		GGT f
		GTA g

		{{ "AAA", '!' }, {"AAT", '?'} };
		*/
		std::cerr<<"generating src/three_byte_to_char.hpp"<<std::endl;
		std::ofstream out ("src/three_byte_to_char.hpp");
		char asci = 0x21;
		out<<"{"<<std::endl<<"\t";
		for (int i = 0; i < 4; i++)//first nt
			{
				words = std::string("{{ '")+all_char[i] + "', '$'" + ", '$'}}" ;
//				std::cerr<<words<<" "<< ++asci<<std::endl;
				out<<"{ '"<<++asci<<"', "<<words<<"}, ";
			}
			out<<std::endl<<"\t";
		asci = 0x25;
		for (int i = 0; i < 4; i++)//first nt
		{
		  for (int j = 0; j < 4; j++)//second nt
			{
				words = std::string("{{ '")+all_char[i] +"', '"+ all_char[j] + "', '$'}}" ;
//				std::cerr<<words<<" "<< ++asci<<std::endl;
				out<<"{ '"<<++asci<<"', "<<words<<"}, ";

			}
			out<<std::endl<<"\t";
		}
		asci = 0x3A;
		for (int i = 0; i < 4; i++)//first nt
		  for (int j = 0; j < 4; j++)//second nt
		  {
			for (int k = 0; k < 4; k++)//third nt
			{
				words = std::string("{{ '")+all_char[i] + "', '" +all_char[j] + "', '"+all_char[k]+"'}}" ;
//				std::cerr<<words<<" "<< ++asci<<std::endl;
				out<<"{ '"<<++asci<<"', "<<words<<"}, ";
			}
			out<<std::endl<<"\t";
		  }
		out<<"{ '"<<(char)0x20<<"', {{ '$', '$', '$'}}"<<"}"<<std::endl;
		out<<"}"<<std::endl;
	}

};

class Word2byte
{
public:
	Word2byte()
	{
		std::string bits;
		std::string words;
		std::vector<std::string> all_char     {"A", "C", "G", "T", "$"};
		std::vector<std::string> all_char_bit {"00", "11", "01", "10", "11"};
		//XX$, X$$, $$$
		/*
		GGG e
		GGT f
		GTA g

		{{ "AAA", '!' }, {"AAT", '?'} };
		*/
		std::cerr<<"generating src/three_char_to_byte.hpp"<<std::endl;
		std::ofstream out ("src/three_char_to_byte.hpp");
		char asci = 0x21;
		out<<"{"<<std::endl<<"\t";
		for (int i = 0; i < 4; i++)//first nt
			{
				words = std::string("{{ '")+all_char[i] + "', '$'" + ", '$'}}" ;
//				std::cerr<<words<<" "<< ++asci<<std::endl;
				out<<"{ "<<words<<", '"<<++asci<<"'}, ";
			}
			out<<std::endl<<"\t";
		asci = 0x25;
		for (int i = 0; i < 4; i++)//first nt
		{
		  for (int j = 0; j < 4; j++)//second nt
			{
				words = std::string("{{ '")+all_char[i] +"', '"+ all_char[j] + "', '$'}}" ;
//				std::cerr<<words<<" "<< ++asci<<std::endl;
				out<<"{ "<<words<<", '"<<++asci<<"'}, ";

			}
			out<<std::endl<<"\t";
		}
		asci = 0x3A;
		for (int i = 0; i < 4; i++)//first nt
		  for (int j = 0; j < 4; j++)//second nt
		  {
			for (int k = 0; k < 4; k++)//third nt
			{
				words = std::string("{{ '")+all_char[i] + "', '" +all_char[j] + "', '"+all_char[k]+"'}}" ;
//				std::cerr<<words<<" "<< ++asci<<std::endl;
				out<<"{ "<<words<<", '"<<++asci<<"'}, ";
			}
			out<<std::endl<<"\t";
	      }
		out<<"{ "<<"{{ '$', '$', '$'}}"<<", '"<< (char)0x20 <<"'}"<<std::endl;
		out<<"}"<<std::endl;
	}

};


class JBit
{
public:
	//std::string seq_;
	std::vector<uint8_t> seq_;
	boost::unordered_map<std::array<char,4>, uint8_t> word2byte;
	boost::unordered_map<uint8_t, std::array<char,4> > byte2word;
	JBit( std::string& original_seq )
	{
		make_table();

		uint32_t tmp_i = 4 - ( (original_seq.length() & 3 ));

		//std::cerr << "tmp_i " << tmp_i << std::endl;
		for(int i=0; i < tmp_i; i++)
		{
			original_seq+="A";
		}

		seq_.reserve(original_seq.length()/4 +1);



		char* seq ( (char*)original_seq.c_str() );
//		for (int j=0; j <1000000; j++)
//		{
		for ( int i = 0; i<original_seq.length(); i+=4)
		{
			//if(i == 1377*4)
			//	std::cerr << (uint32_t) word2byte[  std::array<char, 4>{ *(seq+i), *(seq+i+1), *(seq+i+2), *(seq+i+3) } ] << std::endl;
			seq_.push_back( word2byte[  std::array<char, 4>{ *(seq+i), *(seq+i+1), *(seq+i+2), *(seq+i+3) } ]);
			//if(i == 1377*4)
			//	std::cerr << (uint32_t) seq_[1377] << std::endl;
		}
		//std::cerr << "jseq size " << seq_.size() << " " << (int)seq_.back() <<  std::endl;;

		//std::string result;
		//result.reserve(original_seq.length());
		//for (auto k : seq_)
		//{
		//	result += byte2word[k].data();// = result + byte2word[k][0] +  byte2word[k][1] + byte2word[k][2];
		//}

		//if (result != original_seq)
		//	std::cerr<<"what?"<<std::endl;
//		seq_.clear();
//		}
		//std::cerr<<result<<std::endl;
	}
	void make_table()
	{
		std::string all_char("ACGT");
		for(int i(0); i<256; i++)
		{
			char w1 = all_char[ (i&255) >> 6 ];
			char w2 = all_char[ (i&63) >> 4 ];
			char w3 = all_char[ (i&15) >> 2 ];
			char w4 = all_char[ (i&3) >> 0 ];
			//std::cerr << i <<" : " << w1 << w2 << w3 << w4 << std::endl;
			word2byte.insert( { std::array<char, 4>{w1,w2,w3,w4}, i } );
			byte2word.insert( { i, std::array<char, 4>{w1,w2,w3,w4} } );
		}
	}

};


/*
//std::unordered_map<std::string, char> JBit::word2byte
boost::unordered_map< char, std::array<char, 3> > JBit::byte2word
{
    { '\"',{{ 'A', '$', '$'}}}, { '#', {{ 'C', '$', '$'}}}, { '$', {{ 'G', '$', '$'}}}, { '%', {{ 'T', '$', '$'}}},
	{ '&', {{ 'A', 'A', '$'}}}, { '\'',{{ 'A', 'C', '$'}}}, { '(', {{ 'A', 'G', '$'}}}, { ')', {{ 'A', 'T', '$'}}},
	{ '*', {{ 'C', 'A', '$'}}}, { '+', {{ 'C', 'C', '$'}}}, { ',', {{ 'C', 'G', '$'}}}, { '-', {{ 'C', 'T', '$'}}},
	{ '.', {{ 'G', 'A', '$'}}}, { '/', {{ 'G', 'C', '$'}}}, { '0', {{ 'G', 'G', '$'}}}, { '1', {{ 'G', 'T', '$'}}},
	{ '2', {{ 'T', 'A', '$'}}}, { '3', {{ 'T', 'C', '$'}}}, { '4', {{ 'T', 'G', '$'}}}, { '5', {{ 'T', 'T', '$'}}},
	{ ';', {{ 'A', 'A', 'A'}}}, { '<', {{ 'A', 'A', 'C'}}}, { '=', {{ 'A', 'A', 'G'}}}, { '>', {{ 'A', 'A', 'T'}}},
	{ '?', {{ 'A', 'C', 'A'}}}, { '@', {{ 'A', 'C', 'C'}}}, { 'A', {{ 'A', 'C', 'G'}}}, { 'B', {{ 'A', 'C', 'T'}}},
	{ 'C', {{ 'A', 'G', 'A'}}}, { 'D', {{ 'A', 'G', 'C'}}}, { 'E', {{ 'A', 'G', 'G'}}}, { 'F', {{ 'A', 'G', 'T'}}},
	{ 'G', {{ 'A', 'T', 'A'}}}, { 'H', {{ 'A', 'T', 'C'}}}, { 'I', {{ 'A', 'T', 'G'}}}, { 'J', {{ 'A', 'T', 'T'}}},
	{ 'K', {{ 'C', 'A', 'A'}}}, { 'L', {{ 'C', 'A', 'C'}}}, { 'M', {{ 'C', 'A', 'G'}}}, { 'N', {{ 'C', 'A', 'T'}}},
	{ 'O', {{ 'C', 'C', 'A'}}}, { 'P', {{ 'C', 'C', 'C'}}}, { 'Q', {{ 'C', 'C', 'G'}}}, { 'R', {{ 'C', 'C', 'T'}}},
	{ 'S', {{ 'C', 'G', 'A'}}}, { 'T', {{ 'C', 'G', 'C'}}}, { 'U', {{ 'C', 'G', 'G'}}}, { 'V', {{ 'C', 'G', 'T'}}},
	{ 'W', {{ 'C', 'T', 'A'}}}, { 'X', {{ 'C', 'T', 'C'}}}, { 'Y', {{ 'C', 'T', 'G'}}}, { 'Z', {{ 'C', 'T', 'T'}}},
	{ '[', {{ 'G', 'A', 'A'}}}, { '\\',{{ 'G', 'A', 'C'}}}, { ']', {{ 'G', 'A', 'G'}}}, { '^', {{ 'G', 'A', 'T'}}},
	{ '_', {{ 'G', 'C', 'A'}}}, { '`', {{ 'G', 'C', 'C'}}}, { 'a', {{ 'G', 'C', 'G'}}}, { 'b', {{ 'G', 'C', 'T'}}},
	{ 'c', {{ 'G', 'G', 'A'}}}, { 'd', {{ 'G', 'G', 'C'}}}, { 'e', {{ 'G', 'G', 'G'}}}, { 'f', {{ 'G', 'G', 'T'}}},
	{ 'g', {{ 'G', 'T', 'A'}}}, { 'h', {{ 'G', 'T', 'C'}}}, { 'i', {{ 'G', 'T', 'G'}}}, { 'j', {{ 'G', 'T', 'T'}}},
	{ 'k', {{ 'T', 'A', 'A'}}}, { 'l', {{ 'T', 'A', 'C'}}}, { 'm', {{ 'T', 'A', 'G'}}}, { 'n', {{ 'T', 'A', 'T'}}},
	{ 'o', {{ 'T', 'C', 'A'}}}, { 'p', {{ 'T', 'C', 'C'}}}, { 'q', {{ 'T', 'C', 'G'}}}, { 'r', {{ 'T', 'C', 'T'}}},
	{ 's', {{ 'T', 'G', 'A'}}}, { 't', {{ 'T', 'G', 'C'}}}, { 'u', {{ 'T', 'G', 'G'}}}, { 'v', {{ 'T', 'G', 'T'}}},
	{ 'w', {{ 'T', 'T', 'A'}}}, { 'x', {{ 'T', 'T', 'C'}}}, { 'y', {{ 'T', 'T', 'G'}}}, { 'z', {{ 'T', 'T', 'T'}}},
	{ ' ', {{ '$', '$', '$'}}}
};



boost::unordered_map<std::array<char, 3>, char> JBit::word2byte
{
   { {{ 'A', '$', '$'}}, '\"'},{ {{ 'C', '$', '$'}}, '#'}, { {{ 'G', '$', '$'}}, '$'}, { {{ 'T', '$', '$'}}, '%'},
   { {{ 'A', 'A', '$'}}, '&'}, { {{ 'A', 'C', '$'}}, '\''},{ {{ 'A', 'G', '$'}}, '('}, { {{ 'A', 'T', '$'}}, ')'},
   { {{ 'C', 'A', '$'}}, '*'}, { {{ 'C', 'C', '$'}}, '+'}, { {{ 'C', 'G', '$'}}, ','}, { {{ 'C', 'T', '$'}}, '-'},
   { {{ 'G', 'A', '$'}}, '.'}, { {{ 'G', 'C', '$'}}, '/'}, { {{ 'G', 'G', '$'}}, '0'}, { {{ 'G', 'T', '$'}}, '1'},
   { {{ 'T', 'A', '$'}}, '2'}, { {{ 'T', 'C', '$'}}, '3'}, { {{ 'T', 'G', '$'}}, '4'}, { {{ 'T', 'T', '$'}}, '5'},
   { {{ 'A', 'A', 'A'}}, ';'}, { {{ 'A', 'A', 'C'}}, '<'}, { {{ 'A', 'A', 'G'}}, '='}, { {{ 'A', 'A', 'T'}}, '>'},
   { {{ 'A', 'C', 'A'}}, '?'}, { {{ 'A', 'C', 'C'}}, '@'}, { {{ 'A', 'C', 'G'}}, 'A'}, { {{ 'A', 'C', 'T'}}, 'B'},
   { {{ 'A', 'G', 'A'}}, 'C'}, { {{ 'A', 'G', 'C'}}, 'D'}, { {{ 'A', 'G', 'G'}}, 'E'}, { {{ 'A', 'G', 'T'}}, 'F'},
   { {{ 'A', 'T', 'A'}}, 'G'}, { {{ 'A', 'T', 'C'}}, 'H'}, { {{ 'A', 'T', 'G'}}, 'I'}, { {{ 'A', 'T', 'T'}}, 'J'},
   { {{ 'C', 'A', 'A'}}, 'K'}, { {{ 'C', 'A', 'C'}}, 'L'}, { {{ 'C', 'A', 'G'}}, 'M'}, { {{ 'C', 'A', 'T'}}, 'N'},
   { {{ 'C', 'C', 'A'}}, 'O'}, { {{ 'C', 'C', 'C'}}, 'P'}, { {{ 'C', 'C', 'G'}}, 'Q'}, { {{ 'C', 'C', 'T'}}, 'R'},
   { {{ 'C', 'G', 'A'}}, 'S'}, { {{ 'C', 'G', 'C'}}, 'T'}, { {{ 'C', 'G', 'G'}}, 'U'}, { {{ 'C', 'G', 'T'}}, 'V'},
   { {{ 'C', 'T', 'A'}}, 'W'}, { {{ 'C', 'T', 'C'}}, 'X'}, { {{ 'C', 'T', 'G'}}, 'Y'}, { {{ 'C', 'T', 'T'}}, 'Z'},
   { {{ 'G', 'A', 'A'}}, '['}, { {{ 'G', 'A', 'C'}}, '\\'},{ {{ 'G', 'A', 'G'}}, ']'}, { {{ 'G', 'A', 'T'}}, '^'},
   { {{ 'G', 'C', 'A'}}, '_'}, { {{ 'G', 'C', 'C'}}, '`'}, { {{ 'G', 'C', 'G'}}, 'a'}, { {{ 'G', 'C', 'T'}}, 'b'},
   { {{ 'G', 'G', 'A'}}, 'c'}, { {{ 'G', 'G', 'C'}}, 'd'}, { {{ 'G', 'G', 'G'}}, 'e'}, { {{ 'G', 'G', 'T'}}, 'f'},
   { {{ 'G', 'T', 'A'}}, 'g'}, { {{ 'G', 'T', 'C'}}, 'h'}, { {{ 'G', 'T', 'G'}}, 'i'}, { {{ 'G', 'T', 'T'}}, 'j'},
   { {{ 'T', 'A', 'A'}}, 'k'}, { {{ 'T', 'A', 'C'}}, 'l'}, { {{ 'T', 'A', 'G'}}, 'm'}, { {{ 'T', 'A', 'T'}}, 'n'},
   { {{ 'T', 'C', 'A'}}, 'o'}, { {{ 'T', 'C', 'C'}}, 'p'}, { {{ 'T', 'C', 'G'}}, 'q'}, { {{ 'T', 'C', 'T'}}, 'r'},
   { {{ 'T', 'G', 'A'}}, 's'}, { {{ 'T', 'G', 'C'}}, 't'}, { {{ 'T', 'G', 'G'}}, 'u'}, { {{ 'T', 'G', 'T'}}, 'v'},
   { {{ 'T', 'T', 'A'}}, 'w'}, { {{ 'T', 'T', 'C'}}, 'x'}, { {{ 'T', 'T', 'G'}}, 'y'}, { {{ 'T', 'T', 'T'}}, 'z'},
   { {{ '$', '$', '$'}}, ' '}
};
*/
/*
boost::unordered_map<std::string, char> JBit::word2byte
{
    {"A$$", '\"'},{"C$$", '#'},{"G$$", '$'},{"T$$", '%'},
	{"AA$", '&'},{"AC$", '\''},{"AG$", '('},{"AT$", ')'},
	{"CA$", '*'},{"CC$", '+'},{"CG$", ','},{"CT$", '-'},
	{"GA$", '.'},{"GC$", '/'},{"GG$", '0'},{"GT$", '1'},
	{"TA$", '2'},{"TC$", '3'},{"TG$", '4'},{"TT$", '5'},
	{"AAA", ';'},{"AAC", '<'},{"AAG", '='},{"AAT", '>'},
	{"ACA", '?'},{"ACC", '@'},{"ACG", 'A'},{"ACT", 'B'},
	{"AGA", 'C'},{"AGC", 'D'},{"AGG", 'E'},{"AGT", 'F'},
	{"ATA", 'G'},{"ATC", 'H'},{"ATG", 'I'},{"ATT", 'J'},
	{"CAA", 'K'},{"CAC", 'L'},{"CAG", 'M'},{"CAT", 'N'},
	{"CCA", 'O'},{"CCC", 'P'},{"CCG", 'Q'},{"CCT", 'R'},
	{"CGA", 'S'},{"CGC", 'T'},{"CGG", 'U'},{"CGT", 'V'},
	{"CTA", 'W'},{"CTC", 'X'},{"CTG", 'Y'},{"CTT", 'Z'},
	{"GAA", '['},{"GAC", '\\'},{"GAG", ']'},{"GAT", '^'},
	{"GCA", '_'},{"GCC", '`'},{"GCG", 'a'},{"GCT", 'b'},
	{"GGA", 'c'},{"GGC", 'd'},{"GGG", 'e'},{"GGT", 'f'},
	{"GTA", 'g'},{"GTC", 'h'},{"GTG", 'i'},{"GTT", 'j'},
	{"TAA", 'k'},{"TAC", 'l'},{"TAG", 'm'},{"TAT", 'n'},
	{"TCA", 'o'},{"TCC", 'p'},{"TCG", 'q'},{"TCT", 'r'},
	{"TGA", 's'},{"TGC", 't'},{"TGG", 'u'},{"TGT", 'v'},
	{"TTA", 'w'},{"TTC", 'x'},{"TTG", 'y'},{"TTT", 'z'},
	{"$$$", ' '}
};
*/
#endif
