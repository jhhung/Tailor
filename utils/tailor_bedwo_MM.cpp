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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp> // for split
#include <boost/filesystem.hpp>       // for exists

using namespace std;

const static int Columns_size  = 18; // number of fields for the input
const static int MD_string_col = 9;  // 0-base column for MD string

int main(int argc, char** argv)
{
	std::ios::sync_with_stdio(false);
	istream* in1{&std::cin};
	if (argc > 1 && strcmp(argv[1],"stdin") && strcmp(argv[1],"-"))
	{
		if (!boost::filesystem::exists(argv[1]))
		{
			cerr << "error: file " << argv[1] << " does not exist" << endl;
			return 1;
		}
		in1 = new ifstream{argv[1]};
	}
	vector<string> tokens;
	string line;
	while(getline(*in1, line))
	{
		boost::split(tokens, line, boost::is_any_of("\t"));
		if (tokens.size() != Columns_size)
		{
			cerr << "error: the number of field should be " << Columns_size << "; but got " << tokens.size() << " for line " << line << endl;
			return 2;
		}
		if (tokens[5] != "+") continue;
		if (tokens[MD_string_col] == "*")
		{
			cout << tokens[0];
		}
		else
		{
			string::iterator begin1 {tokens[MD_string_col].begin()};
			while(*begin1 <= '9' && *begin1 >= '0') ++begin1;
			cout << tokens[0] << '.' << string {tokens[MD_string_col].begin(), ++begin1};
		}
		// cout << '\t' << stoi(tokens[11]) - stoi(tokens[1]) << '\t' << stoi(tokens[2]) - stoi(tokens[12]);
		for (int i = 1; i < MD_string_col; ++i) cout << '\t' << tokens[i];
		cout << '\n';
	}
	if (in1 != &std::cin) static_cast<ifstream*>(in1)->close();
}