#include<iostream>
#include<fstream>
#include<string>
#include <algorithm>

bool cop(std::string &s)
{
	bool is_N(false);
	std::transform(s.begin(), s.end(), s.begin(),
		[&is_N](char c)
		{
			switch (c)
			{
				case 'A':
					c = 'T';
					break;
				case 'C':
					c = 'G';
					break;
				case 'G':
					c = 'C';
					break;
				case 'T':
					c = 'A';
					break;
				case 'U':
					c = 'A';
					break;
				case 'N':
					is_N = true;
					break;
			}
			return c;
		}
	);
	return is_N;
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		std::cerr << "usage: ./fastq_trimer trim_5p_len trim_3p_len is_rc[0|1] input_file output_file" << std::endl;
		return 0;
	}
	int trim5_len = std::stoi(argv[1]); //3 len barcode
	int trim3_len = std::stoi(argv[2]); //30 len adapter
	int is_rc = std::stoi(argv[3]); // is reverse complement
	std::string infn (argv[4]); //input file name
	std::string outfn (argv[5]); //output file name
	
	std::cout << "trim5_len: " << trim5_len << " trim3_len: " << trim3_len << " is_rc: " << is_rc << std::endl;
	std::ifstream in(infn);
	std::ofstream out(outfn);
	
	std::string fastq("");
	bool is_N(false);
	std::string line;
	
	int fastq_count = 0;
	
	while(!in.eof())
	{
		//fastq 1st line
		std::getline(in, line);
		if(line == "")
			break;
		if(fastq_count == 0)
			fastq += line;
		else
			fastq += "\n"+line;
		
		//fastq 2st line
		std::getline(in, line);
		line = line.substr(trim5_len, line.size() - trim5_len - trim3_len);
		if(is_rc)
		{
			is_N = cop(line);
			std::reverse(line.begin(), line.end());
		}
		fastq += "\n"+line;
		
		//fastq 3st line
		std::getline(in, line);
		fastq += "\n"+line;
		
		//fastq 4st line
		std::getline(in, line);
		line = line.substr(trim5_len, line.size() - trim5_len - trim3_len);
		if(is_rc)
		{
			std::reverse(line.begin(), line.end());
		}
		fastq += "\n"+line;
		
		if(!is_N)
		{
			out << fastq;
			is_N = false;
			fastq = "";
		}
		++fastq_count;
	}
	std::cout << "fastq count: " << fastq_count << std::endl;
	in.close();
	out.close();
	return 0;
}
