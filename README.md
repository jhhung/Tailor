###Tailor

Tailor program is a Burrows–Wheeler transform based fast short read aligner (like bowtie/BWA/SOAP2) specialized in discovering trimming and tailing events for small silencing RNAs from Next Generation Sequencing data. 
Two shell based pipelines for miRNA or small RNA in general are included in the **utils** directory.

Tailor is released under GPLv2. 

###INSTALL

1. Install the dependencies
	- 1.1 Relative recent C++ compiler that support most features of C++11. We recommend [GCC](http://gcc.gnu.org/).
	- 1.2 [Boost](http://www.boost.org/users/download/)
	- 1.3 [CMake](http://www.cmake.org/)

2. Get the latest version of the software

	git clone git@github.com:jhhung/Tailor.git

3. Enter the folder Tailor and type:

	cmake .
   
    Set enviromental variable "BOOST_ROOT" to the directory of boost if CMake cannot find boost automatically;
    Set enviromental variable "CC" and "CXX" to the gcc/g++ compiler you want to use.	
	
4. Compile the software by typing:

	make

5. troubleshooting
	- 1.1. If you got linker error, it is possible that the default library in the lib/ is not suitable to your platform.
	  There are two libraries available, one is for Mac OSX one is for Linux, rename the one that fit the best to "libabwt_table.a",
          and retype 
	make
	
###USAGE

1.	Build genomic index (similar to bowtie-build)

	tailor build
	tailor build -i genome.fa -p genome
	
2.  Mapping 

	tailor map -p genome -n 8 -i smallRNA.fq

###Citing Tailor

*not yet

###Contact
	Jui-Hung Hung <juihunghung@gmail.com>
	Chou Min-Te <poi5305@gmail.com>
	Bo W Han <bowhan@me.com>