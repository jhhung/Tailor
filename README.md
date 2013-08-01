Tailor
======

Tailor program is a Burrows–Wheeler transform based fast short read aligner (like bowtie/BWA/SOAP2) specialized in discovering trimming and tailing events for small silencing RNAs from Next Generation Sequencing data. 

##INSTALL
=======
1. Install the dependencies
	- 1.1 Relative recent C++ compiler that support most features of C++11. We recommend [GCC](http://gcc.gnu.org/).
	- 1.2 [Boost](http://www.boost.org/users/download/)
	- 1.3 [CMake](http://www.cmake.org/)

2. Get the latest version of the software

	`git clone git@github.com:jhhung/Tailor.git`

3. Enter the folder Tailor and type:

	`cmake .`
   
   - [x] Set enviromental variable "BOOST_ROOT" to the directory of boost if CMake cannot find boost automatically;
   - [x] Set enviromental variable "CC" and "CXX" to the gcc/g++ compiler you want to use.	
	
4. Compile the software by typing:

	`make`
	
##USAGE
=====

1.	Build genomic index (similar to bowtie-build)

	`tailor build `
	
2.  Mapping 

	`tailor map`
