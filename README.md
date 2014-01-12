Tailor
=========
Tailor program is a Burrows–Wheeler transform based fast short read aligner (like bowtie/BWA/SOAP2) specialized in discovering trimming and tailing events for small silencing RNAs from Next Generation Sequencing data. 
Two shell based pipelines for miRNA or small RNA in general are included in the **utils** directory.

Tailor is released under GPLv2 with additional restriction so that is only applicable to individuals and non-profits and that any for-profit company must purchase a different license.

##INSTALL

*Only 64 bits systems are able to compile and run Tailor. 

#### Run the binary directly without instaillation (if you are lucky)

Try the precompiled binaries first, most of the linux systems should be able to run Tailor without any troubles.

Please find the binary that suits your platform:

```
bin/tailor_debian_x64  // for ubuntu/fedara/...
bin/tailor_redhat_x64  // for centos/redhat/...
bin/tailor_mac_x64     // for OSX
```

Or you can find them in the release tab in this page or at this link:
https://github.com/jhhung/Tailor/releases

You can rename it to "Tailor" for your convinience.
If you don't know which one is good for your platform, just pick one and try it out.
If unfortunately, none of them works, please see below to build a binary for your box.

#### Install the dependencies

- 1.1 Relative recent C++ compiler that support most features of C++11. We recommend [GCC](http://gcc.gnu.org/).
- 1.2 [Boost](http://www.boost.org/users/download/)
- 1.3 [CMake](http://www.cmake.org/)

#### Get the latest version of the software

```
git clone git@github.com:jhhung/Tailor.git
```

#### Enter the folder Tailor and:

- Set enviromental variable "BOOST_ROOT" to the directory of boost if CMake cannot find boost automatically;
- Set enviromental variable "CC" and "CXX" to the gcc/g++ compiler you want to use.	

```
cmake .
```
   
	
#### Compile the software by typing:

```
make
```

#### troubleshooting
- If you got linker error, it is possible that the default library in the lib/ is not suitable to your platform. 
 There are two libraries available, one is for Mac OSX one is for Linux, rename the one that fit the best to "libabwt_table.a",
 and retype 

```
make
```
	
##USAGE

#### Build genomic index (similar to bowtie-build)

```
tailor build -i genome.fa -p genome
```

#### Mapping 

```
tailor map -p genome -n 8 -i smallRNA.fq
```

####microRNA pipeline

- Modify the $PATH to include the Tailor/utils directory

```
run_miRNA_tailing_pipeline.sh -i reads.fq -m miRNA.fa -p hairpin.fa
```

- miRNA.fa contains the mature miRNA sequences of a specific organism (from mirBase)
- hairpin.fa contains the hairpin miRNA sequences of a specific organism (from mirBase)


#### general small RNA pipeline

- The general idea of this pipeline is to perform genomic mapping of small RNA to genome with Tailor
then the reads will be assigned to different genomic feature using BEDtools
- You will need to modify `intersect_all_tailor.sh` to include the BED files you would like to include
examples have been given in that script

```
run_tailing_pipeline.sh -i reads.fq -s mouse/fly 
```

- -s mouse/fly choose the organism you are working on
- In `run_tailing_pipeline.sh` and `intersect_all_tailor.sh`, there are "case" switch to choose different behavior according to the organism. 
- You will need to modify this according to your need. Examples are given in the scripts. 

##Citing Tailor
* not yet

##Contact
	Jui-Hung Hung <juihunghung@gmail.com>
	Chou Min-Te <poi5305@gmail.com>
	Bo W Han <bowhan@me.com>
