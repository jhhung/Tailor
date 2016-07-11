Tailor
=========
Tailor is a Burrows–Wheeler transform (BWT) based fast aligner specialized in discovering non-templated addition of nucleotide to the 3' end of nucleic acids (a.k.a., tailing) from Next Generation Sequencing data.  

Tailor is released under GPLv2 with additional restriction so that is only applicable to individuals and non-profits and that any for-profit company must purchase a different license.   

A Shell based pipelines are provided using fastq as input and produce publication quality figures.    

## INSTALL
*Only 64 bits systems are able to compile and run Tailor.
    
### Run the binary directly without installation 
Please try the precompiled binaries first, most of the linux systems should be able to run Tailor without any troubles.
```bash
bin/tailor_v.1.0_linux_static 		# for linux
```
Or you can find it in the release tab in this page or at this [link](https://github.com/jhhung/Tailor/releases).

### Compile from the source code
#### Install the dependencies
- 1.1 Relative recent C++ compiler that support most features of C++11. We recommend [GCC](http://gcc.gnu.org/).
- 1.2 [Boost](http://www.boost.org/users/download/)
- 1.3 [CMake](http://www.cmake.org/)

#### Get the latest version of the software
```bash
git clone git@github.com:jhhung/Tailor.git
```

#### Enter the folder Tailor and:
- Set environmental variable `$BOOST_ROOT` to the directory of boost if CMake cannot find boost automatically;
- Set environmental variable `$CC` and `$CXX` to the gcc/g++ compiler you want to use.	
```bash
cmake .
```

#### Compile the software by typing:
```bash
make
```

#### troubleshooting
- If you got linker error, it is possible that the default library in the lib/ is not suitable to your platform. 
 There is one library available for Linux, rename the one that fit the best to "libabwt_table.a",
 and retype 
```bash
make
```
	
## USAGE
***
### tailor

#### Build genomic index (similar to bowtie-build)
Options:
```bash
-h [ --help ]         display this help message and exit
-i [ --input ] arg    The input fasta file.
-p [ --prefix ] arg   Prefix of index file to generate.
-f [ --force ]        Overwrite the existing index files if they already exist.
```
Example:
```
tailor build -i genome.fa -p genome
```

#### Mapping 
```bash

-h [ --help ]                 display this help message and exit
-i [ --input ] arg            Input fastq file
-p [ --index ] arg            Prefix of the index
-o [ --output ] arg (=stdout) Output SAM file, stdout by default 
-n [ --thread ] arg (=1)      Number of thread to use; if the number is larger than the core available, it will be adjusted automatically
-l [ --minLen ] arg (=18)     minimal length of exact match (prefix match) allowed
-v [ --mismatch ]             to allow mismatch in the middle of the query
```
Examples:
```
tailor map -p genome -n 8 -i smallRNA.fq -o smallRNA.sam         # no mismatch
tailor map -p genome -n 8 -i smallRNA.fq -o smallRNA.sam   -v    # allow one mismatch
```

***
### tailing pipeline
```bash
# reads.fq: input reads in fastq format: 
# dm3.fa: genome sequence in fasta format; it will generate index in the index folder of tailor directory, if it doesn't exist
# genomic_feature_file; used to generate figures for different genomic features (exon, intron...). See our example in the annotation folder to make such file
# using 24 CPUs
# use PHRED score 20 as filter: only reads with every base equal to to higher than 20 pass the filter and enters the pipeline
run_miRNA_tailing_pipeline.sh \ 
	-i reads.fq  \ 
	-g dm3.fa \ 
	-t genomic_feature_file \ 
	-o output_dir \ 
	-c 24 \ 
	-q 20
```

## Download

***
#### Indexes

- Tailor indexes

```bash
# You can find all the pre-bulit indexes in:
http://www.jhhlab.tw/Tailor/index/

# Human:
http://www.jhhlab.tw/Tailor/index/hg18.tar.gz
http://www.jhhlab.tw/Tailor/index/hg19.tar.gz

# Mouse:
http://www.jhhlab.tw/Tailor/index/mm9.tar.gz
http://www.jhhlab.tw/Tailor/index/mm10.tar.gz

# for downloading
lftp -c "pget -n 4 http://www.jhhlab.tw/Tailor/index/hg18.tar.gz"
```

***
#### Speed test files for the publication
- You can download all the related files for the speed test from the link
```
http://www.jhhlab.tw/Tailor/speed_test_samples/
```
- And the links of original data for non-tailed and tailed reads
```
http://www.jhhlab.tw/Tailor/speed_test_samples/Drosophila_melanogaster.2m.fq
http://www.jhhlab.tw/Tailor/speed_test_samples/Drosophila_melanogaster.all.randomeTailed.fq
```
- And the speed test log (3 times)
```
http://www.jhhlab.tw/Tailor/speed_test_samples/test_speed.log
http://www.jhhlab.tw/Tailor/speed_test_samples/test_speed.log2
http://www.jhhlab.tw/Tailor/speed_test_samples/test_speed.log3
```
- And the speed log for bowtie tailing
```
http://www.jhhlab.tw/Tailor/speed_test_samples/tailing.log
http://www.jhhlab.tw/Tailor/speed_test_samples/tailing.log2
http://www.jhhlab.tw/Tailor/speed_test_samples/tailing.log3
```
- All scripts for speed test can be found in the `utils` directory

## Citing Tailor
* Chou, M.-T., Han, B. W., Hsiao, C.-P., Zamore, P. D., Weng, Z., and Hung, J.-H. (2015). [Tailor: a computational framework for detecting non-templated tailing of small silencing RNAs.](http://www.ncbi.nlm.nih.gov/pubmed/26007652) *Nucleic Acids Res*. 43, e109.

#### to cite annotations included in `Tailor`
* If you have used the fly piRNA cluster annotation, please cite:
    * Brennecke, J., Aravin, A. A., Stark, A., Dus, M., Kellis, M., Sachidanandam, R., and Hannon, G. J. (2007). [Discrete small RNA-generating loci as master regulators of transposon activity in Drosophila.](http://www.ncbi.nlm.nih.gov/pubmed/17346786) *Cell* 128, 1089-1103.
* If you have used the mouse piRNA cluster annotation, please cite:
    * Li, X. Z., Roy, C. K., Dong, X., Bolcun-Filas, E., Wang, J., Han, B. W., Xu, J., Moore, M. J., Schimenti, J. C., Weng, Z., and Zamore, P. D. (2013). [An ancient transcription factor initiates the burst of piRNA production during early meiosis in mouse testes.](http://www.ncbi.nlm.nih.gov/pubmed/23523368) *Mol Cell* 50, 67-81.

## Contact
```bash
    Jui-Hung Hung <juihunghung `at` gmail.com>
    Bo W Han <bowhan `at` me.com>
    Chiung-Po Hsiao <restart0216s `at` gmail.com>
```
