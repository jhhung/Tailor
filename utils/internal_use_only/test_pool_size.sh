#! /bin/bash -x

# This script is associated with the Tailor software and 
# only used for internal testing purpose and to repeat the 
# result in the published paper. 

for POOL_SIZE in 1000 2000 4000 6000 8000 10000 20000;
do 
	/home/hanb/bin/c++  -DPOOLSIZE=$POOL_SIZE -DmyALPHABET -Ofast -std=c++11 -lrt -I/home/hanb/nearline/boost_1_54_0    -o CMakeFiles/bin/tailor.dir/src/main.cpp.o -c /home/hanb/src/Tailor/src/main.cpp
	make VERBOSE=1
	echo $POOL_SIZE
	time tailor map -i ~/scratch/tailor/speedTest/Phil.hen1Muts.trimmed.p10.fq -p /home/hanb/nearline/Drosophila_melanogaster/UCSC/dm3/Sequence/tailorIndex/dm3 -n 8 -o ~/scratch/tailor/speedTest/Phil.hen1Muts.trimmed.p10.t.sam
done
