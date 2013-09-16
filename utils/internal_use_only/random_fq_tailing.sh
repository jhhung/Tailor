#! /bin/bash -x

# This script is associated with the Tailor software and 
# only used for internal testing purpose and to repeat the 
# result in the published paper. 

# randomly append 1-3 nucleotide to the end of fastq

gawk 'BEGIN{nt[0]="A"; nt[1]="C"; nt[2]="G"; nt[3]="T"; srand(systime() + PROCINFO["pid"]); } {ranLen=int(rand()*100)%4; ranTail=""; for (i=1;i<=ranLen;++i) ranTail=ranTail""nt[int(rand()*100)%4]; print "@"ranTail; getline; printf "%s%s", $0, ranTail; getline; getline; printf "\n+\n%s", $0; for (i=1;i<=ranLen;++i) printf "I"; printf "\n"}' $1 > ${1%.f[aq*]}.randomeTailed.fq
