#! /bin/bash -x

# This script is associated with the Tailor software and 
# only used for internal testing purpose and to repeat the 
# result in the published paper. 

# map zebrafish

#bin
bowtie="/home/andy/andy/bowtie-1.0.0/bowtie";
bwa="/home/andy/andy/bwa-0.7.4/bwa";
tailor="/home/andy/andy/Tailor/bin/tailor";

# tailor index
TAILOR_INDEX=~/andy/Tailor/utils/internal_use_only/test_index/Drosophila_melanogaster/tailorIndex/dm3
# bowtie index
BOWTIE_INDEX=~/andy/Tailor/utils/internal_use_only/test_index/Drosophila_melanogaster/BowtieIndex/genome
# bwa index
BWA_INDEX=~/andy/Tailor/utils/internal_use_only/test_index/Drosophila_melanogaster/BWAIndex/genome.fa

FQ=Drosophila_melanogaster.all.randomeTailed.fq
for CPU in 2 4 8 12 24; do 
	time $bowtie -p $CPU -S -v 0 -a $BOWTIE_INDEX $FQ > ${FQ%fq}b1.sam #2>${FQ%fq}b1.log
	time $bowtie -p $CPU -S -v 0 --best --strata -a $BOWTIE_INDEX $FQ > ${FQ%fq}b1.sam #2>${FQ%fq}b1.log
	time ( $bwa aln -t $CPU -o0 -R1000000 $BWA_INDEX $FQ > ${FQ%fq}sai && $bwa samse $BWA_INDEX ${FQ%fq}sai $FQ > ${FQ%fq}bwa.sam ) 
	time $tailor map -i $FQ -p $TAILOR_INDEX -n $CPU -o ${FQ%fq}t.sam
done

FQ=Drosophila_melanogaster.2m.fq
for CPU in 2 4 8 12 24; do
	time $bowtie -p $CPU -S -v 0 -a $BOWTIE_INDEX $FQ > ${FQ%fq}b1.sam #2>${FQ%fq}b1.log
	time $bowtie -p $CPU -S -v 0 --best --strata -a $BOWTIE_INDEX $FQ > ${FQ%fq}b1.sam #2>${FQ%fq}b1.log
	time ( $bwa aln -t $CPU -o0 -R1000000 $BWA_INDEX $FQ > ${FQ%fq}sai && $bwa samse $BWA_INDEX ${FQ%fq}sai $FQ > ${FQ%fq}bwa.sam ) 
	time $tailor map -i $FQ -p $TAILOR_INDEX -n $CPU -o ${FQ%fq}t.sam
done

exit;

# map fruitfly

# tailor index
TAILOR_INDEX=/home/hanb/nearline/Drosophila_melanogaster/UCSC/dm3/Sequence/tailorIndex/dm3
# bowtie index
BOWTIE_INDEX=/home/hanb/nearline/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome
# bwa index
BWA_INDEX=/home/hanb/nearline/Drosophila_melanogaster/UCSC/dm3/Sequence/BWAIndex/genome.fa

FQ=Phil.hen1Hets.trimmed.p10.fq 
for CPU in 2 4 8 12 24; do 
	time bowtie -p $CPU -S -v 0 --best --strata -a $BOWTIE_INDEX $FQ > ${FQ%fq}b1.sam 2>${FQ%fq}b1.log
	time ( bwa aln -t $CPU -o0 -R1000000 $BWA_INDEX $FQ > ${FQ%fq}sai && bwa samse $BWA_INDEX ${FQ%fq}sai $FQ > ${FQ%fq}bwa.sam ) 
	time tailor map -i $FQ -p $TAILOR_INDEX -n $CPU -o ${FQ%fq}t.sam
done

FQ=Phil.hen1Muts.trimmed.p10.fq
for CPU in 2 4 8 12 24; do 
	time bowtie -p $CPU -S -v 0 --best --strata -a $BOWTIE_INDEX $FQ > ${FQ%fq}b1.sam 2>${FQ%fq}b1.log
	time ( bwa aln -t $CPU -o0 -R1000000 $BWA_INDEX $FQ > ${FQ%fq}sai && bwa samse $BWA_INDEX ${FQ%fq}sai $FQ > ${FQ%fq}bwa.sam ) 
	time tailor map -i $FQ -p $TAILOR_INDEX -n $CPU -o ${FQ%fq}t.sam
done
