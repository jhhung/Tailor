#! /bin/bash -x

# This script is associated with the Tailor software and 
# only used for internal testing purpose and to repeat the 
# result in the published paper. 

# map zebrafish

# tailor index
TAILOR_INDEX=/home/hanb/nearline/Danio_rerio/UCSC/danRer7/Sequence/tailorIndex/danRer7
# bowtie index
BOWTIE_INDEX=/home/hanb/nearline/Danio_rerio/UCSC/danRer7/Sequence/BowtieIndex/genome
# bwa index
BWA_INDEX=/home/hanb/nearline/Danio_rerio/UCSC/danRer7/Sequence/BWAIndex/genome.fa

FQ=Ketting.hen1Hets.trimmed.p10.fq 
for CPU in 2 4 8 12 24; do 
	time bowtie -p $CPU -S -v 0 --best --strata -a $BOWTIE_INDEX $FQ > ${FQ%fq}b1.sam 2>${FQ%fq}b1.log
	time ( bwa aln -t $CPU -o0 -R1000000 $BWA_INDEX $FQ > ${FQ%fq}sai && bwa samse $BWA_INDEX ${FQ%fq}sai $FQ > ${FQ%fq}bwa.sam ) 
	time tailor map -i $FQ -p $TAILOR_INDEX -n $CPU -o ${FQ%fq}t.sam
done

FQ=Ketting.hen1Muts.trimmed.p10.fq
for CPU in 2 4 8 12 24; do 
	time bowtie -p $CPU -S -v 0 --best --strata -a $BOWTIE_INDEX $FQ > ${FQ%fq}b1.sam 2>${FQ%fq}b1.log
	time ( bwa aln -t $CPU -o0 -R1000000 $BWA_INDEX $FQ > ${FQ%fq}sai && bwa samse $BWA_INDEX ${FQ%fq}sai $FQ > ${FQ%fq}bwa.sam ) 
	time tailor map -i $FQ -p $TAILOR_INDEX -n $CPU -o ${FQ%fq}t.sam
done

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
