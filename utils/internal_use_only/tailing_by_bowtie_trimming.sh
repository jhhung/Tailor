#! /bin/bash -x

# This script is associated with the Tailor software and 
# only used for internal testing purpose and to repeat the 
# result in the published paper. 

# This is a scripts implemented to find untemplated addition to the 3' end of small RNA using iteration between bowtie and trimming

INPUT_FQ=${1}
BOWTIE_INDEX=/share/biocore/genome/bowtie/hg19
SAM=${1%.f[aq]*}.hg19.bowtie.sam
SAM_HEADER=""
# keep iterating until the no reads left in 
TIME_RUN=0
UN_MAPPED_FQ=${INPUT_FQ}"."${TIME_RUN}
while [ -s $INPUT_FQ ]
do
	bowtie -S $SAM_HEADER -p 8 -v 0 -a -3 $TIME_RUN --un $UN_MAPPED_FQ $BOWTIE_INDEX $INPUT_FQ >> $SAM && \
	TIME_RUN=$((TIME_RUN+1)) && \
	SAM_HEADER="--sam-nohead" && \
	INPUT_FQ=$UN_MAPPED_FQ && \
	UN_MAPPED_FQ=${1}"."${TIME_RUN}
done
samtools view -bS $SAM > ${SAM%sam}bam
