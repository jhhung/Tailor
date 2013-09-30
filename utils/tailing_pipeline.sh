#!/bin/bash -x

################
# Major Config #
################
# using small RNA pipeline
export PIPELINE_DIRECTORY=/home/hanb/nearline/small_RNA_Pipeline
# set PATH to be aware of pipeline/bin; this bin directory needs to be searched first
export PATH=${PIPELINE_DIRECTORY}/bin:$PATH
# index location
TAILOR_INDEX_DM3=/home/hanb/nearline/Drosophila_melanogaster/UCSC/dm3/Sequence/tailorIndex/dm3
TAILOR_INDEX_MM9=/home/hanb/nearline/Mus_musculus/UCSC/mm9/Sequence/tailorIndex/mm9

#########
# USAGE #
#########
# usage function
usage() {
echo -en "\e[1;36m"
cat << EOF

usage: $0 -i input_file.fq -s fly|mouse -o output_directory[current directory] -c cpu[8] 

This is small RNA pipeline specificly for discovering tailing events
by using deep sequencing data.

Please email Bo.Han@umassmed.edu for any questions or bugs. 
Thanks for using it. 

OPTIONS:
	-h      Show this message
	-i      Input file in fastq format, with full directory
	-o      Output directory, default: current directory
	-s      Species [  mouse | fly ]
	-c      Number of CPUs to use, default: 8

EOF
echo -en "\e[0m"
}

# taking options
while getopts "hi:c:o:s:" OPTION
do
	case $OPTION in
		h)
			usage && exit 1
		;;
		i)
			FQ=$OPTARG
		;;
		o)
			OUTDIR=$OPTARG
		;;
		c)
			CPU=$OPTARG
		;;
		s)
			ORGANISM=$OPTARG
		;;
		?)
			usage && exit 1
		;;
	esac
done

# if FQ or ORGANISM is undefined, print out usage and exit
if [[ -z $FQ ]] || [[ -z $ORGANISM ]] 
then
	usage && exit 1
fi

#######################
# Function definition #
#######################
# count reads for bed2 format
function bedwc {
	awk '{a[$7]=$4}END{COUNTER=0; for (b in a){COUNTER+=a[b]} printf "%d" , COUNTER}' $1
}
export -f bedwc 
# produce length distribution for bed2 format
function bed2lendis {
	awk '{a[$7]=$4}END{m=0; for (b in a){c[length(b)]+=a[b]; if (length(b)>m) m=length(b)} for (d=1;d<=m;++d) {print d"\t"(c[d]?c[d]:0)}}'  $1 | sort -k1,1n 
}
export -f bed2lendis

##################
# File/Dir check #
##################
# check whether input file exist
[ ! -f $FQ ] && echo -e "\e[1;31mError: Cannot find input file $FQ\e[0m" && exit 1
# if CPU is undefined or containing non-numeric char, then use 8
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
# get full path of the FQ and use the full path when getting the insert file.
FULLPATH_FQ=`readlink -f $FQ`
# use $FQ later only as filename
FQ=${FQ##*/}
# insert file
INSERT=${FQ%.f[qa]*}.insert
# log file name
LOG=${FQ%.f[qa]*}.log
# OUTPUT directory defined? otherwise use current directory
[ ! -z $OUTDIR ] || OUTDIR=$PWD
# test wether we need to create the new directory
mkdir -p "${OUTDIR}" || echo -e "\e[1;31mWarning: Cannot create directory ${OUTDIR}. Using the direcory of input fastq file\e[0m"
# enter destination direcotry, since later it will be ran in this directory, there is NO need to add OUTPUT DIR in the following codes
cd ${OUTDIR} || echo -e "\e[1;31mError: Cannot access directory ${OUTDIR}... Exiting...\e[0m" || exit 1
# test writtability
touch ${OUTDIR}/.writting_permission && rm -rf ${OUTDIR}/.writting_permission || echo -e "\e[1;31mError: Cannot write in directory ${OUTDIR}... Exiting...\e[0m" || exit 1

#############
# Variables #
#############
# step counter
STEP=1
# unique job id
JOBUID=`echo ${FQ##*/} | md5sum | cut -d" " -f1`
# formating date for log
ISO_8601='%Y-%m-%d %H:%M:%S %Z'

case ${ORGANISM} in
mouse)
	TAILOR_INDEX=$TAILOR_INDEX_MM9
	GENOME_FA=/home/hanb/nearline/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa
;;
fly)
	TAILOR_INDEX=$TAILOR_INDEX_DM3
	GENOME_FA=/home/hanb/nearline/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa
;;
*)
	echo -e "\e[1;31mError: Unknown orgnanism... Currently only mouse/fly is supported\e[0m"
	exit 2
;;
esac

##############################
# beginning running pipeline #
##############################
# if log file already exists, then "resume" it. otherwise "begin" it.
[ ! -f $LOG ] && \
echo -e "`date "+$ISO_8601"`\tbeginning running Zamore Lab small RNA pipeline version $smallRNA_Pipeline_Version in one lib mode"   | tee -a $LOG || \
echo -e "`date "+$ISO_8601"`\t...resuming running Zamore Lab small RNA pipeline version $smallRNA_Pipeline_Version in one lib mode" | tee -a $LOG 

# converting fastq to insert format
echo -e "`date "+$ISO_8601"`\tconverting fastq format into insert format" | tee -a $LOG
[ ! -f .status.${STEP}.fq2insert ] && \
	fastq2insert ${FULLPATH_FQ} ${INSERT} && \
	touch .status.${STEP}.fq2insert
STEP=$((STEP+1))

# converting insert to naive fastq
echo -e "`date "+$ISO_8601"`\tconverting insert format into naive fastq format" | tee -a $LOG
[ ! -f .status.${STEP}.insert2fastq ] && \
	awk '{qual=""; for (i=1;i<=length($1);++i) qual=qual"I"; printf "@%s_%s\n%s\n+\n%s\n", $1,$2,$1,qual}' ${INSERT} > ${FQ}2 && \
	touch .status.${STEP}.insert2fastq
STEP=$((STEP+1))

# tailor mapping
echo -e "`date "+$ISO_8601"`\ttailor mapping" | tee -a $LOG
[ ! -f .status.${STEP}.tailor_mapping ] && \
	tailor map -i ${FQ}2 -p $TAILOR_INDEX -n $CPU > ${FQ%.f[qa]*}.tailor.sam && \
	samtools view -bS ${FQ%.f[qa]*}.tailor.sam | tee ${FQ%.f[qa]*}.tailor.bam | bedtools bamtobed -i - > ${FQ%.f[qa]*}.tailor.bed && \
touch .status.${STEP}.tailor_mapping
STEP=$((STEP+1))

# converting bed to bed2
echo -e "`date "+$ISO_8601"`\tconverting bed to bed2" | tee -a $LOG
[ ! -f .status.${STEP}.converting2 ] && \
	awk 'BEGIN{OFS="\t"}{if (ARGIND==1) {++ct[$4]} else {split ($4,arr,"_"); print $1,$2,$3,arr[2],ct[$4],$6,arr[1],$5}}' ${FQ%.f[qa]*}.tailor.bed ${FQ%.f[qa]*}.tailor.bed > ${FQ%.f[qa]*}.tailor.bed3 && \
touch .status.${STEP}.converting2
STEP=$((STEP+1))

# separating unique and multip 
UNIQUE_BED=${FQ%.f[qa]*}.tailor.unique.bed3
MULTIP_BED=${FQ%.f[qa]*}.tailor.multip.bed3
echo -e "`date "+$ISO_8601"`\tseparating unique and multip" | tee -a $LOG
[ ! -f .status.${STEP}.sep_unique_multip ] && \
	awk 'BEGIN{FS="\t"; OFS="\t"}{if ($5==1	) print $0 >> "/dev/stdout"; else print $0 >> "/dev/stderr"}' ${FQ%.f[qa]*}.tailor.bed3
		1 > $UNIQUE_BED \
		2 > $MULTIP_BED && \
touch .status.${STEP}.sep_unique_multip
STEP=$((STEP+1))

# intersecting to genomic features
#@ using small RNA pipeline
echo -e "`date "+$ISO_8601"`\tdoing intersecting" | tee -a $LOG
[ ! -f .status.${STEP}.intersecting ] && \
	intersect_all_tailor.sh \
	$UNIQUE_BED \
	$MULTIP_BED \
	${FQ%.f[qa]*}.tailor.tailing.summary \
	$OUT \
	$ORGANISM \
	$GENOME_FA \
	$CPU && \
touch .status.${STEP}.intersecting
STEP=$((STEP+1))






