#!/bin/bash -x

##########
# Config #
##########
export PIPELINE_DIRECTORY=$(dirname `readlink -f $0`)
export PATH=${PIPELINE_DIRECTORY}/bin:$PATH
export TAILOR_INDEX=$PIPELINE_DIRECTORY/indexes
export VERSION=1.0.0
#########
# USAGE #
#########
# usage function
usage() {
cat << EOF

+++++++++++++++++++++++++
+Tailor general pipeline+
+++++++++++++++++++++++++

A generalized pipeline to analyze tailing events from Next Generation Sequencing.

usage: $0 \ 
	-i input_file.fq \ 
	-g dm3.fa \ 
	-o output_directory[current directory] \ 
	-c cpu[8] 

OPTIONS:
<required>
	-h      Show this message
	-i      Input file in fastq format; it is highly recommend to filter reads by Phred score
	-g      Genome fasta file. The pipeline will automatically generate index with the prefix if it does not already exist
<optional>
	-o      Output directory, default <current working directory>
	-c      Number of CPUs to use, default <8>
	-q      The Phred score used to filter the reads. One read needs to have all its bases no less than this value to pass, default <20>
	
EOF
}

while getopts "hi:g:o:c:q:" OPTION
do
	case $OPTION in
		h)	usage && exit 0 ;;
		i)	export	INPUT_FQ=`readlink -f $OPTARG` ;;
		g)	INDEX_FA=`readlink -f $OPTARG` 
			INDEX_FA_NAME=`basename $INDEX_FA`
			INDEX_PREFIX=${INDEX_FA_NAME%.fa*}
			INDEX=$TAILOR_INDEX/$INDEX_PREFIX
		;;
		o)	export	OUTDIR=`readlink -f $OPTARG` ;;
		c)	export	CPU=$OPTARG ;;
		q)	MIN_PHRED=$OPTARG ;;
		?)	usage && exit 1 ;;
	esac
done

#######################
# Function definition #
#######################
function echo2 {
COLOR_GREEN="\e[32;40m";
COLOR_RED_BOLD="\e[31;1m"; 
COLOR_MAGENTA="\e[35;40m"; 
COLOR_END="\e[0m"; 	
ISO_8601='%Y-%m-%d %H:%M:%S %Z'
case $2 in 
	error)		echo -e $COLOR_RED_BOLD"[`date "+$ISO_8601"`] Error: $1${COLOR_END}" && exit 1 ;;
	warning)	echo -e $COLOR_MAGENTA"[`date "+$ISO_8601"`] Warning: $1${COLOR_END}" ;;
	*)			echo -e $COLOR_GREEN"[`date "+$ISO_8601"`] $1${COLOR_END}";;
esac
}
export -f echo2

function bedwc {
	awk '{a[$7]=$4}END{COUNTER=0; for (b in a){COUNTER+=a[b]} printf "%d" , COUNTER}' $1
}
export -f bedwc 

function bed2lendis {
	awk '{a[$7]=$4}END{m=0; for (b in a){c[length(b)]+=a[b]; if (length(b)>m) m=length(b)} for (d=1;d<=m;++d) {print d"\t"(c[d]?c[d]:0)}}'  $1 | sort -k1,1n 
}
export -f bed2lendis

##################
# File/Dir check #
##################
[ -z $INPUT_FQ ] && echo2 "missing -i for input fastq file" "error"
[ -z $INDEX_FA ] && echo2 "missing -g for reference fasta file" "error"
[ ! -f $INPUT_FQ ] && echo2 "cannot file $INPUT_FQ" "error"
[ -z "${CPU##*[!0-9]*}" ] && export CPU=8 && echo2 "using 8 CPUs" "warning"
[ -z "${MIN_PHRED}" ] && export MIN_PHRED=20 && echo2 "using 20 as minimal phred score allowed" "warning"
[ -z $OUTDIR ] && export OUTDIR=$PWD
mkdir -p "${OUTDIR}" || echo2 "Cannot create directory ${OUTDIR}. Using the direcory of input fastq file" "warning"
cd ${OUTDIR} || echo2 "Cannot access directory ${OUTDIR}..." "error"
touch .writting_permission  || echo2 "Cannot write in directory ${OUTDIR}..." "error"
rm -rf .writting_permission

#############
# Variables #
#############
STEP=1
FQ=`basename $INPUT_FQ`
export PREFIX=${FQ%.f[aq]*}
JOBUID=`echo ${FQ##*/} | md5sum | cut -d" " -f1`
INSERT=${FQ%.f[qa]*}.insert

##########
# folder #
##########
MAPPING_DIR=mapping && mkdir -p $MAPPING_DIR
FEATURES_DIR=genomic_feature && mkdir -p $FEATURES_DIR

##############################
# beginning running pipeline #
##############################
echo2 "Begin running tailing pipeline version $VERSION"

echo2 "Checking phred version" 
PHRED_SCORE=`perl $PIPELINE_DIRECTORY/bin/SolexaQA_piper.pl ${INPUT_FQ}`
case ${PHRED_SCORE} in
solexa)		OFFSET=59 ;; # Solexa+64, raw reads typically (-5, 40)
illumina)	OFFSET=64 ;; # Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
sanger)		OFFSET=33 ;; # Phred+33,  raw reads typically (0, 40) (http://en.wikipedia.org/wiki/FASTQ_format)
*)			echo2 "unable to determine the fastq version. Using sanger..." "warning"
			OFFSET=33 ;;
esac

echo2 "Filtering the fastq by the phred score and pool reads with same sequence"
[ ! -f .status.${STEP}.phred_filter_and_pool ] && \
phred_filter_and_pool \
	-i	$INPUT_FQ \
	-q	$MIN_PHRED \
	-s	$OFFSET \
	-o	${PREFIX}.p${MIN_PHRED}.fq \
	2> ${PREFIX}.p${MIN_PHRED}.phred_filter.log
touch .status.${STEP}.phred_filter_and_pool
STEP=$((STEP+1))

echo2 "Building the index if not exist"
tailor build -i $INDEX_FA -p $INDEX

echo2 "Mapping the input fastq to the genome reference" 
[ ! -f .status.${STEP}.tailor_mapping ] && \
	tailor map \
		-i ${PREFIX}.p${MIN_PHRED}.fq \
		-p $INDEX \
		-n $CPU \
		2> $MAPPING_DIR/${PREFIX}.tailor.log | \
	tailor_sam_to_bed \
	> $MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.bed2 && \
touch .status.${STEP}.tailor_mapping
STEP=$((STEP+1))

# # intersecting to genomic features
# echo2 "Assigning reads to different genomic structures"
# [ ! -f .status.${STEP}.intersecting ] && \
# 	intersect_to_genomic_feature_tailor.sh \
# 		$MAPPING_DIR/${PREFIX}.tailor.bed2 && \
# 	touch .status.${STEP}.intersecting
# STEP=$((STEP+1))






