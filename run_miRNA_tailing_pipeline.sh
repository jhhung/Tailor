#! /bin/bash -x

# This script is small RNA tailing analysis pipeline 
# associated with the Tailor software.

################
# Major Config #
################
VERSION="1.0.0"
TAILOR_PATH=
# function to check whether current shell can find all the software/programs needed to finish running the pipeline
function checkExist {
	echo -ne "\e[1;32m\"${1}\" is using: \e[0m" && which "$1"
	[[ $? != 0 ]] && echo -e "\e[1;36mError: cannot find software/function ${1}! Please make sure that you have installed the pipeline correctly.\nExiting...\e[0m" && \
	exit 1
}

echo -e "\e[1;35mTesting required softwares/scripts:\e[0m"
checkExist "echo" # used to print message
checkExist "tee" # used to store to file while keeping piping
checkExist "date" # used to get time
checkExist "basename" # used to get basename
checkExist "dirname" # used to get diretory name
checkExist "mkdir" # used to creat new directory
checkExist "readlink" # used to find the real location of a file
checkExist "rm" # used to delete file/directory
checkExist "mv" # used to move or rename file
checkExist "sort" # used to sort alphnum
checkExist "touch" # used to create new empty file
checkExist "md5sum" # used to create a uid so that the status file for different run will not mess up eath other
checkExist "awk" # a programming language
checkExist "samtools" # tools to process sam/bam file
checkExist "bedtools" # tools to process bed formatted file
checkExist "tailor" # tools to perform mapping
echo -e "\e[1;35mDone with testing required softwares/scripts, starting pipeline...\e[0m"


#########
# USAGE #
#########
# usage function
usage() {
echo -en "\e[1;36m"
cat << EOF

usage: $0 -i input_file.fq -p hairpin.fa -m mature.fa -o output_directory[current directory] -c cpu[8] 

This is a small RNA tailing pipeline associated with the software Tailor (git@github.com:jhhung/Tailor.git).

OPTIONS:
	-h      Show this message
	-i      Input file in fastq
	-m      Mature miRNA fasta file from miRBase
	-p      Hairpin fasta file from miRBase
	-o      Output directory, default: current directory
	-c      Number of CPUs to use, default: 8

EOF
echo -en "\e[0m"
}

# taking options
while getopts "hi:p:o:c:m:" OPTION
do
	case $OPTION in
		h)
			usage && exit 1
		;;
		i)
			INPUT_FQ=`readlink -f $OPTARG`
		;;
		o)
			OUTDIR=$OPTARG
		;;
		c)
			CPU=$OPTARG
		;;
		p)
			HAIRPIN_FA=`readlink -f $OPTARG`
		;;
		m)
			MATURE_FA=`readlink -f $OPTARG`
		;;
		?)
			usage && exit 1
		;;
	esac
done

# if INPUT_FQ or HAIRPIN_FA is undefined, print out usage and exit
if [[ -z $INPUT_FQ ]] || [[ -z $HAIRPIN_FA ]] || [[ -z $MATURE_FA ]]
then
	usage && exit 1
fi

# check file status
[ ! -s $INPUT_FQ ]  && echo -e "\e[1;31mError: cannot open $INPUT_FQ \e[0m" && exit 1
[ ! -s $HAIRPIN_FA ] && echo -e "\e[1;31mError: cannot open $HAIRPIN_FA\e[0m" && exit 1
[ ! -s $MATURE_FA ] && echo -e "\e[1;31mError: cannot open $MATURE_FA\e[0m" && exit 1 

# if CPU is undefined or containing non-numeric char, then use 8
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8

# OUTPUT directory defined? otherwise use current directory
[ ! -z $OUTDIR ] || OUTDIR=$PWD

# test wether we need to create the new directory
mkdir -p "${OUTDIR}"
if [[ $? == 1 ]]; then echo -e "\e[1;31mWarning: Cannot create directory ${OUTDIR}. Using the current direcory\e[0m" && OUTDIR=$PWD ;fi  

# enter destination direcotry
cd ${OUTDIR}
if [[ $? == 1 ]]; then echo -e "\e[1;31mError: Cannot access directory ${OUTDIR}... Exiting...\e[0m" && exit 1 ;fi

# test writtability
touch .writting_permission && rm -rf .writting_permission
if [[ $? == 1 ]]; then echo -e "\e[1;31mError: Cannot write in directory ${OUTDIR}... Exiting...\e[0m"  && exit 1 ;fi

#############
# Variables #
#############
# Prefix
PREFIX=`basename $INPUT_FQ` && PREFIX=${PREFIX%.f[qa]*}

# index configuration
INDEX=${HAIRPIN_FA%.fa*}
INDEX_PREFIX=`basename $HAIRPIN_FA`
INDEX_FOLDER=`dirname $HAIRPIN_FA`
MATURE_BASENAME=`basename $MATURE_FA`
# test if the user has privilege to write in the hairpin.fa folder
touch ${INDEX_FOLDER}/.writting_permission && rm -rf ${INDEX_FOLDER}/.writting_permission || INDEX=$PWD/$INDEX_PREFIX

# step counter
STEP=1

# unique job id
JOBUID=`echo ${PREFIX} | md5sum | cut -d" " -f1`

# formating date for log
ISO_8601='%Y-%m-%d %H:%M:%S %Z'

# table to store the basic statistics of the library (genomic mappability). PS: counts of each genomic features will be stored in .summary file, produced by intersect_all.sh 
TABLE=${PREFIX}.stats

# log file to store the timing
LOG=${PREFIX}.log

##############################
# beginning running pipeline #
##############################
# if log file already exists, then "resume" it. otherwise "begin" it.
[ ! -f $LOG ] && \
echo -e "`date "+$ISO_8601"`\tbeginning Tailor pipeline for small RNA version $VERSION"   | tee -a $LOG || \
echo -e "`date "+$ISO_8601"`\tresuming Tailor pipeline for small RNA version $VERSION" | tee -a $LOG

# build index, if an index has already exist, tailor will not rewrite it, until using -f option
echo -e "`date "+$ISO_8601"`\tbuilding hairpin index if not exist" | tee -a $LOG
[ ! -f .${JOBUID}.status.build_hairpin_index ] && \
	tailor build -i $HAIRPIN_FA -p $INDEX && \
	touch .${JOBUID}.status.build_hairpin_index
STEP=$((STEP+1))

# map annotated mature miRNA sequence to hairpin index and determine the annotated coordinate
echo -e "`date "+$ISO_8601"`\tmap annoated mature miRNA to hairpin to determine the annoated coordiate" | tee -a $LOG
[ ! -f .${JOBUID}.status.find_annotated_coordinates ] && \
	awk '{printf "@%s\n", substr($1,2); getline; l=length($1); printf "%s\n+\n", $1; for (i=1;i<=l;++i) printf "%s","I"; printf "\n" }' $MATURE_FA > ${MATURE_BASENAME%.fa*}.fq && \
	tailor map -i ${MATURE_BASENAME%.fa*}.fq -p $INDEX -n $CPU | \
	samtools view -bS - | \
	bedtools bamtobed -i - > \
	${PREFIX}.annotated_coordinates.bed && \
	rm -rf ${MATURE_BASENAME%.fa*}.fq && \
	touch .${JOBUID}.status.find_annotated_coordinates
STEP=$((STEP+1))

# perform the mapping
echo -e "`date "+$ISO_8601"`\tmap small RNA reads to the index" | tee -a $LOG
[ ! -f .${JOBUID}.status.tailor_mapping ] && \
	tailor map -i $INPUT_FQ -p $INDEX -n $CPU -o ${PREFIX}.tailor.sam && \
	touch .${JOBUID}.status.tailor_mapping
STEP=$((STEP+1))

# separating perfect match and tailing match and count their length distribution
echo -e "`date "+$ISO_8601"`\tseparating perfect mapping and prefix mapping and generate bed" | tee -a $LOG
[ ! -f .${JOBUID}.status.separate_perfect_and_tail ] && \
	awk 'BEGIN \
	{ \
		FS="\t"; OFS="\t"; \
		getline; \
		while (substr($1,1,1)=="@") \
		{ \
			print $0 >> "/dev/stdout"; print $0 >> "/dev/stderr"; \
			getline; \
		} \
		if ($5<255) print $0 >> "/dev/stderr"; \
		else print $0 >> "/dev/stdout"; \
	} \
	{ \
		if ($5<255) print $0 >> "/dev/stderr"; \
		else print $0 >> "/dev/stdout"; \
	}' ${PREFIX}.tailor.sam \
	1> ${PREFIX}.tailor.perfectMatch.sam \
	2> ${PREFIX}.tailor.tailedMatch.sam && \
	samtools view -bS ${PREFIX}.tailor.perfectMatch.sam | tee ${PREFIX}.tailor.perfectMatch.bam | bedtools bamtobed -i - | tee ${PREFIX}.tailor.perfectMatch.bed | awk 'BEGIN{minLen=18; maxLen=30}{l=$3-$2; if (l<minLen) minLen=l; if (l>maxLen) maxLen=l; ++ct[l]}END{for (i=minLen;i<=maxLen;++i){printf "%d\t%d\n", i, ct[i]?ct[i]:0}}' > ${PREFIX}.tailor.perfectMatch.lenDis && \
	samtools view -bS ${PREFIX}.tailor.tailedMatch.sam  | tee ${PREFIX}.tailor.tailedMatch.bam  | bedtools bamtobed -i - | tee ${PREFIX}.tailor.tailedMatch.bed  | awk 'BEGIN{minLen=18; maxLen=30}{l=$3-$2; if (l<minLen) minLen=l; if (l>maxLen) maxLen=l; ++ct[l]}END{for (i=minLen;i<=maxLen;++i){printf "%d\t%d\n", i, ct[i]?ct[i]:0}}' > ${PREFIX}.tailor.tailedMatch.lenDis
	touch .${JOBUID}.status.separate_perfect_and_tail
STEP=$((STEP+1))

# rename the reads based on the mature name;adjust the 5' and 3' end positions
echo -e "`date "+$ISO_8601"`\treannoated the bed file based on miRBase annoated mature name" | tee -a $LOG
[ ! -f .${JOBUID}.status.reannotate ] && \
	cat ${PREFIX}.tailor.perfectMatch.bed ${PREFIX}.tailor.tailedMatch.bed | \
	bedtools intersect -wo -a - -b ${PREFIX}.annotated_coordinates.bed -f 0.75 | \
	awk 'BEGIN{FS="\t";OFS="\t"}{print $10,$2-$8,$3-$9,$4,$5,$6}' \
	> ${PREFIX}.reannoated.bed && \
	touch .${JOBUID}.status.reannotate
STEP=$((STEP+1))

# getting statistics on each mature miRNA
echo -e "`date "+$ISO_8601"`\tgetting statistics" | tee -a $LOG
[ ! -f .${JOBUID}.status.statistics ] && \
	echo -ne "miRNA\ttotal counts\tperfect match counts\tprefix match counts\n" > ${PREFIX}.count && \
	awk '{if ($5==255) ++pm[$1]; else ++pf[$1]; ++all[$1]; }END{for (mir in all) {printf "%s\t%d\t%d\t%d\n", mir, all[mir], pm[mir]?pm[mir]:0, pf[mir]?pf[mir]:0}}' ${PREFIX}.reannoated.bed >> ${PREFIX}.count && \
	touch .${JOBUID}.status.statistics
STEP=$((STEP+1))

################
# draw figures #
################
checkExist "Rscript" # R, a programing language
checkExist "draw_lendis.R" # drawing reads
draw_lendis.R ${PREFIX}.tailor.perfectMatch.lenDis ${PREFIX}.tailor.tailedMatch.lenDis $PREFIX





