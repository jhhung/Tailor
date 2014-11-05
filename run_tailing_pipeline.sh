#!/bin/bash

# Tailor, a BWT-based aligner for non-templated RNA tailing
# Copyright (C) 2014 Min-Te Chou, Bo W Han, Jui-Hung Hung
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

##########
# Config #
##########
export PIPELINE_DIRECTORY=$(dirname `readlink -f $0`)
export PATH=${PIPELINE_DIRECTORY}/bin:$PATH
export TAILOR_INDEX=$PIPELINE_DIRECTORY/indexes
export VERSION=1.1.0

#########
# USAGE #
#########
# usage function
usage() {
cat << EOF

+++++++++++++++++++
+ Tailor pipeline +
+++++++++++++++++++

Tailor pipeline to analyze tailing events from Next Generation Sequencing.
It requires 
1. the input reads in FastQ format;
2. genome fasta sequence in FastA format; it will generate Tailor index in a sub-directory of this pipeline;
   The index will be named as the md5 value of the ABSOLUTE path of the fasta file; 
   Once the index for a specific Fasta file has been built, it will not be built next time you use it (unless you rename the file or move it);
   To avoid collision, for a new genome, only run one job at a time.
3. genomic_feature_file; this file stores information on the name of different genomic structures. please refer to the annotation/dm3.genomic_features for the format;
   We have included the feature files for several genomes in the annotation directory.
4. [ optional ] microRNA hairpin sequence in FastA file; We provide a script to obtain mature and hairpin sequence from miRbase: obtain_miRNA.sh
5. [ optional ] microRNA mature sequence in FastA format; We provide a script to obtain mature and hairpin sequence from miRbase: obtain_miRNA.sh
6. [ optional ] output directory
7. [ optional ] number of threads to use

usage: $0 \ 
	-i input_file.fq \ 
	-g dm3.fa \ 
	-t genomic_feature_file [ annotation/dm3.genomic_features ] \ 
	-H hairpin.fa [ optional ] \ 
	-M mature.fa [ optional ]  \ 
	-o output_directory [ current directory ] \ 
	-c cpu[ 8 ] 

OPTIONS:
<required>
	-h      Show this message
	-i      Input file in fastq format; already has its adaptor and baracode removed
	-g      Genome fasta file. The pipeline will automatically generate index with the prefix if it does not already exist
	-t      Files to store the genomic features. See annotation folder for examples.
<optional>
	-H      microRNA hairpin sequence in fasta format. The pipeline will automatically generate index.
	-M      microRNA mature sequence in fasta format. It is used to annotation the 5' and 3' ends of miRNA.
	-o      Output directory, default <current working directory>
	-c      Number of CPUs to use, default <8>
	-q      The Phred score used to filter the reads. One read needs to have all its bases no less than this value to pass, default <20>
	-v      Allow mismatches in mapping
EOF
}

while getopts "hi:g:o:c:q:t:H:M:v" OPTION
do
	case $OPTION in
		h)	usage && exit 0 ;;
		i)	export INPUT_FQ=`readlink -f $OPTARG` ;;
		g)	INDEX_FA=`readlink -f $OPTARG`
			INDEX_UID=`echo ${INDEX_FA} | md5sum | cut -d" " -f1`
			INDEX_FA_NAME=`basename $INDEX_FA`
			INDEX_PREFIX=${INDEX_FA_NAME%.fa*}
			INDEX=$TAILOR_INDEX/$INDEX_UID
		;;
		t)	export GENOMIC_FEATURE_FILE=$OPTARG ;;
		H)	HAIRPIN_INDEX_FA=`readlink -f $OPTARG` 
			HAIRPIN_INDEX_FA_NAME=`basename $HAIRPIN_INDEX_FA`
			HAIRPIN_INDEX_PREFIX=${HAIRPIN_INDEX_FA_NAME%.fa*}
			HAIRPIN_INDEX=$TAILOR_INDEX/$HAIRPIN_INDEX_PREFIX
		;;
		M)	MATURE_FA=`readlink -f $OPTARG` ;;
		o)	export OUTDIR=`readlink -f $OPTARG` ;;
		c)	export CPU=$OPTARG ;;
		v)	ALLOW_MISMATCH="-v 1" ;;
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
[ ! -s $INPUT_FQ ] && echo2 "cannot file file $INPUT_FQ" "error"
[ -z "${CPU##*[!0-9]*}" ] && export CPU=8 && echo2 "using 8 CPUs" "warning"
[ -z "${MIN_PHRED}" ] && export MIN_PHRED=20 && echo2 "using 20 as minimal phred score allowed" "warning"
[ -z "$GENOMIC_FEATURE_FILE" ] && \
	echo2 "-t option unspecified, using ${PIPELINE_DIRECTORY}/annotation/dm3.genomic_features" "warning" && \
	export GENOMIC_FEATURE_FILE=${PIPELINE_DIRECTORY}/annotation/dm3.genomic_features

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
export PDF_DIR=pdfs && mkdir -p $PDF_DIR

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
[ ! -f .${JOBUID}.status.${STEP}.phred_filter_and_pool ] && \
phred_filter_and_pool \
	-i	$INPUT_FQ \
	-q	$MIN_PHRED \
	-s	$OFFSET \
	-o	${PREFIX}.p${MIN_PHRED}.fq \
	2> ${PREFIX}.p${MIN_PHRED}.phred_filter.log
touch .${JOBUID}.status.${STEP}.phred_filter_and_pool
STEP=$((STEP+1))

echo2 "Building the index if not exist"
tailor build -i $INDEX_FA -p $INDEX # will raise warning if the index exists

echo2 "Mapping the input fastq to the genome reference" 
[ ! -f .${JOBUID}.status.${STEP}.tailor_mapping ] && \
	tailor map \
		-i ${PREFIX}.p${MIN_PHRED}.fq \
		-p $INDEX \
		-n $CPU \
		$ALLOW_MISMATCH \
		2> $MAPPING_DIR/${PREFIX}.tailor.log | \
	tee $MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.sam | \
	tailor_sam_to_bed \
	> $MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.bed2 && \
touch .${JOBUID}.status.${STEP}.tailor_mapping
STEP=$((STEP+1))

# bed2 format specification:
# f1-6: same as ordinary bed
# f7: original sequence as in input
# f8: sequence of the tail (* if no tail)
# f9: length of tail 
# f10: MAPQ string for mismatches (* if not mismatch)
# tailor_bed2_counter.py does not consider f10 (mismatches)
echo2 "Draw overall length distribution with tailing information"
[ ! -f .${JOBUID}.status.${STEP}.draw_overall_lendis ] && \
	python $PIPELINE_DIRECTORY/bin/tailor_bed2_counter.py \
		$MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.bed2 \
		1> $MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.sum \
		2> $MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.single_nt_sum && \
	Rscript --slave $PIPELINE_DIRECTORY/bin/draw_tailor_lendis.R \
		$MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.single_nt_sum \
		$PDF_DIR/${PREFIX}.p${MIN_PHRED}.pdf \
		${PREFIX} \
		"total" && \
touch .${JOBUID}.status.${STEP}.draw_overall_lendis
STEP=$((STEP+1))

echo2 "Assigning reads to different genomic structures"
[ ! -f .${JOBUID}.status.${STEP}.intersecting ] && \
	bash $DEBUG intersect_to_genomic_feature_tailor.sh \
		$MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.bed2 && \
	touch .${JOBUID}.status.${STEP}.intersecting
STEP=$((STEP+1))

#############################
# microRNA tailing analysis #
#############################
if [ ! -z $HAIRPIN_INDEX_FA ]; then

	[ -z $MATURE_FA ]	&& echo2 "missing -M for mirBase mature miRNA in fasta format" "error"
	[ ! -s $MATURE_FA ]	&& echo2 "cannot find file $MATURE_FA" "error"
	[ ! -s $HAIRPIN_INDEX_FA ] && echo2 "cannot fine file $HAIRPIN_INDEX_FA" "error"
	
	ANNOTATION_DIR=annotate_mature_miRNA && mkdir -p $ANNOTATION_DIR
	MAPPING_DIR=hairpin_mapping && mkdir -p $MAPPING_DIR
	BALLOON_DIR=balloon_plot_individual && mkdir -p $BALLOON_DIR
	
	echo2 "Building index, if not already exist"
	tailor build -i $HAIRPIN_INDEX_FA -p $HAIRPIN_INDEX 2>/dev/null
	
	echo2 "Map annoated mature miRNA to hairpin to determine the annoated coordiate"
	[ ! -f .${JOBUID}.status.${STEP}.find_annotated_coordinates ] && \
		awk '{printf "@%s\n", substr($1,2); getline; l=length($1); printf "%s\n+\n", $1; for (i=1;i<=l;++i) printf "%s","I"; printf "\n" }' $MATURE_FA > $ANNOTATION_DIR/mature.fq && \
		tailor map \
			-i $ANNOTATION_DIR/mature.fq \
			-p $HAIRPIN_INDEX \
			-n $CPU | \
		samtools view -bS - | \
		bedtools_tailor bamtobed -i - > \
		$ANNOTATION_DIR/mature.annotate_coordinates.bed && \
		rm -rf $ANNOTATION_DIR/mature.fq && \
		touch .${JOBUID}.status.${STEP}.find_annotated_coordinates
	STEP=$((STEP+1))
	
	echo2 "Mapping the input fastq to the hairpin reference" 
	[ ! -f .${JOBUID}.status.${STEP}.tailor_mapping ] && \
		tailor map \
			-i ${PREFIX}.p${MIN_PHRED}.fq \
			-p $HAIRPIN_INDEX \
			-n $CPU \
			$ALLOW_MISMATCH \
			2> $MAPPING_DIR/${PREFIX}.hairpin.tailor.log | \
		tailor_sam_to_bed \
		> $MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.hairpin.bed2 && \
	touch .${JOBUID}.status.${STEP}.tailor_mapping
	STEP=$((STEP+1))
	
	echo2 "Draw overall length distribution with tailing information"
	[ ! -f .${JOBUID}.status.${STEP}.draw_overall_lendis ] && \
		python $PIPELINE_DIRECTORY/bin/tailor_bed2_counter.py \
			$MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.hairpin.bed2 \
			1> $MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.hairpin.sum \
			2> $MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.hairpin.single_nt_sum && \
		Rscript --slave $PIPELINE_DIRECTORY/bin/draw_tailor_lendis.R \
			$MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.hairpin.single_nt_sum \
			$PDF_DIR/${PREFIX}.p${MIN_PHRED}.hairpin.pdf \
			${PREFIX} \
			"total miRNA" && \
	touch .${JOBUID}.status.${STEP}.draw_overall_lendis
	STEP=$((STEP+1))

	echo2 "Adjust the coordinate according to the annotated mature miRNA"
	[ ! -f .${JOBUID}.status.${STEP}.reannotate ] && \
		bedtools_tailor intersect -wo -s -a $MAPPING_DIR/${PREFIX}.p${MIN_PHRED}.hairpin.bed2 -b $ANNOTATION_DIR/mature.annotate_coordinates.bed -f 0.75 | \
		tee $MAPPING_DIR/${PREFIX}.relative.bed.wo | \
		awk 'BEGIN{FS="\t";OFS="\t"}{print $13,$11-$2,$3-$12,$4,$5,$6,$7,$8,$9}' \
		> $MAPPING_DIR/${PREFIX}.relative.bed && \
		touch .${JOBUID}.status.${STEP}.reannotate
	STEP=$((STEP+1))

	echo2 "Draw balloon plot for tailing"
	[ ! -f .${JOBUID}.status.${STEP}.draw_balloon ] && \
		Rscript --slave $PIPELINE_DIRECTORY/bin/draw_tailor_balloon.R  \
			$MAPPING_DIR/${PREFIX}.relative.bed \
			$CPU \
			$PREFIX \
			$BALLOON_DIR && \
		gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PDF_DIR/${PREFIX}.p${MIN_PHRED}.balloon.pdf $BALLOON_DIR/*miRNATailingBalloonPlot.pdf && \
		touch .${JOBUID}.status.${STEP}.draw_balloon
	STEP=$((STEP+1))

fi
