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

export MATURE_URL='ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz'
export HAIRPIN_URL='ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz'

#########
# USAGE #
#########
usage() {
cat << EOF

This script obtain mature and hairpin sequence from miRbase and extract sequences for a specific organism.
The output files will be stored in $TAILOR_INDEX

usage: $0 -g dme
 
OPTIONS:
<required>
	-g      three letter string to indicate the organism; for example: 
	          dme for fly
	          mmu for mouse
	          hsa for human

EOF
}

while getopts "hg:" OPTION
do
	case $OPTION in
		h)	usage && exit 0 ;;
		g)	export ORGANISM_STR=$OPTARG ;;
		?)	usage && exit 1 ;;
	esac
done

[[ ! -n $ORGANISM_STR ]] && echo "[Error]: please provide -g option" && exit 1
############
# Download #
############
which wget 1>/dev/null && \
	wget -q -O - $MATURE_URL  | gunzip | faUtoT.py > $TAILOR_INDEX/mature.fa && \
	wget -q -O - $HAIRPIN_URL | gunzip | faUtoT.py > $TAILOR_INDEX/hairpin.fa
[[ $? != 0 ]] && \
which curl && \
	curl $MATURE_URL  | gunzip | faUtoT.py > $TAILOR_INDEX/mature.fa && \
	curl $HAIRPIN_URL | gunzip | faUtoT.py > $TAILOR_INDEX/hairpin.fa
[ ! -s $TAILOR_INDEX/mature.fa -o ! -s $TAILOR_INDEX/hairpin.fa ] && echo "[Error]: fail to obtain fasta file from miRbase" && exit 1

####################
# extract organism #
####################
extract_org_from_fa.py $TAILOR_INDEX/mature.fa  $ORGANISM_STR > $TAILOR_INDEX/${ORGANISM_STR}.mature.fa
extract_org_from_fa.py $TAILOR_INDEX/hairpin.fa $ORGANISM_STR > $TAILOR_INDEX/${ORGANISM_STR}.hairpin.fa




