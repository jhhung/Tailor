#!/bin/bash
senseBed2=$1
antisenseBed2=$2
mergedBed=${antisenseBed2%.AS.bed2}.lendis
awk '{if (ARGIND==1) {s[$1]=$2} else {as[$1]=$2}}END{m=0.0;for (l in s){l=l+0;if (l>m) m=l} for (l in as){l=l+0; if (l>m) m=l} for (len=1;len<=m;++len) {print len"\t"(s[len]?s[len]:0)"\t"(as[len]?as[len]:0)}}' ${senseBed2}.lendis ${antisenseBed2}.lendis > ${mergedBed} && \
Rscript ${PIPELINE_DIRECTORY}/bin/draw_lendis2.R ${mergedBed} ${mergedBed##*/}