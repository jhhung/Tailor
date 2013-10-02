#! /bin/bash -x

###############
#configuration#
###############
case ${4} in
human)
# @! THIS IS SOMETHING YOU NEED TO CHANGE
	FOLDER=$PIPELINE_DIRECTORY/common_files/mm9/UCSC_BEDS
	declare -a TARGETS=( \
	)
;;
## mouse specific intersectBed target files
mouse)
# @! THIS IS SOMETHING YOU NEED TO CHANGE
	FOLDER=$PIPELINE_DIRECTORY/common_files/mm9/UCSC_BEDS
	declare -a TARGETS=( \
	)
;;
## fly specific intersectBed target files
fly)
# @! THIS IS SOMETHING YOU NEED TO CHANGE
	FOLDER=$PIPELINE_DIRECTORY/common_files/dm3/UCSC_BEDS
# get rid of rtRNA first
	rtRNA=$FOLDER/repeat_mask.bed.rtRNA
# @! EXAMPLES OF HOW TO ADD NEW FEATURE (NEED TO BE COMPATIBLE WITH BEDTOOLS)
	
	FLY_flyBase_GENE=$FOLDER/UCSC.flyBase.Genes.bed
	FLY_REPEATMASKER=$FOLDER/repeat_mask.bed
	declare -a TARGETS=( \
	"FLY_flyBase_GENE" \
	"FLY_REPEATMASKER" )
;;
*)
	echo "unknown orgnanism... currently only mouse/fly is supported"
	exit 2
;;
esac

#####################
#variable assignment#
#####################
UNIQ_BED=${1%bed2}x_rpmk_rtRNA.bed2
PURE_MULTI_BED=${2%bed2}x_rpmk_rtRNA.bed2
OUT=$3
FA=$5
[ ! -z $6 ] && CPU=$6 || CPU=8
LO_RANGE=16
HI_RANGE=50
LENGTH_RANGE=`seq ${LO_RANGE} ${HI_RANGE}`

# get rid of rRNA mappers
echo "bedtools intersect -v -wa -a $1 -b $rtRNA > ${UNIQ_BED}" > x_rpmk_rtRNA.para
echo "bedtools intersect -v -wa -a $2 -b $rtRNA > ${PURE_MULTI_BED}" >> x_rpmk_rtRNA.para
ParaFly -c x_rpmk_rtRNA.para -CPU $CPU -failed_cmds x_rpmk_rtRNA.para.failed_commands

##############
#print header#
##############
echo -ne "Sample\tTotal_Perfect_Unique_Reads\tTotal_Perfect_Unique_Species\tTotal_Perfect_Multiple_Reads\tTotal_Perfect_Multiple_Species\tTotal_Prefix_Unique_Reads\tTotal_Prefix_Unique_Species\tTotal_Prefix_Multiple_Reads\tTotal_Prefix_Multiple_Species\t" > $OUT;
for t in ${TARGETS[@]}
do \
	echo -ne "${t}_perfect_uniq_reads\t${t}_perfect_uniq_species\t${t}_prefix_uniq_reads\t${t}_prefix_uniq_species\t${t}_sense_perfect_uniq_reads\t${t}_sense_perfect_uniq_species\t${t}_sense_prefix_uniq_reads\t${t}_sense_prefix_uniq_species\t${t}_antisense_perfect_uniq_reads\t${t}_antisense_perfect_uniq_species\t${t}_antisense_prefix_uniq_reads\t${t}_antisense_prefix_uniq_species\t" >> $OUT;
	echo -ne "${t}_perfect_multi_reads\t${t}_perfect_multi_species\t${t}_prefix_multi_reads\t${t}_prefix_multi_species\t${t}_sense_perfect_multi_reads\t${t}_sense_perfect_multi_species\t${t}_sense_prefix_multi_reads\t${t}_sense_prefix_multi_species\t${t}_antisense_perfect_multi_reads\t${t}_antisense_perfect_multi_species\t${t}_antisense_prefix_multi_reads\t${t}_antisense_prefix_multi_species\t" >> $OUT;
done
echo -ne "\n" >> $OUT;

##################
#calculate counts#
##################
echo -ne "${UNIQ_BED%%.bed*}\t" >> $OUT;
TOTAL_PERFECT_UNIQ_READS=`awk '{if ($8==255) a+=$4}END{printf "%d", a}' $UNIQ_BED`
TOTAL_PERFECT_UNIQ_SPECIES=`awk '{if ($8==255) a[$7]=1}END{printf "%d", length(a)}' $UNIQ_BED`
TOTAL_PREFIX_UNIQ_READS=`awk '{if ($8!=255) a+=$4}END{printf "%d", a}' $UNIQ_BED`
TOTAL_PREFIX_UNIQ_SPECIES=`awk '{if ($8!=255) a[$7]=1}END{printf "%d", length(a)}' $UNIQ_BED`
echo -ne "$TOTAL_PERFECT_UNIQ_READS\t$TOTAL_PERFECT_UNIQ_SPECIES\t$TOTAL_PREFIX_UNIQ_READS\t$TOTAL_PREFIX_UNIQ_SPECIES\t" >> $OUT;
awk '{ if ($8==255)perfect_ct[$7]=$4; else prefix_ct[$7]=$4 }END{perfect=0; prefix=0; for (seq in perfect_ct){perfect+=perfect_ct[seq]}; for (seq in prefix_ct){prefix+=prefix_ct[seq]}; printf "%d\t%d\t%d\t%d\t", perfect, length(perfect_ct), prefix, length(prefix_ct)}' $PURE_MULTI_BED >> $OUT

##############
#intersecting#
##############
# -wa -wb uses too much space... so change back to -u -wa
parafly_file="intersect1".para && \
rm -rf $parafly_file
# processing data
for t in ${TARGETS[@]}
do \
	echo "awk '\$8 == 255' $UNIQ_BED | bedtools intersect -f 0.5 -wa -u -a - -b  ${!t} > ${UNIQ_BED%bed2}${t}.perfect.bed2 && bed2lendis ${UNIQ_BED%bed2}${t}.perfect.bed2 > ${UNIQ_BED%bed2}${t}.perfect.bed2.lendis" >> $parafly_file ; 
	echo "awk '\$8 != 255' $UNIQ_BED | bedtools intersect -f 0.5 -wa -u -a - -b  ${!t} > ${UNIQ_BED%bed2}${t}.prefix.bed2  && bed2lendis ${UNIQ_BED%bed2}${t}.prefix.bed2 > ${UNIQ_BED%bed2}${t}.prefix.bed2.lendis" >> $parafly_file ; 
	echo "awk '\$8 == 255' $PURE_MULTI_BED | bedtools intersect -f 0.5 -wa -u -a - -b  ${!t} > ${PURE_MULTI_BED%bed2}${t}.perfect.bed2 && bed2lendis ${PURE_MULTI_BED%bed2}${t}.perfect.bed2 > ${PURE_MULTI_BED%bed2}${t}.perfect.bed2.lendis" >> $parafly_file ; 
	echo "awk '\$8 != 255' $PURE_MULTI_BED | bedtools intersect -f 0.5 -wa -u -a - -b  ${!t} > ${PURE_MULTI_BED%bed2}${t}.prefix.bed2 && bed2lendis ${PURE_MULTI_BED%bed2}${t}.prefix.bed2 > ${PURE_MULTI_BED%bed2}${t}.prefix.bed2.lendis" >> $parafly_file ; 
done
# running ParaFly if no jobs has been ran (no .completed file) or it has ran but has some failed (has .failed_commands)
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi

parafly_file="intersect2".para && \
rm -rf $parafly_file
# processing data
for t in ${TARGETS[@]}
do \
	echo "bedtools intersect -f 0.5 -wa -u -a ${UNIQ_BED%bed2}${t}.perfect.bed2       -b  ${!t} -s > ${UNIQ_BED%bed2}${t}.perfect.S.bed2        && bed2lendis ${UNIQ_BED%bed2}${t}.perfect.S.bed2        > ${UNIQ_BED%bed2}${t}.perfect.S.bed2.lendis"               >> $parafly_file ;
	echo "bedtools intersect -f 0.5 -wa -u -a ${UNIQ_BED%bed2}${t}.perfect.bed2       -b  ${!t} -S > ${UNIQ_BED%bed2}${t}.perfect.AS.bed2       && bed2lendis ${UNIQ_BED%bed2}${t}.perfect.AS.bed2       > ${UNIQ_BED%bed2}${t}.perfect.AS.bed2.lendis"              >> $parafly_file ;
	echo "bedtools intersect -f 0.5 -wa -u -a ${PURE_MULTI_BED%bed2}${t}.perfect.bed2 -b  ${!t} -s > ${PURE_MULTI_BED%bed2}${t}.perfect.S.bed2  && bed2lendis ${PURE_MULTI_BED%bed2}${t}.perfect.S.bed2  > ${PURE_MULTI_BED%bed2}${t}.perfect.S.bed2.lendis"         >> $parafly_file ;
	echo "bedtools intersect -f 0.5 -wa -u -a ${PURE_MULTI_BED%bed2}${t}.perfect.bed2 -b  ${!t} -S > ${PURE_MULTI_BED%bed2}${t}.perfect.AS.bed2 && bed2lendis ${PURE_MULTI_BED%bed2}${t}.perfect.AS.bed2 > ${PURE_MULTI_BED%bed2}${t}.perfect.AS.bed2.lendis"        >> $parafly_file ;
done
# running ParaFly if no jobs has been ran (no .completed file) or it has ran but has some failed (has .failed_commands)
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi

##################
#calculate counts#
##################
for t in ${TARGETS[@]}
do \
	TOTAL_PERFECT_UNIQ_READS=`awk '{if ($8==255) a+=$4}END{printf "%d", a}' ${UNIQ_BED%bed2}${t}.perfect.bed2`
	TOTAL_PERFECT_UNIQ_SPECIES=`awk '{if ($8==255) a[$7]=1}END{printf "%d", length(a)}' ${UNIQ_BED%bed2}${t}.perfect.bed2`
	TOTAL_PREFIX_UNIQ_READS=`awk '{if ($8!=255) a+=$4}END{printf "%d", a}' ${UNIQ_BED%bed2}${t}.prefix.bed2`
	TOTAL_PREFIX_UNIQ_SPECIES=`awk '{if ($8!=255) a[$7]=1}END{printf "%d", length(a)}' ${UNIQ_BED%bed2}${t}.prefix.bed2`
	echo -ne $TOTAL_PERFECT_UNIQ_READS"\t"$TOTAL_PERFECT_UNIQ_SPECIES"\t"$TOTAL_PREFIX_UNIQ_READS"\t"$TOTAL_PREFIX_UNIQ_SPECIES"\t" >> $OUT;
	
	TOTAL_PERFECT_S_UNIQ_READS=`awk '{if ($8==255) a+=$4}END{printf "%d", a}' ${UNIQ_BED%bed2}${t}.perfect.S.bed2 `
	TOTAL_PERFECT_S_UNIQ_SPECIES=`awk '{if ($8==255) a[$7]=1}END{printf "%d", length(a)}' ${UNIQ_BED%bed2}${t}.perfect.S.bed2`
	TOTAL_PREFIX_S_UNIQ_READS=`awk '{if ($8!=255) a+=$4}END{printf "%d", a}' ${UNIQ_BED%bed2}${t}.prefix.S.bed2`
	TOTAL_PREFIX_S_UNIQ_SPECIES=`awk '{if ($8!=255) a[$7]=1}END{printf "%d", length(a)}' ${UNIQ_BED%bed2}${t}.prefix.S.bed2`
	echo -ne $TOTAL_PERFECT_S_UNIQ_READS"\t"$TOTAL_PERFECT_S_UNIQ_SPECIES"\t"$TOTAL_PREFIX_S_UNIQ_READS"\t"$TOTAL_PREFIX_S_UNIQ_SPECIES"\t" >> $OUT;
	
	TOTAL_PERFECT_AS_UNIQ_READS=`awk '{if ($8==255) a+=$4}END{printf "%d", a}' ${UNIQ_BED%bed2}${t}.perfect.AS.bed2 `
	TOTAL_PERFECT_AS_UNIQ_SPECIES=`awk '{if ($8==255) a[$7]=1}END{printf "%d", length(a)}' ${UNIQ_BED%bed2}${t}.perfect.AS.bed2`
	TOTAL_PREFIX_AS_UNIQ_READS=`awk '{if ($8!=255) a+=$4}END{printf "%d", a}' ${UNIQ_BED%bed2}${t}.prefix.AS.bed2`
	TOTAL_PREFIX_AS_UNIQ_SPECIES=`awk '{if ($8!=255) a[$7]=1}END{printf "%d", length(a)}' ${UNIQ_BED%bed2}${t}.prefix.AS.bed2`
	echo -ne $TOTAL_PERFECT_AS_UNIQ_READS"\t"$TOTAL_PERFECT_AS_UNIQ_SPECIES"\t"$TOTAL_PREFIX_AS_UNIQ_READS"\t"$TOTAL_PREFIX_AS_UNIQ_SPECIES"\t" >> $OUT;
	
	awk '{ if ($8==255)perfect_ct[$7]=$4; else prefix_ct[$7]=$4 }END{perfect=0; prefix=0; for (seq in perfect_ct){perfect+=perfect_ct[seq]}; for (seq in prefix_ct){prefix+=prefix_ct[seq]}; printf "%d\t%d\t%d\t%d\t", perfect, length(perfect_ct), prefix, length(prefix_ct)}' ${PURE_MULTI_BED%bed2}${t}.perfect.bed2    >> $OUT;
	awk '{ if ($8==255)perfect_ct[$7]=$4; else prefix_ct[$7]=$4 }END{perfect=0; prefix=0; for (seq in perfect_ct){perfect+=perfect_ct[seq]}; for (seq in prefix_ct){prefix+=prefix_ct[seq]}; printf "%d\t%d\t%d\t%d\t", perfect, length(perfect_ct), prefix, length(prefix_ct)}' ${PURE_MULTI_BED%bed2}${t}.perfect.S.bed2  >> $OUT;
	awk '{ if ($8==255)perfect_ct[$7]=$4; else prefix_ct[$7]=$4 }END{perfect=0; prefix=0; for (seq in perfect_ct){perfect+=perfect_ct[seq]}; for (seq in prefix_ct){prefix+=prefix_ct[seq]}; printf "%d\t%d\t%d\t%d\t", perfect, length(perfect_ct), prefix, length(prefix_ct)}' ${PURE_MULTI_BED%bed2}${t}.perfect.AS.bed2 >> $OUT;
	rm -rf ${UNIQ_BED%bed2}${t}.bed2  ${PURE_MULTI_BED%bed2}${t}.bed2 ;
done

################################
#draw information content graph#
################################
parafly_file="information_content".para
for t in ${TARGETS[@]}
do \
	touch ${UNIQ_BED%bed2}${t}.perfect.S.bed2 && touch ${UNIQ_BED%bed2}${t}.perfect.AS.bed2 && echo "plot_length_S_AS.sh ${UNIQ_BED%bed2}${t}.perfect.S.bed2 ${UNIQ_BED%bed2}${t}.perfect.AS.bed2" >> $parafly_file	
	touch ${UNIQ_BED%bed2}${t}.prefix.S.bed2 && touch ${UNIQ_BED%bed2}${t}.prefix.AS.bed2 && echo "plot_length_S_AS.sh ${UNIQ_BED%bed2}${t}.prefix.S.bed2 ${UNIQ_BED%bed2}${t}.prefix.AS.bed2" >> $parafly_file	
done
# running ParaFly if no jobs has been ran (no .completed file) or it has ran but has some failed (has .failed_commands)
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi

# transpose the table
ruby -lane 'BEGIN{$arr=[]}; $arr.concat([$F]); END{$arr.transpose.each{ |a| puts a.join ("\t") } }' -F"\t" $OUT > ${OUT}.t  && \
mv ${OUT}.t ${OUT} 

###########################
#print length distribution#
###########################
# print header
rm -rf ${OUT}.lendis
for i in `echo $LENGTH_RANGE`;do echo -ne $i"\t" >> ${OUT}.lendis; done
# append newline to the end of lendis
echo >> ${OUT}.lendis
echo >> ${OUT}.lendis;
echo >> ${OUT}.lendis
echo >> ${OUT}.lendis;
echo >> ${OUT}.lendis;

for t in ${TARGETS[@]}
do \
	for lendis in \
		${UNIQ_BED%bed2}${t}.perfect.bed2.lendis \
		${UNIQ_BED%bed2}${t}.perfect.S.bed2.lendis \
		${UNIQ_BED%bed2}${t}.perfect.AS.bed2.lendis \
		${UNIQ_BED%bed2}${t}.prefix.bed2.lendis \
		${UNIQ_BED%bed2}${t}.prefix.S.bed2.lendis \
		${UNIQ_BED%bed2}${t}.prefix.AS.bed2.lendis \
		${PURE_MULTI_BED%bed2}${t}.perfect.bed2.lendis \
		${PURE_MULTI_BED%bed2}${t}.perfect.S.bed2.lendis \
		${PURE_MULTI_BED%bed2}${t}.perfect.AS.bed2.lendis \
		${PURE_MULTI_BED%bed2}${t}.prefix.bed2.lendis \
		${PURE_MULTI_BED%bed2}${t}.prefix.S.bed2.lendis \
		${PURE_MULTI_BED%bed2}${t}.prefix.AS.bed2.lendis;
	do \
		echo "awk '{l[\$1]=\$2}END{for (i=$LO_RANGE;i<=$HI_RANGE;++i){printf \"%d\\n\",l[i]?l[i]:0}}' ${lendis}" | sh | ruby -lane 'BEGIN{$arr=[]}; $arr.concat([$F]); END{$arr.transpose.each{ |a| puts a.join ("\t") } }' -F"\t" >> ${OUT}.lendis && \
		echo >> ${OUT}.lendis;
	done
done

#####################
#compress everything#
#####################
paste ${OUT} ${OUT}.lendis  | awk 'BEGIN{getline;print;head=$0;nf=NF}{if (NF==nf) print}' > ${OUT}.with_lendis 
parafly_file="compress".para
for t in ${TARGETS[@]}
do \
	echo "gzip ${UNIQ_BED%bed2}${t}.perfect.bed2" >> $parafly_file
	echo "gzip ${UNIQ_BED%bed2}${t}.perfect.S.bed2" >> $parafly_file
	echo "gzip ${UNIQ_BED%bed2}${t}.perfect.AS.bed2" >> $parafly_file
	echo "gzip ${UNIQ_BED%bed2}${t}.prefix.bed2" >> $parafly_file
	echo "gzip ${UNIQ_BED%bed2}${t}.prefix.S.bed2" >> $parafly_file
	echo "gzip ${UNIQ_BED%bed2}${t}.prefix.AS.bed2" >> $parafly_file
	echo "gzip ${PURE_MULTI_BED%bed2}${t}.perfect.bed2" >> $parafly_file
	echo "gzip ${PURE_MULTI_BED%bed2}${t}.perfect.S.bed2" >> $parafly_file
	echo "gzip ${PURE_MULTI_BED%bed2}${t}.perfect.AS.bed2" >> $parafly_file
	echo "gzip ${PURE_MULTI_BED%bed2}${t}.prefix.bed2" >> $parafly_file
	echo "gzip ${PURE_MULTI_BED%bed2}${t}.prefix.S.bed2" >> $parafly_file
	echo "gzip ${PURE_MULTI_BED%bed2}${t}.prefix.AS.bed2" >> $parafly_file
done
# running ParaFly if no jobs has been ran (no .completed file) or it has ran but has some failed (has .failed_commands)
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi

###########
#join pdfs#
###########
PDF_NAMES=""
for t in ${TARGETS[@]}
do \
	[ -f ${UNIQ_BED%bed2}${t}.perfect.lendis.pdf ] && PDF_NAMES=${PDF_NAMES}" "${UNIQ_BED%bed2}${t}.perfect.lendis.pdf
	[ -f ${UNIQ_BED%bed2}${t}.prefix.lendis.pdf ] && PDF_NAMES=${PDF_NAMES}" "${UNIQ_BED%bed2}${t}.prefix.lendis.pdf
	[ -f ${UNIQ_BED%bed2}${t}.perfect.ping-pong.pdf ] && PDF_NAMES=${PDF_NAMES}" "${UNIQ_BED%bed2}${t}.perfect.ping-pong.pdf
	[ -f ${UNIQ_BED%bed2}${t}.perfect.S.bed2.reads.5end_60.percentage.pdf ]  && PDF_NAMES=${PDF_NAMES}" "${UNIQ_BED%bed2}${t}.perfect.S.bed2.reads.5end_60.percentage.pdf
	[ -f ${UNIQ_BED%bed2}${t}.prefix.S.bed2.reads.5end_60.percentage.pdf ]  && PDF_NAMES=${PDF_NAMES}" "${UNIQ_BED%bed2}${t}.prefix.S.bed2.reads.5end_60.percentage.pdf
	[ -f ${UNIQ_BED%bed2}${t}.perfect.AS.bed2.reads.5end_60.percentage.pdf ] && PDF_NAMES=${PDF_NAMES}" "${UNIQ_BED%bed2}${t}.perfect.AS.bed2.reads.5end_60.percentage.pdf
	[ -f ${UNIQ_BED%bed2}${t}.prefix.AS.bed2.reads.5end_60.percentage.pdf ] && PDF_NAMES=${PDF_NAMES}" "${UNIQ_BED%bed2}${t}.prefix.AS.bed2.reads.5end_60.percentage.pdf
done
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${1}.features.pdf ${PDF_NAMES} && rm -rf $PDF_NAMES


