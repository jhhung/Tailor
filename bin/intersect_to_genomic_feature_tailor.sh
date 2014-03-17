
TOTAL_BED=$1
UNIQ_BED=${TOTAL_BED%bed2}uniq.bed2
awk '$5==1' $TOTAL_BED > $UNIQ_BED
OUT=${TOTAL_BED%bed2}summary

[ -z $GENOMIC_FEATURE_FILE ] && echo2 "variable GENOMIC_FEATURE_FILE unspecified" "error"
[ -s $GENOMIC_FEATURE_FILE ] || echo2 "file $GENOMIC_FEATURE_FILE is not-exist or empty" "error"
. $GENOMIC_FEATURE_FILE

echo2 "Printing header"
echo -ne "Sample\tTotal_Perfect_Unique_Reads\tTotal_Perfect_Unique_Species\tTotal_Perfect_All_Reads\tTotal_Perfect_All_Species\tTotal_Prefix_Unique_Reads\tTotal_Prefix_Unique_Species\tTotal_Prefix_All_Reads\tTotal_Prefix_All_Species\t" > $OUT;
for t in ${TARGETS[@]}
do \
	echo -ne "${t}_perfect_uniq_reads\t${t}_perfect_uniq_species\t${t}_prefix_uniq_reads\t${t}_prefix_uniq_species\t${t}_sense_perfect_uniq_reads\t${t}_sense_perfect_uniq_species\t${t}_sense_prefix_uniq_reads\t${t}_sense_prefix_uniq_species\t${t}_antisense_perfect_uniq_reads\t${t}_antisense_perfect_uniq_species\t${t}_antisense_prefix_uniq_reads\t${t}_antisense_prefix_uniq_species\t" >> $OUT;
	echo -ne "${t}_perfect_all_reads\t${t}_perfect_all_species\t${t}_prefix_all_reads\t${t}_prefix_all_species\t${t}_sense_perfect_all_reads\t${t}_sense_perfect_all_species\t${t}_sense_prefix_all_reads\t${t}_sense_prefix_all_species\t${t}_antisense_perfect_all_reads\t${t}_antisense_perfect_all_species\t${t}_antisense_prefix_all_reads\t${t}_antisense_prefix_all_species\t" >> $OUT;
done
echo -ne "\n" >> $OUT;

echo2 "Calculating overall statistics"
echo -ne "`basename ${UNIQ_BED%.bed2}`\t" >> $OUT;
TOTAL_PERFECT_UNIQ_READS=`awk   '{if ($9==0) a+=$4}  END{printf "%d", a}' $UNIQ_BED`
TOTAL_PERFECT_UNIQ_SPECIES=`awk '{if ($9==0) a[$7]=1}END{printf "%d", length(a)}' $UNIQ_BED`
TOTAL_PREFIX_UNIQ_READS=`awk    '{if ($9!=0) a+=$4}  END{printf "%d", a}' $UNIQ_BED`
TOTAL_PREFIX_UNIQ_SPECIES=`awk  '{if ($9!=0) a[$7]=1}END{printf "%d", length(a)}' $UNIQ_BED`
echo -ne "$TOTAL_PERFECT_UNIQ_READS\t$TOTAL_PERFECT_UNIQ_SPECIES\t$TOTAL_PREFIX_UNIQ_READS\t$TOTAL_PREFIX_UNIQ_SPECIES\t" >> $OUT;
awk '{ if ($9==0)perfect_ct[$7]=$4; else prefix_ct[$7]=$4 }END{perfect=0; prefix=0; for (seq in perfect_ct){perfect+=perfect_ct[seq]}; for (seq in prefix_ct){prefix+=prefix_ct[seq]}; printf "%d\t%d\t%d\t%d\t", perfect, length(perfect_ct), prefix, length(prefix_ct)}' $TOTAL_BED >> $OUT

echo2 "Intersecting BED file to different genomic feature"
parafly_file="intersect1".para && \
rm -rf $parafly_file
for t in ${TARGETS[@]}
do \
	echo "bedtools intersect -f 0.99 -wa -u -a ${TOTAL_BED} -b ${!t} > ${TOTAL_BED%bed2}${t}.bed2" >> $parafly_file
done
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi
rm -rf ${parafly_file}*

echo2 "Draw figures for different genomic structure"
parafly_file="draw_fig".para
for t in ${TARGETS[@]}
do \
		echo "python $PIPELINE_DIRECTORY/bin/tailor_bed2_counter.py ${TOTAL_BED%bed2}${t}.bed2 1> ${TOTAL_BED%bed2}${t}.sum 2> ${TOTAL_BED%bed2}${t}.single_nt_sum && Rscript --slave $PIPELINE_DIRECTORY/bin/draw_tailor_lendis.R ${TOTAL_BED%bed2}${t}.single_nt_sum $PDF_DIR/`basename ${TOTAL_BED%bed2}${t}.single_nt_sum`.pdf ${PREFIX} ${t}" >> $parafly_file
done
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi
rm -rf ${parafly_file}*

echo2 "Joining pdf"
PDF_NAMES=""
for t in ${TARGETS[@]}
do \
	[ -f $PDF_DIR/`basename ${TOTAL_BED%bed2}${t}.single_nt_sum`.pdf ] && PDF_NAMES=${PDF_NAMES}" "$PDF_DIR/`basename ${TOTAL_BED%bed2}${t}.single_nt_sum`.pdf
done
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PDF_DIR/`basename ${TOTAL_BED%bed2}features.pdf` ${PDF_NAMES} && rm -rf $PDF_NAMES
