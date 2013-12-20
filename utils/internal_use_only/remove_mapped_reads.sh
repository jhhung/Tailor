
bowtie="/home/andy/andy/bowtie-1.0.0/bowtie"



remove_unmapped()
{
	name=$1
	tail=$2
	n=$3
	
	INPUT_FQ=$name"."$tail".randomeTailed.fq"
	UN_MAPPED_FQ=$name"."$tail".randomeTailed.unmapped.fq"
	BOWTIE_INDEX="/home/andy/andy/Tailor/utils/internal_use_only/test_index/Drosophila_melanogaster/BowtieIndex/genome"
	
	$bowtie -S $SAM_HEADER -p 4 -v 0 -a -3 $n --un $UN_MAPPED_FQ.tmp $BOWTIE_INDEX $INPUT_FQ > /dev/null
	
	head -n 2000000 $UN_MAPPED_FQ.tmp > $UN_MAPPED_FQ
	
	#rm $UN_MAPPED_FQ.tmp
	
	wc -l $UN_MAPPED_FQ
}

name="Drosophila_melanogaster"

remove_unmapped $name "1" "0"
remove_unmapped $name "2" "1"
remove_unmapped $name "3" "2"
remove_unmapped $name "4" "3"

cat $name.*.randomeTailed.unmapped.fq > $name.all.randomeTailed.fq