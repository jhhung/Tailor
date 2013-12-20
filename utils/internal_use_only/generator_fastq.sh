ArtificialFastqGenerator="/home/andy/andy/Tailor/utils/internal_use_only/ArtificialFastqGenerator";

name="Danio_rerio_fastq";
name="Drosophila_melanogaster";

./seqtk/seqtk sample -s100 $ArtificialFastqGenerator/$name.1.fastq 12000000 > $ArtificialFastqGenerator/$name.fq

sed -n '1,3000000p' $ArtificialFastqGenerator/$name.fq > $name.1.fq
sed -n '3000001,6000000p' $ArtificialFastqGenerator/$name.fq > $name.2.fq
sed -n '6000001,9000000p' $ArtificialFastqGenerator/$name.fq > $name.3.fq
sed -n '9000001,12000000p' $ArtificialFastqGenerator/$name.fq > $name.4.fq


gawk 'BEGIN{nt[0]="A"; nt[1]="C"; nt[2]="G"; nt[3]="T"; srand(systime() + PROCINFO["pid"]); C=0} {ranLen=1; ranTail=""; for (i=1;i<=ranLen;++i) ranTail=ranTail""nt[int(rand()*100)%4]; print "@"C++"."(ranTail==""?"p":ranTail); getline; printf "%s%s", $0, ranTail; getline; getline; printf "\n+\n%s", $0; for (i=1;i<=ranLen;++i) printf "I"; printf "\n"}' $name.1.fq > $name.1.randomeTailed.fq

gawk 'BEGIN{nt[0]="A"; nt[1]="C"; nt[2]="G"; nt[3]="T"; srand(systime() + PROCINFO["pid"]); C=0} {ranLen=2; ranTail=""; for (i=1;i<=ranLen;++i) ranTail=ranTail""nt[int(rand()*100)%4]; print "@"C++"."(ranTail==""?"p":ranTail); getline; printf "%s%s", $0, ranTail; getline; getline; printf "\n+\n%s", $0; for (i=1;i<=ranLen;++i) printf "I"; printf "\n"}' $name.2.fq > $name.2.randomeTailed.fq

gawk 'BEGIN{nt[0]="A"; nt[1]="C"; nt[2]="G"; nt[3]="T"; srand(systime() + PROCINFO["pid"]); C=0} {ranLen=3; ranTail=""; for (i=1;i<=ranLen;++i) ranTail=ranTail""nt[int(rand()*100)%4]; print "@"C++"."(ranTail==""?"p":ranTail); getline; printf "%s%s", $0, ranTail; getline; getline; printf "\n+\n%s", $0; for (i=1;i<=ranLen;++i) printf "I"; printf "\n"}' $name.3.fq > $name.3.randomeTailed.fq

gawk 'BEGIN{nt[0]="A"; nt[1]="C"; nt[2]="G"; nt[3]="T"; srand(systime() + PROCINFO["pid"]); C=0} {ranLen=4; ranTail=""; for (i=1;i<=ranLen;++i) ranTail=ranTail""nt[int(rand()*100)%4]; print "@"C++"."(ranTail==""?"p":ranTail); getline; printf "%s%s", $0, ranTail; getline; getline; printf "\n+\n%s", $0; for (i=1;i<=ranLen;++i) printf "I"; printf "\n"}' $name.4.fq > $name.4.randomeTailed.fq

