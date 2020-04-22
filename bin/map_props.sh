#/!/bin/sh
##/hps/nobackup2/flicek/user/cander21/LCE20190618/bin/map_props.sh
##this collection of scripts will count up the number of MA and BA sites in the original set of strelka variants. It will also output 100 iterations of the same number of sites from across DEN-sensitive control sites, outputing the number of MA and BA sites. The average proportion of control sites with ALT alleles, considering both BA and MA sites, versus those without any, is reported for each iteration.
##example use: map_props.sh /hps/nobackup2/flicek/user/cander21/LCE20190618/MA/mpileup/${STRAIN}/ 90791_N2 /hps/nobackup2/flicek/user/cander21/LCE20190618/REF/${STRAIN}.20180918.familialFilter.bed /hps/nobackup2/flicek/user/cander21/LCE20190618/MA/controls/${STRAIN}_dysplastic_den.arc_filtered_control_sites.bed.gz
NOD=$1
STRAIN=$2
ALLSITES=$3
CONTROL=$4
MAP=$5

BAS=$(cat ma/${STRAIN}/${NOD}/${NOD}.allele_count.bed | bedtools intersect -a stdin -b ma/${STRAIN}/${NOD}/${NOD}.bed | grep "BA" | wc -l)
MAS=$(cat ma/${STRAIN}/${NOD}/${NOD}.allele_count.bed | bedtools intersect -a stdin -b ma/${STRAIN}/${NOD}/${NOD}.bed | grep "MA" | wc -l)
SUM=$((BAS+MAS))
PROPM=$(echo $MAS | X=$SUM awk -F"\t" '{print ($1/ENVIRON["X"])}')

zcat $ALLSITES | bedtools subtract -a stdin -b ma/${STRAIN}/${NOD}/${NOD}.bed > ma/${STRAIN}/${NOD}/tmp.bed

for x in {1..100} 
 do
 shuf -n $SUM ma/${STRAIN}/${NOD}/tmp.bed | bedtools intersect -a stdin -b ma/${STRAIN}/${NOD}/${NOD}.allele_count.bed -wb > ma/${STRAIN}/${NOD}/tmp2.bed
 BAC=$(grep "BA" ma/${STRAIN}/${NOD}/tmp2.bed | wc -l)
 MAC=$(grep "MA" ma/${STRAIN}/${NOD}/tmp2.bed | wc -l)
 SUC=$((BAC+MAC))
 PROPC=$(echo $SUC | X=$SUM awk -F"\t" '{print ($1/ENVIRON["X"])}')
 echo $BAC $MAC $PROPC >> $CONTROL
 unset PROPC
 unset MAC
 unset BAC 
 unset SUC 
done 

PROPC2=$(cat  $CONTROL | awk '{sum+=$3} END {print sum/NR}')
rm ma/${STRAIN}/${NOD}/tmp2.bed && rm ma/${STRAIN}/${NOD}/tmp.bed && rm $CONTROL
echo $NOD $MAS $BAS $PROPM $PROPC2 > $MAP
###endofscript###
