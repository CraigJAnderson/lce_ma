##the following shell script sets up the coordinates for plotting multiallelism across the genome, before doing so in R
#/!/bin/sh
##MA_SCE_plot.sh
CHR_BED=$1
STRAIN=$2
NOD=$3
PATH_TO_FNOD=$4 ## /hps/nobackup2/flicek/user/mst/lce/nodules/93133_N1.segs.T.A.fnod
PATH_TO_DRCR=$5 ##/hps/nobackup2/flicek/user/mst/lce/nodules/93133_N1.T.A.fnodDrcr

##paste together SCE segments with new names and /hps/nobackup2/flicek/user/mst/lce/nodules/*.T.A.fnodDrcr #fnod file with the Drcr asymmetry metric (scales 1 to -1) as the final column.
paste <(python bin/make_fnod_beds.py ${PATH_TO_FNOD}/${NOD}.segs.T.A.fnod | sed '$d') <(awk -F"," '{print $3,$6,$15}' ${PATH_TO_DRCR}/$NOD}.T.A.fnodDrcr) | awk '{print $1"\t"$2"\t"$3"\t"$4":"$7}' > ma/${STRAIN}/${NOD}/${NOD}.fnod_bed

while read line ; do BAS=$(echo ${line} | sed 's/ /\t/g' | bedtools intersect -a stdin -b ma/${STRAIN}/${NOD}/${NOD}.allele_count_bed -wb | bedtools intersect -a stdin -b ma/${STRAIN}/${NOD}/${NOD}.variant_bed |grep BA | wc -l) ; MAS=$(echo ${line} | sed 's/ /\t/g' | bedtools intersect -a stdin -b ma/${STRAIN}/${NOD}/${NOD}.allele_count_bed -wb | bedtools intersect -a stdin -b ma/${STRAIN}/${NOD}/${NOD}.variant_bed | grep MA | wc -l) ;SUM=$((BAS+MAS)) ; PROPM=$(echo $MAS | X=$SUM awk -F"\t" '{print ($1/ENVIRON["X"])}') ; echo $line $PROPM | sed 's/_/\t/g' | sed 's/:/\t/g' | awk '{if (NF == 6) {print $1"\t"$2"\t"$3"\t"$4"\t1\t"$5"\t"$6} else if (NF == 7) {print $0}}' | sed 's/\t/ /g' >> ma/${STRAIN}/${NOD}/${NOD}.sce_bed ; done < ma/${STRAIN}/${NOD}/${NOD}.fnod_bed

VAR=$(awk '{$4 = $3 ; total += $4; $4 = total}1' ${CHR_BED} | awk '{ print prev"\t"$0} { prev = $0 }' | sed 1d | sed '$ d' | awk '{print "if ($1 == "$5") {print $1,$2+"$4",$3+"$4",$4,$5,$6,$7} else"}'| xargs)
VAR2=${VAR::-4}
GENOME_LENGTH=$(awk '{$4 = $3 ; total += $4; $4 = total}1' ${CHR_BED} | tail -n 1 |sed 's/ //g')
awk -F" " '{if ($1 == 1) {print $1,$2,$3,$4,$5,$6,$7} else $VAR2 }' ma/${STRAIN}/${NOD}/${NOD}.sce_bed | X=$GENOME_LENGTH awk '{$8 = ($2/ENVIRON["X"]) ; $9 = ($3/ENVIRON["X"]) ; print $0}' > ma/${STRAIN}/${NOD}/${NOD}.sce_bed2
##midpoints and end point of reference if it doesn't already exist
if [ ! -f bin/${STRAIN}_chr_prop ]
then
 X=$GENOME_LENGTH awk '{$4 = $3 ;$5 = $3; total += $4; $4 = ((total-($3/2))/ENVIRON["X"])/2 ; total += $5; $5 = (total/ENVIRON["X"])/2}1' $CHR_BED > bin/${STRAIN}_chr_prop
fi

unset VAR && unset VAR2 && unset GENOME_LENGTH

Rscript --vanilla MA_SCE_plot.R ${STRAIN} ${NOD} 
