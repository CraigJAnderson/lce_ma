#/!/bin/sh
#calc_timing_props.sh
STRAIN=$1
NOD=$2
TUMOUR_INFO=$3 ##/hps/nobackup2/flicek/user/mst/lce/nodules/c3h.tumourCollateInfo.tab

MA=$(awk -F" " '{ if ($7 >= 0.04) print $0}' ma/${STRAIN}/${NOD}/${NOD}.sce_bed2 | wc -l)
BA=$(awk -F" " '{ if ($7 < 0.04) print $0}' ma/${STRAIN}/${NOD}/${NOD}.sce_bed2 | wc -l)
SUM=$((BA+MA))
PROPM=$(echo $MA | X=$SUM awk -F"\t" '{print ($1/ENVIRON["X"])}')
if [ $PROPM == 1.0 ] ; then VAR1=1 ; else VAR1=0 ; fi
paste <(echo $NOD $PROPM $VAR1) <(grep -e "^${NOD}," ${TUMOUR_INFO} | awk -F"," '{print $28}' ) | awk -F" " '{if ($4 ~ /Braf/) {print $0"\t1"} else {print $0"\t0"}}' | awk -F" " '{if ($4 ~ /Egfr/) {print $0"\t1"} else {print $0"\t0"}}' | awk -F" " '{if ($4 ~ /Hras/) {print $0"\t1"} else {print $0"\t0"}}' | awk -F" " '{if ($4 ~ /Kras/) {print $0"\t1"} else {print $0"\t0"}}' > ma/${STRAIN}/${NOD}/${NOD}.timing
