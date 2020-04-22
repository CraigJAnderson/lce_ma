#snakefile: /hps/nobackup2/flicek/user/cander21/dojo/SM/snakefile
#import numpy as np

configfile : "config.yaml"

#include : "rules/test.smk"

import pandas as pd
sample_names = pd.read_table(config["SAMPLE_LIST"],sep=" ",header=None)
sample_names.columns = ['bam','nod','NA']
samples= list(sample_names.nod)
bam_names = []
for x in list(sample_names.bam):
 tmp = config["PATH_TO_BAM"]+x+".bam"
 bam_names.append(tmp)

sample_linked_bams = dict(zip(samples, bam_names))

rule all:
 input: 
  SAMPLE_NAMES = expand("geno/{strain}/{nod}_arc_filtered.nodMat", nod=samples, strain=config["STRAIN"]),
  CONTROL_POS = expand("ma/{strain}/{strain}_den.arc_filtered_control_sites.pos", strain=config["STRAIN"]),
  CONTROL_BED_GZ = expand("ma/{strain}/{strain}_den.arc_filtered_control_sites.bed.gz", strain=config["STRAIN"]),
  MPILEUP_OUT = expand("ma/{strain}/{nod}/{nod}.mpileup", nod=samples, strain=config["STRAIN"]),
  ALLELE_COUNT = expand("ma/{strain}/{nod}/{nod}.allele_count.bed", nod=samples, strain=config["STRAIN"]),
  PROPORTION_MULTIALLELIC_CONTROL_SITES = expand("ma/{strain}/{nod}/{nod}.control_map", nod=samples, strain=config["STRAIN"]),
  PROPORTION_MULTIALLELIC = expand("ma/{strain}/{nod}/{nod}.map", nod=samples, strain=config["STRAIN"]),
  ALL_PROPORTIONS = expand("ma/{strain}/{strain}.all_map", strain=config["STRAIN"]),
  NOD_VARIANT_BED = expand("ma/{strain}/{nod}/{nod}.bed", nod=samples, strain=config["STRAIN"]),
  PROXIMAL_VARIANTS = expand("ma/{strain}/{nod}/{nod}.prox.pos", nod=samples, strain=config["STRAIN"]),
  PROXIMAL_VARIANT_MPILEUP = expand("ma/{strain}/{nod}/{nod}.prox.mpileup", nod=samples, strain=config["STRAIN"]),
  ALLELIC_COMBINATINS = expand("ma/{strain}/{nod}/{nod}.MA_combos", nod=samples, strain=config["STRAIN"]),
  COMBINATION_CLASS_COUNT = expand("ma/{strain}/{nod}/{nod}.combos_txt", nod=samples, strain=config["STRAIN"]),
#  ALL_COMBOS = expand("ma/{strain}/{strain}.all_combos.txt", nod=samples, strain=config["STRAIN"]),

ruleorder : filter_variants > collect_control_coordinates > genotype > allele_count > proportion_multiallelic

rule filter_variants:
 output: 
  "geno/{strain}/{nod}_arc_filtered.nodMat"
 params:
  config["VARIANTS"]
 shell:
  """grep -e ',{wildcards.nod},' {params} | awk -F, '{{if (($19 == "") && ($20 == 0) && ($34 == "" || $34 > 1e-4 ) && ($10 != "")) print $0}}' > {output} """

rule collect_control_coordinates:
 input:
  VARIANTS =expand("geno/"+"{strain}"+"/{nod}_arc_filtered.nodMat", nod=samples, strain=config["STRAIN"])
 output:
  CONTROL_POS = "ma/{strain}/{strain}_den.arc_filtered_control_sites.pos",
  CONTROL_BED = "ma/{strain}/{strain}_den.arc_filtered_control_sites.bed",
  CONTROL_BED_GZ = "ma/{strain}/{strain}_den.arc_filtered_control_sites.bed.gz"
 shell:
  """ cat {input.VARIANTS} | awk '{{print $1"\t"$2}}' FS="[_,]" | sort -k 1,1 -k2,2n -u > {output.CONTROL_POS} ; """
  """ cat {input.VARIANTS} | awk '{{print $1"\t"$2"\t"$2+1}}' FS="[_,]" | sort -k 1,1 -k2,2n -u > {output.CONTROL_BED} ; """
  """ bgzip -c {output.CONTROL_BED} > {output.CONTROL_BED_GZ} """

rule genotype:
 input:
  CONTROL_POS = "ma/{strain}/{strain}_den.arc_filtered_control_sites.pos",
 params:
  GENOME = config["GENOME"],
  BAM = lambda wildcards : sample_linked_bams[wildcards.nod]
 output:
  MPILEUP = "ma/{strain}/{nod}/{nod}.mpileup",
 shell:
  """samtools mpileup {params.BAM} -l {input.CONTROL_POS} -f {params.GENOME} -Q 20 --skip-indels --ff DUP | awk -v FS="\t" -v OFS="\t" '{{gsub(/\./,",",$5)}}1' | sed 's/[+-]100[ACGTNacgtn]{{100}}//g' | sed 's/[+-]99[ACGTNacgtn]{{99}}//g' | sed 's/[+-]98[ACGTNacgtn]{{98}}//g' | sed 's/[+-]97[ACGTNacgtn]{{97}}//g' | sed 's/[+-]96[ACGTNacgtn]{{96}}//g' | sed 's/[+-]95[ACGTNacgtn]{{95}}//g' | sed 's/[+-]94[ACGTNacgtn]{{94}}//g' | sed 's/[+-]93[ACGTNacgtn]{{93}}//g' | sed 's/[+-]92[ACGTNacgtn]{{92}}//g' | sed 's/[+-]91[ACGTNacgtn]{{91}}//g' | sed 's/[+-]90[ACGTNacgtn]{{90}}//g' | sed 's/[+-]89[ACGTNacgtn]{{89}}//g' | sed 's/[+-]88[ACGTNacgtn]{{88}}//g' | sed 's/[+-]87[ACGTNacgtn]{{87}}//g' | sed 's/[+-]86[ACGTNacgtn]{{86}}//g' | sed 's/[+-]85[ACGTNacgtn]{{85}}//g' | sed 's/[+-]84[ACGTNacgtn]{{84}}//g' | sed 's/[+-]83[ACGTNacgtn]{{83}}//g' | sed 's/[+-]82[ACGTNacgtn]{{82}}//g' | sed 's/[+-]81[ACGTNacgtn]{{81}}//g' | sed 's/[+-]80[ACGTNacgtn]{{80}}//g' | sed 's/[+-]79[ACGTNacgtn]{{79}}//g' | sed 's/[+-]78[ACGTNacgtn]{{78}}//g' | sed 's/[+-]77[ACGTNacgtn]{{77}}//g' | sed 's/[+-]76[ACGTNacgtn]{{76}}//g' | sed 's/[+-]75[ACGTNacgtn]{{75}}//g' | sed 's/[+-]74[ACGTNacgtn]{{74}}//g' | sed 's/[+-]73[ACGTNacgtn]{{73}}//g' | sed 's/[+-]72[ACGTNacgtn]{{72}}//g' | sed 's/[+-]71[ACGTNacgtn]{{71}}//g' | sed 's/[+-]70[ACGTNacgtn]{{70}}//g' | sed 's/[+-]69[ACGTNacgtn]{{69}}//g' | sed 's/[+-]68[ACGTNacgtn]{{68}}//g' | sed 's/[+-]67[ACGTNacgtn]{{67}}//g' | sed 's/[+-]66[ACGTNacgtn]{{66}}//g' | sed 's/[+-]65[ACGTNacgtn]{{65}}//g' | sed 's/[+-]64[ACGTNacgtn]{{64}}//g' | sed 's/[+-]63[ACGTNacgtn]{{63}}//g' | sed 's/[+-]62[ACGTNacgtn]{{62}}//g' | sed 's/[+-]61[ACGTNacgtn]{{61}}//g' | sed 's/[+-]60[ACGTNacgtn]{{60}}//g' | sed 's/[+-]59[ACGTNacgtn]{{59}}//g' | sed 's/[+-]58[ACGTNacgtn]{{58}}//g' | sed 's/[+-]57[ACGTNacgtn]{{57}}//g' | sed 's/[+-]56[ACGTNacgtn]{{56}}//g' | sed 's/[+-]55[ACGTNacgtn]{{55}}//g' | sed 's/[+-]54[ACGTNacgtn]{{54}}//g' | sed 's/[+-]53[ACGTNacgtn]{{53}}//g' | sed 's/[+-]52[ACGTNacgtn]{{52}}//g' | sed 's/[+-]51[ACGTNacgtn]{{51}}//g' | sed 's/[+-]50[ACGTNacgtn]{{50}}//g' | sed 's/[+-]49[ACGTNacgtn]{{49}}//g' | sed 's/[+-]48[ACGTNacgtn]{{48}}//g' | sed 's/[+-]47[ACGTNacgtn]{{47}}//g' | sed 's/[+-]46[ACGTNacgtn]{{46}}//g' | sed 's/[+-]45[ACGTNacgtn]{{45}}//g' | sed 's/[+-]44[ACGTNacgtn]{{44}}//g' | sed 's/[+-]43[ACGTNacgtn]{{43}}//g' | sed 's/[+-]42[ACGTNacgtn]{{42}}//g' | sed 's/[+-]41[ACGTNacgtn]{{41}}//g' | sed 's/[+-]40[ACGTNacgtn]{{40}}//g' | sed 's/[+-]39[ACGTNacgtn]{{39}}//g' | sed 's/[+-]38[ACGTNacgtn]{{38}}//g' | sed 's/[+-]37[ACGTNacgtn]{{37}}//g' | sed 's/[+-]36[ACGTNacgtn]{{36}}//g' | sed 's/[+-]35[ACGTNacgtn]{{35}}//g' | sed 's/[+-]34[ACGTNacgtn]{{34}}//g' | sed 's/[+-]33[ACGTNacgtn]{{33}}//g' | sed 's/[+-]32[ACGTNacgtn]{{32}}//g' | sed 's/[+-]31[ACGTNacgtn]{{31}}//g' | sed 's/[+-]30[ACGTNacgtn]{{30}}//g' | sed 's/[+-]29[ACGTNacgtn]{{29}}//g' | sed 's/[+-]28[ACGTNacgtn]{{28}}//g' | sed 's/[+-]27[ACGTNacgtn]{{27}}//g' | sed 's/[+-]26[ACGTNacgtn]{{26}}//g' | sed 's/[+-]25[ACGTNacgtn]{{25}}//g' | sed 's/[+-]24[ACGTNacgtn]{{24}}//g' | sed 's/[+-]23[ACGTNacgtn]{{23}}//g' | sed 's/[+-]22[ACGTNacgtn]{{22}}//g' | sed 's/[+-]21[ACGTNacgtn]{{21}}//g' | sed 's/[+-]20[ACGTNacgtn]{{20}}//g' | sed 's/[+-]19[ACGTNacgtn]{{19}}//g' | sed 's/[+-]18[ACGTNacgtn]{{18}}//g' | sed 's/[+-]17[ACGTNacgtn]{{17}}//g' | sed 's/[+-]16[ACGTNacgtn]{{16}}//g' | sed 's/[+-]15[ACGTNacgtn]{{15}}//g' | sed 's/[+-]14[ACGTNacgtn]{{14}}//g' | sed 's/[+-]13[ACGTNacgtn]{{13}}//g' | sed 's/[+-]12[ACGTNacgtn]{{12}}//g' | sed 's/[+-]11[ACGTNacgtn]{{11}}//g' | sed 's/[+-]10[ACGTNacgtn]{{10}}//g' | sed 's/[+-]9[ACGTNacgtn]{{9}}//g' | sed 's/[+-]8[ACGTNacgtn]{{8}}//g' | sed 's/[+-]7[ACGTNacgtn]{{7}}//g' | sed 's/[+-]6[ACGTNacgtn]{{6}}//g' | sed 's/[+-]5[ACGTNacgtn]{{5}}//g' | sed 's/[+-]4[ACGTNacgtn]{{4}}//g' | sed 's/[+-]3[ACGTNacgtn]{{3}}//g' | sed 's/[+-]2[ACGTNacgtn]{{2}}//g' | sed 's/[+-]1[ACGTNacgtn]{{1}}//g' | awk -v FS="\t" -v OFS="\t" '{{gsub(/$/,"",$5)}}1' | sed 's/\*/,/g' | sed 's/\^.//g' | sed 's/\$//g' | awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' > {output.MPILEUP} """

rule allele_count:
 input:
  expand("ma/{strain}/{nod}/{nod}.mpileup", nod=samples, strain=config["STRAIN"])
 output:
  "ma/{strain}/{nod}/{nod}.allele_count.bed"
 shell:
  """python bin/mp_allele_counter.py {input} > {output}"""

rule proportion_multiallelic:
 input:
  VARIANTS = expand("geno/{strain}/{nod}_arc_filtered.nodMat", nod=samples, strain=config["STRAIN"]),
  CONTROL_SITES = expand("ma/{strain}/{strain}_den.arc_filtered_control_sites.bed.gz", strain=config["STRAIN"])
 output:
  VARIANT_BED = "ma/{strain}/{nod}/{nod}.bed",
  PROPORTION_MULTIALLELIC = "ma/{strain}/{nod}/{nod}.map",
  PROPORTION_MULTIALLELIC_CONTROL_SITES = "ma/{strain}/{nod}/{nod}.control_map"
 shell:
  """ cat {input.VARIANTS} | awk -F",|_" '{{print $1"\t"$2"\t"$2+1}}' | sort -k 1,1 -k2,2n > {output.VARIANT_BED} ; """
  """ bin/map_props.sh {wildcards.nod} {wildcards.strain} {input.CONTROL_SITES} {output.PROPORTION_MULTIALLELIC_CONTROL_SITES} {output.PROPORTION_MULTIALLELIC} """

rule all_proportions:
 input:
  expand("ma/{strain}/{nod}/{nod}.map", nod=samples, strain=config["STRAIN"])
 output:
  "ma/{strain}/{strain}.all_map"
 shell:
  """ cat {input} >> {output} """

rule multiallelic_combinations:
 params:
  GENOME = config["GENOME"],
  BAM = lambda wildcards : sample_linked_bams[wildcards.nod]
 output:
  PROXIMAL_VARIANTS = "ma/{strain}/{nod}/{nod}.prox.pos",
  PROXIMAL_VARIANT_MPILEUP = "ma/{strain}/{nod}/{nod}.prox.mpileup",
  ALLELIC_COMBINATINS = "ma/{strain}/{nod}/{nod}.MA_combos"
 shell:
  """bin/combo_validation.sh {wildcards.strain} {params.BAM} {wildcards.nod} {params.GENOME} """

rule allelic_combination_grader:
 input:
  "ma/{strain}/{nod}/{nod}.MA_combos"
 output:
  IND_COMBOS = "ma/{strain}/{nod}/{nod}.combos_txt",
 shell:
  """ VAR=$(egrep $'\t100\t|\t111\t|\t101\t|\t110\t' {input} | awk '{{if ($5 <=150) {{print $0}}}}' | awk '{{if ($5 >=3) {{print $0}}}}' | awk -F"\t" '{{print $2" "$4" "$8}}' | awk '{{ gsub("[(),0-9\047 ]","",$0); print $0 }}'| sed 's/.\{{2\}}/&,/g'| sed 's/./&,/1'| sed -r 's/,/ /' | sed -r 's/,/ /' | sed 's/,$//g' | python bin/allelic_combination_grader.py | cat - bin/allelic_combination_grades.txt | sed '/^$/d' | sort -k2,2 | uniq -c | egrep -v " 0$" | awk '{{print $1}}' | tr '\n' ' ' | awk '{{print $1-1"\t"$2-1"\t"$3-1"\t"$4-1"\t"$5-1"\t"(($3-1)+($4-1)+($5-1))/(($1-1)+($2-1)+($3-1)+($4-1)+($5-1))}}') ; VAR2=$(grep "{wildcards.nod} " ma/c3h/c3h.all_map | awk '{{print $4}}') ; echo {wildcards.nod} $VAR $VAR2 | awk -F" " '{{ if ($7 != -1) print $0}}' >> {output.IND_COMBOS} ; """
