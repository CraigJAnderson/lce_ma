#snakefile: /hps/nobackup2/flicek/user/cander21/dojo/SM/snakefile
#import numpy as np

configfile : "config_c3h.yaml"

#include : "rules/test.smk"

import pandas as pd
sample_names = pd.read_table(config["SAMPLE_LIST"],sep=" ",header=None)
sample_names.columns = ['bam','nod']
samples= list(sample_names.nod)
bam_names = []
for x in list(sample_names.bam):
 tmp = config["PATH_TO_BAM"]+x+".bam"
 bam_names.append(tmp)

sample_linked_bams = dict(zip(samples, bam_names))

rule all:
 input: 
  CONTROL_POS = expand("ma/{strain}/{strain}.arc_filtered_control_sites.pos", strain=config["STRAIN"]),
  CONTROL_BED = expand("ma/{strain}/{strain}.arc_filtered_control_sites.bed", strain=config["STRAIN"]),
  ALL_PROPORTIONS = expand("ma/{strain}/{strain}.multi_allelic_proportions", strain=config["STRAIN"]),
  ALLELIC_COMBINATIONS = expand("ma/{strain}/{nod}/{nod}.MA_combos", nod=samples, strain=config["STRAIN"]),
  PLOT_MA = expand("ma/{strain}/{nod}/{nod}.ma_plot.pdf", nod=samples, strain=config["STRAIN"]),
  PLOT_TIMING = expand("ma/{strain}/{strain}.timing_plot.pdf", strain=config["STRAIN"])

ruleorder : filter_variants > collect_control_coordinates > genotype > allele_count > proportion_multiallelic > all_proportions > multiallelic_combinations > allelic_combination_grader

rule filter_variants:
 output: 
  "geno/{strain}/{nod}_arc_filtered.nodMat"
 params:
  config["VARIANTS"]
 shell:
  """grep -e ',{wildcards.nod},' {params} | awk -F, '{{if (($19 == "") && ($20 == 0) && ($34 == "" || $34 > 1e-4 ) && ($10 != "")) print $0}}' > {output} """

rule collect_control_coordinates:
 input:
  VARIANTS = expand("geno/{strain}/{nod}_arc_filtered.nodMat", nod=samples, strain=config["STRAIN"])
 output:
  CONTROL_POS = "ma/{strain}/{strain}.arc_filtered_control_sites.pos",
  CONTROL_BED = "ma/{strain}/{strain}.arc_filtered_control_sites.bed",
 shell:
  """ cat {input.VARIANTS} | awk '{{print $1"\t"$2}}' FS="[_,]" | sort -k 1,1 -k2,2n -u > {output.CONTROL_POS} ; """
  """ cat {input.VARIANTS} | awk '{{print $1"\t"$2"\t"$2+1}}' FS="[_,]" | sort -k 1,1 -k2,2n -u > {output.CONTROL_BED} ; """

rule genotype:
 input:
  CONTROL_POS = "ma/{strain}/{strain}.arc_filtered_control_sites.pos",
 params:
  GENOME = config["GENOME"],
  BAM = lambda wildcards : sample_linked_bams[wildcards.nod]
 output:
  MPILEUP = "ma/{strain}/{nod}/{nod}.mpileup",
 shell:
  """bin/genotype.sh {params.BAM} {input.CONTROL_POS} {params.GENOME} {output.MPILEUP}"""

rule allele_count:
 input:
  "ma/{strain}/{nod}/{nod}.mpileup"
 output:
  "ma/{strain}/{nod}/{nod}.allele_count_bed"
 shell:
  """python bin/mp_allele_counter.py {input} > {output}"""

rule proportion_multiallelic:
 input:
  ALLELE_COUNT = "ma/{strain}/{nod}/{nod}.allele_count_bed",
  MPILEUP = "ma/{strain}/{nod}/{nod}.mpileup",
  VARIANTS = "geno/{strain}/{nod}_arc_filtered.nodMat",
  CONTROL_SITES = "ma/{strain}/{strain}.arc_filtered_control_sites.bed"
 output:
  VARIANT_BED = "ma/{strain}/{nod}/{nod}.variant_bed",
  PROPORTION_MULTIALLELIC = "ma/{strain}/{nod}/{nod}.map",
 shell:
  """ awk -F",|_" '{{print $1"\t"$2"\t"$2+1}}' {input.VARIANTS} | sort -k 1,1 -k2,2n > {output.VARIANT_BED} ; """
  """ bin/map_props.sh {wildcards.nod} {wildcards.strain} {input.CONTROL_SITES} ma/{wildcards.strain}/{wildcards.nod}/{wildcards.nod}.control_map {output.PROPORTION_MULTIALLELIC} """

rule all_proportions:
 input:
  expand("ma/{strain}/{nod}/{nod}.map", nod=samples, strain=config["STRAIN"])
 output:
  "ma/{strain}/{strain}.multi_allelic_proportions"
 shell:
  """ cat {input} > {output} """

rule multiallelic_combinations:
 input:
  ALLELE_COUNT_BED = "ma/{strain}/{nod}/{nod}.allele_count_bed",
  VARIANT_BED = "ma/{strain}/{nod}/{nod}.variant_bed"
 params:
  GENOME = config["GENOME"],
  BAM = lambda wildcards : sample_linked_bams[wildcards.nod]
 output:
  PROXIMAL_VARIANTS = "ma/{strain}/{nod}/{nod}.prox_pos",
  PROXIMAL_VARIANT_MPILEUP = "ma/{strain}/{nod}/{nod}.prox_mpileup",
  ALLELIC_COMBINATIONS = "ma/{strain}/{nod}/{nod}.MA_combos"
 shell:
  """bin/combo_validation.sh {wildcards.strain} {params.BAM} {wildcards.nod} {params.GENOME} """

rule allelic_combination_grader:
 input:
  MA_COMBINATIONS = "ma/{strain}/{nod}/{nod}.MA_combos",
  MAP = "ma/{strain}/{nod}/{nod}.map"
 output:
  IND_COMBOS = "ma/{strain}/{nod}/{nod}.combos_txt",
 shell:
  """ VAR=$(egrep $'\t100\t|\t111\t|\t101\t|\t110\t' {input.MA_COMBINATIONS} | awk '{{if ($5 <=150) {{print $0}}}}' | awk '{{if ($5 >=3) {{print $0}}}}' | awk -F"\t" '{{print $2" "$4" "$8}}' | awk '{{ gsub("[(),0-9\\047 ]","",$0); print $0 }}'| sed 's/.\{{2\}}/&,/g'| sed 's/./&,/1'| sed -r 's/,/ /' | sed -r 's/,/ /' | sed 's/,$//g' | python bin/allelic_combination_grader.py | cat - bin/allelic_combination_grades.txt | sed '/^$/d' | sort -k2,2 | uniq -c | egrep -v " 0$" | awk '{{print $1}}' | tr '\n' ' ' | awk '{{print $1-1"\t"$2-1"\t"$3-1"\t"$4-1"\t"$5-1"\t"(($3-1)+($4-1)+($5-1))/(($1-1)+($2-1)+($3-1)+($4-1)+($5-1))}}') ; VAR2=$(grep '{wildcards.nod} ' {input.MAP} | awk '{{print $4}}') ; echo {wildcards.nod} $VAR $VAR2 | awk -F" " '{{ if ($7 != -1) print $0}}' > {output.IND_COMBOS} """

rule multiallelism_SCE:
 input:
  VARIANT_BED = "ma/{strain}/{nod}/{nod}.variant_bed",
  STRAIN_BED = "bin/{strain}_chr.bed",
  ALLELE_COUNT_BED = "ma/{strain}/{nod}/{nod}.allele_count_bed"
 output:
  PLOT_MA = "ma/{strain}/{nod}/{nod}.ma_plot.pdf",
  FNOD_BED = "ma/{strain}/{nod}/{nod}.fnod_bed",
  SCE_BED = "ma/{strain}/{nod}/{nod}.sce_bed"
 params:
  PATH_TO_FNOD = config["PATH_TO_FNOD"],
  PATH_TO_DRCR = config["PATH_TO_DRCR"]
 shell:
  """ bin/MA_SCE_plot.sh bin/{wildcards.strain}_chr.bed {wildcards.strain} {wildcards.nod} {params.PATH_TO_FNOD} {params.PATH_TO_DRCR}"""

rule calc_timing_props:
 input:
  TUMOUR_INFO = "/hps/nobackup2/flicek/user/cander21/lce20191105/samples/allStrainsDRIVERS.LCE_isdWGSs.tab",
  SCE_MA2 = "ma/{strain}/{nod}/{nod}.sce_bed",
 output:
  "ma/{strain}/{nod}/{nod}.timing"
 shell:
  """ bin/calc_timing_props.sh {wildcards.strain} {wildcards.nod} {input.TUMOUR_INFO} """

rule plot_driver_times:
 input:
  expand("ma/{strain}/{nod}/{nod}.timing", nod=samples , strain=config["STRAIN"])
 output:
  "ma/{strain}/{strain}.timing_plot.pdf"
 shell:
  """ cat {input} | sed 's/ /\t/g' | Rscript --vanilla bin/plot_driver_times.R {wildcards.strain} """
