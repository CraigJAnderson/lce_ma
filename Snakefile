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
  MPILEUP_OUT = expand("ma/{strain}/{nod}/{nod}.mpileup", nod=samples, strain=config["STRAIN"])
  
ruleorder : filter_variants > collect_control_coordinates > genotype

rule filter_variants:
 output: 
  "geno/{strain}/{nod}_arc_filtered.nodMat"
 params:
  config["VARIANTS"]
 shell:
  """grep -e ',{wildcards.nod},' {params} | awk -F, '{{if (($19 == "") && ($20 == 0) && ($34 == "" || $34 > 1e-4 ) && ($10 != "")) print $0}}' > {output} """

rule collect_control_coordinates:
 input:
  expand("geno/"+"{strain}"+"/{nod}_arc_filtered.nodMat", nod=samples, strain=config["STRAIN"])
 output:
  CONTROL_POS = "ma/{strain}/{strain}_den.arc_filtered_control_sites.pos",
  CONTROL_BED = "ma/{strain}/{strain}_den.arc_filtered_control_sites.bed",
  CONTROL_BED_GZ = "ma/{strain}/{strain}_den.arc_filtered_control_sites.bed.gz"
 shell:
  """ cat {input} | awk '{{print $1"\t"$2}}' FS="[_,]" | sort -k 1,1 -k2,2n -u > {output.CONTROL_POS} ; """
  """ cat {input} | awk '{{print $1"\t"$2"\t"$2+1}}' FS="[_,]" | sort -k 1,1 -k2,2n -u > {output.CONTROL_BED} ; """
  """ bgzip -c {output.CONTROL_BED} > {output.CONTROL_BED_GZ} """

rule genotype:
 input:
  CONTROL_POS = "ma/{strain}/{strain}_den.arc_filtered_control_sites.pos",
 params:
  GENOME = config["GENOME"],
  BAM = lambda wildcards : sample_linked_bams[wildcards.nod]
 output:
  "ma/{strain}/{nod}/{nod}.mpileup"
 shell:
  """samtools mpileup {params.BAM} -l {input.CONTROL_POS} -f {params.GENOME} -Q 20 --skip-indels --ff DUP > {output} """
