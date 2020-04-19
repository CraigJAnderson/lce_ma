ruleorder : filter_variants > make_genotpe_coords

import pandas as pd
sample_names = pd.read_table(config["SAMPLE_LIST"],sep=" ",header=None)
sample_names.columns = ['bam','nod','NA']
samples= list(sample_names.nod)

rule filter_variants:
 output:
  GENO_STRAIN_FIL = "geno/{wildcards.strain}/{wildcards.nod}_arc_filtered.nodMat"
 shell:
  """grep -e ",{{wildcards.nod}," config["VARIANTS"] | awk -F, '{if (($19 == "") && ($20 == 0) && ($34 == "" || $34 > 1e-4 ) && ($10 != "")) print $0}}' > {output.GENO_STRAIN_FIL}"""

rule make_genotpe_coords: ##why necessary? So controls are derived with respect to the sample list
 input : 
  expand("geno/"+"{strain}"+"/{nod}_arc_filtered.nodMat", nod=samples, strain=config["STRAIN"])
 output : 
  "geno/c3h/filtered.nodMat"
 shell:
  "cat {input} >> {output}"


#  CONTROL_POS = "ma/{strain}/{strain}_den.arc_filtered_control_sites.pos",
#  CONTROL_BED = "ma/{strain}/{strain}_den.arc_filtered_control_sites.bed",
#  CONTROL_BED_GZ = "ma/{strain}/{strain}_den.arc_filtered_control_sites.bed.gz"
# shell:
#  """ cat {input} | awk '{{print $1"\t"$2}}' FS="[_,]" | sort -k 1,1 -k2,2n -u > {CONTROL_POS} ; """
#  """ cat {input} | awk '{{print $1"\t"$2"\t"$2+1}}' FS="[_,]" | sort -k 1,1 -k2,2n -u > {CONTROL_BED} ; """
#  """ bgzip {CONTROL_BED} """
