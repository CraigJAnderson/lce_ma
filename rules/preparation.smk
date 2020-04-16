from pathlib import Path
ruleorder: setup_file_system > filter_variants

rule setup_file_system:
 input:
  SAMPLE_LIST = "/hps/nobackup2/flicek/user/cander21/dojo/c3h_bam_name.txt", ##list of samples in variant files and bams, repsectively.
  WRK_DIR = "/hps/nobackup2/flicek/user/cander21/dojo/SM" ##work dir
 run:
  import os
  os.chdir({WRK_DIR})
  dir_tree = ['ma','geno','geno/c3h','geno/ma']
  for dir in dir_tree:
   os.mkdir(dir)
  filepath= {SAMPLE_LIST}
  with open(filepath) as names:
   line = names.readline()
   count = 1
   while line:
    os.mkdir(nod_dir)
    count += 1

rule filter_variants:
 input: 
  "/hps/nobackup2/flicek/user/cander21/lce20191105/geno/c3h.projectedFilteredMu.tab"
 output:
  "geno/"+config["STRAIN"]+"/{nod}-arc-filtered.nodMat"
                                           
#   GENO_STRAIN_FIL = "{G_S}{nod}-arc-filtered.nodMat"
#  CONTROL_POS = "ma/{STRAIN}/{STRAIN}_den.arc_filtered_control_sites.pos",
#  CONTROL_BED = "ma/{STRAIN}/{STRAIN}_den.arc_filtered_control_sites.bed",
#  CONTROL_BED_GZ = "{CONTROL_BED}.gz"
 params:
  NOD = "{nod}"
 shell:
  "cat {input} | grep -e \",{params.NOD},\" | awk -F, '{{if (($19 == \"\") && ($20 == 0) && ($34 == \"\" || $34 > 1e-4 ) && ($10 != \"\")) print $0}}' > {output} ;" ##EXPAND ON DETAILS
#  "cat {ALL_VARIANTS} | awk '{print $1\"\\t\"$2}' FS=\"[_,]\" | sort -k 1,1 -k2,2n -u > {CONTROL_POS} ;"
#  "cat {ALL_VARIANTS} | awk '{print $1\"\\t\"$2\"\\t\"$2+1}' FS=\"[_,]\" | sort -k 1,1 -k2,2n -u > {CONTROL_BED} ;"
#  "bgzip {CONTROL_BED} ;"
#  "tabix {CONTROL_BED_GZ}"
