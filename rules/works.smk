samples = ["93133_N1", "89985_N1"]

rule all_filter_variants:
 input:
  expand("{nod}.txt", nod=samples)

rule filter_variants:
 output:
  "{nod}.txt"
 shell:
  "grep {wildcards.nod} c3h.projectedFilteredMu.tab > {output} "
