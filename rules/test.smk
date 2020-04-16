rule filter_variants:
 output:
  "{nod}.txt"
 shell:
  "grep {wildcards.nod} c3h.projectedFilteredMu.tab > {output} "
