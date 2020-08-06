The following was used to calculate multi-allelism in WGS from liver nodules derived from C3H mice exposed to a single dose of diethylnitrosamine: https://www.nature.com/articles/s41586-020-2435-1

How to run on LSF Scheduler

<pre>
bsub -J over_lord -W 72:00 -o log.out -e log.error 'snakemake -pr -j 800 --cluster-config cluster.json --cluster "bsub -W {cluster.run_time} -M {cluster.memory} -R {cluster.resources} -e {cluster.error} -o {cluster.output} "' 
</pre>

The idea is that one can provide a different config file for independent cohorts using the --configfile command, which defaults to config_c3h.yaml.


Output files and fields:
<pre>
<sample>.map
The counts of each of the types of variants in a sample, reporting the total proportion of multi-allelic variants, as well as the null estimate. This is based upon the mean of 100 random draws of variant positions seen in the population, exclusing sample specific variant sites
nodule,MA.count,BA.count,proportion.multiallelic.sites,mean.proportion.control.sites.with.unexpected.additional.allele
<sample> 4584 26989 0.145187 0.00127799

<sample>.allele_count_bed
A file of genomic coordinates that also provide information as to whether a variant is multi-allelic (MA) or biallelic (BA), as well as the counts of each of the alleles.
chromosome,start.bp,end.bp,Annotation[Varianttype.Acount.Ccount.Gcount.Tcount]
1	140510	140511	BA,26,8,0,0

<sample>.MA_combos
Identifies proximal variants to define their haplotypes
Location.and.mutation.class.allele1,Reference.allele1,location.and.mutation.class.allele2,Reference.allele2,bp.distance,combination.code,number.of.combinations,specific.haplotype.counts
1:18994504_A/C  A       1:18994839.0_A/T        A       335.0   011     3       (('A,A', 8), ('C,A', 2), ('A,T', 1))
combination codes are as follows:
100: only evidence of ALT alleles in phase
101: evidence of ALT alleles in phase, as well as the REF1 and ALT2
110: evidence of ALT alleles in phase, as well as the ALT1 and REF2
011: ALT alleles occur out of phase
111: evidence of ALT alleles occurring in phase, as well as out of phase
NOTE- there are other varieties, such as 010, 001 and 000, these reflect instance where we can't equivically place one, the other or both ALT alleles with respect to each other and are therefore not included in subsequent analyses.

<sample>.combos_txt
Grades proximal combinations of markers as those the proportion of those impossible under a classical model of clonal expansion.
Nodule,grade1,grade2,grade3,grade4,grade5,(grade3+4+5/1+2+3+4+5) ##the script should have output the proportion of multiallelic sites on the end, but this isn't happening for some reason.
<sample> 71 14 0 2 5 0.076087
NOTES ON GRADES:
each pair of markers is given a grade defined by a 4 digit code referring to: ALT1+ALT2,REF1+REF2,REF1+ALT2,ALT1+REF2.
Grade 1 = 1,1,0,0
Grade 2 = 1,1,1,0 & 1,1,0,1
Grade 3 = >=1,1,0,>=1 & >=1,1,>=1,0
Grade 4 = 1,1,1,1
Grade 5 = >=1,1,>=1,>=1 or >=1,1,0,0
NOTES: only considers the following mutation haplotypes occuring as pairs of variants between 3 and 150 bp of each other: 100,111,101,110

<sample>.fnod_bed
A genomic coordinate file of sister chromatid exchange locations and their strand asymmetry 
chromosome,start.bp,end.bp,segment.specific.identifier.with.strand.asymmetry.score
1       22780   88690149        1_1:-0.02432676881172

<sample>.sce_bed
as the fnod_bed, but the proportion of multi-allelic variants is calculated for each segment
chromosome,start.bp,end.bp,chromosome,segment.specific.identifier,strand.asymmetry.score,proportion.multi-allelic.sites
1 22780 88690149 1 1 -0.02432676881172 0.0483019

<sample>.sce_bed2
as the sce_bed, but the relative coordinates for plotting are given.
chromosome,start.bp,end.bp,chromosome,segment.specific.identifier,strand.asymmetry.score,proportion.multi-allelic.sites,relative.start.site,relative.end.site
1 22780 88690149 1 1 -0.02432676881172 0.0483019 9.11447e-06 0.0354857

<sample>.timing
the proportion of segments that are multiallelic is used to infer the generation after exposure that a clonal population initiated tumouregenesis. I provide annotations to easily perform a Fisher's exact test to see if it's earlier or later than expected. We set the threshold proportion for calling a segment multi-allelic at 0.04
nodule,proportion.multialleic.segments,timing.subgroup[1=early,0=late],drivers,braf[1=true,0=false],egfr[1=true,0=false],hras[1=true,0=false],kras[1=true,0=false]
<sample> 0.638889 0     Kras;Ctnnb1     0       0       0       1
</pre>
