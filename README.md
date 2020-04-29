How to run on LSF Scheduler

<pre>
bsub -J over_lord -W 72:00 -o log.out -e log.error 'snakemake -pr -j 800 --cluster-config cluster.json --cluster "bsub -W {cluster.run_time} -M {cluster.memory} -R {cluster.resources} -e {cluster.error} -o {cluster.output} "' 
</pre>

The idea is that one can provide a different config file for independent cohorts using the --configfile command, which defaults to config_c3h.yaml.

