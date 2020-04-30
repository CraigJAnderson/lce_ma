#!/bin/sh
# properties = {"rule": "genotype", "local": false, "input": ["ma/c3h/c3h.arc_filtered_control_sites.pos"], "output": ["ma/c3h/91613_N1/91613_N1.mpileup"], "wildcards": ["c3h", "91613_N1"], "params": {"GENOME": "/nfs/leia/research/flicek/user/fnc-odompipe/lce-collaboration/genome_metadata/genomes/mus_musculus/C3H_HeJ_v1_S288C/C3H_HeJ_v1_S288C.fa", "BAM": "/nfs/leia/research/flicek/user/fnc-odompipe/lce-collaboration/repository/merged_alignments/alignments_20171018/do9419_91613_N1.bam"}, "log": [], "threads": 1, "resources": {}, "jobid": 2317, "cluster": {"run_time": "25:00", "memory": "15000", "resources": "\"rusage[mem=15000]\"", "output": "/hps/nobackup2/flicek/user/cander21/dojo/SM/logs/", "error": "/hps/nobackup2/flicek/user/cander21/dojo/SM/logs/"}}
cd /hps/nobackup2/flicek/user/cander21/lce20191105/share && \
/hps/nobackup2/flicek/user/cander21/bin/conda/envs/lce_ma/bin/python -m snakemake ma/c3h/91613_N1/91613_N1.mpileup --snakefile /hps/nobackup2/flicek/user/cander21/lce20191105/share/Snakefile \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /hps/nobackup2/flicek/user/cander21/lce20191105/share/.snakemake/tmp.dhwovhw5 ma/c3h/c3h.arc_filtered_control_sites.pos --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
 --configfile /hps/nobackup2/flicek/user/cander21/lce20191105/share/config_c3h.yaml -p --nocolor \
--notemp --quiet --no-hooks --nolock --printshellcmds  --force-use-threads  --allowed-rules genotype  && touch "/hps/nobackup2/flicek/user/cander21/lce20191105/share/.snakemake/tmp.dhwovhw5/2317.jobfinished" || (touch "/hps/nobackup2/flicek/user/cander21/lce20191105/share/.snakemake/tmp.dhwovhw5/2317.jobfailed"; exit 1)

