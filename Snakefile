#snakefile: /hps/nobackup2/flicek/user/cander21/dojo/SM/snakefile
import os
#import numpy as np

configfile : "config.yaml"

include : "rules/test.smk"

import pandas as pd
sample_names = pd.read_table("c3h_bam_name.txt",sep=" ",header=None)
sample_names.columns = ['bam','nod','NA']
samples= list(sample_names.nod)

rule all:
 input: expand("{nod}.txt", nod=samples)
