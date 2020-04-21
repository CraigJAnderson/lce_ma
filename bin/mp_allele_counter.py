##python 3
##mp_allele_counter.py
##mpileup multiallelic genotyper, takes cleaned up mpileup output and designates a site as MA or BA, before outputting a bed file with the class of variant and genotype counts for each nucleotide
import sys
import pandas as pd
import collections
from collections import Counter
import numpy as np

#this makes any print assignments for debugging easy to read, then the next disables a complaint about copying data, which we don't mind
pd.set_option('display.max_rows', 5)
pd.options.mode.chained_assignment = None

#import mp file
colnames = ['a','b','c','d','e']
reads = pd.read_csv(sys.argv[1],sep="\t", names=colnames, header=None, dtype = {"a":object, "b":np.int64, "c":object, "d":np.int64, "e":object})

reads['countvars']=reads.e.str.replace(',', '').str.upper().apply(collections.Counter)
reads['ref']= reads.e.str.count(r'[,]')
reads = reads.reset_index()

##remove empties and those with very low alt alleles
keep_list = []
for x, row in reads.iterrows():
 if (len(row.countvars.values()) != 0) :
  if (sum(row.countvars.values()) >= 2) :
   keep_list.append(x)
  else :
   continue
 else :
  continue

reads = reads[reads.index.isin(keep_list)].reset_index(drop=True)

##apply threshold of minimum reads:
thresh = int(2)
for x, row in reads.iterrows():
 if reads['countvars'].iloc[x]['A'] < thresh:
  del reads['countvars'].iloc[x]['A']
 if reads['countvars'].iloc[x]['C'] < thresh:
  del reads['countvars'].iloc[x]['C']
 if reads['countvars'].iloc[x]['G'] < thresh:
  del reads['countvars'].iloc[x]['G']
 if reads['countvars'].iloc[x]['T'] < thresh:
  del reads['countvars'].iloc[x]['T']

keep_list = []
for x, row in reads.iterrows():
 if (len(row.countvars.values()) != 0) :
  if (row.ref >= 10) :
   keep_list.append(x)
  else :
   continue
 else :
  continue

reads = reads[reads.index.isin(keep_list)].reset_index(drop=True)

reads = pd.concat([reads, reads['countvars'].apply(pd.Series)], axis=1)
reads = reads.fillna(0)
if 'A' not in reads:
 reads['A'] = float(0.0)

if 'C' not in reads:
 reads['C'] = float(0.0)

if 'G' not in reads:
 reads['G'] = float(0.0)

if 'T' not in reads:
 reads['T'] = float(0.0)

reads.c = reads.c.str.title()
reads = reads.rename(columns = {'c':'REF'})
reads = reads.reset_index()
df = reads[['a', 'b', 'ref', 'REF', 'A', 'C', 'G', 'T']]

####transfer ref count to ref base
for x,y in df.iterrows():
 tmp = df['REF'].iloc[x]
 df[tmp].iloc[x] = df['ref'].iloc[x]

df['type'] = "Nan"
for x, row in df.iterrows():
 nuc_list = ["A","C","G","T"]
 nuc_list.remove(df.REF[x])
 ##annotate multiallelic
 if ((df.loc[x,nuc_list] == 0.0) .sum()) < 2 :
  df.at[x,'type'] = "MA"
 elif ((df.loc[x,nuc_list] == 0.0) .sum()) == 2 :
  df.at[x,'type'] = "BA"
 else:
  continue

df = df[df['type'].str.contains("MA|BA")]

keep_list = []
for x, row in df.iterrows():
 nuc_list = ["A","C","G","T"]
 nuc_list.remove(df.REF[x])
 ##keep if VAF is >0.025
 if (((df.loc[x,nuc_list]).sum())/(((df.loc[x,nuc_list]).sum())+(df.loc[x,df.REF[x]]).sum()))> 0.025 :
  keep_list.append(x)
 else:
   pass

df = df[df.index.isin(keep_list)].reset_index(drop=True)
df[['A','C','G','T']] = df[['A','C','G','T']].astype(np.int64)
##outputbed
df['nuc'] = df["type"].map(str) + "," + df["A"].map(str) + "," +df["C"].map(str)+ "," +df["G"].map(str)+ "," +df["T"].map(str)
df.astype({'b': 'int64'}).dtypes
df["c"] = df.b+1
df = df[['a', 'b', 'c', 'nuc']]

from io import StringIO
output = StringIO()
df.to_csv(output,sep="\t",header=False, index=False)
output.seek(0)
print(output.read())
##end_script
