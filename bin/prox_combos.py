##This is a script to take all the pairs of markers that lie within 1 kb of each other. It phases each read that covers the positions and makes a record of the genotypes, counting the different combinations
##/hps/nobackup2/flicek/user/cander21/LCE20190618/bin/prox_combos.py
#python 3
import pandas as pd
import sys
import numpy as np
import pandas as pd
import sys
import getopt
import collections
from collections import OrderedDict
from collections import Counter
import statistics
import random
import math

pd.options.mode.chained_assignment = None
reads = pd.read_csv(sys.argv[1],sep="\t",header=None)
reads.columns = ['a','b','c','d','e','g']
reads = reads.sort_values(['a', 'b'], ascending=[True, True])
reads = reads.dropna().reset_index(drop=True)

#get only the read id's with a variant on
alt = []
for x, row in reads.iterrows():
 y = [pos for pos, char in enumerate(row['e']) if char != ","]
 outx=[]
 for items in y:
  outx.append(list(row.g.split(","))[items])
 alt.append(",".join(list(outx)))


#get only the read id's with no variant on
ref = []
for x, row in reads.iterrows():
 y = [pos for pos, char in enumerate(row['e']) if char == ","]
 outy=[]
 for items in y:
  outy.append(list(row.g.split(","))[items])
 ref.append(",".join(list(outy)))

reads['ALT'] = pd.DataFrame(alt)
reads['REF'] = pd.DataFrame(ref)

for x, row in reads.iterrows():
 reads.c[x] = reads.c[x].upper()
 REF = reads.c[x]
 reads.e[x] = reads.e[x].replace(',', REF).upper()

match = reads[:]
match = match.drop([0])
match = match.reset_index()
match.columns = ['index2','a_2','b_2','c_2','d_2','e_2','g_2','ALT2','REF2']
result = pd.concat([reads, match], axis=1)#, join_axes=[reads.index])
result = result.dropna()
keepers = []
for x, row in result.iterrows(): 
 if ((row.b_2 - row.b) <= 1000):
  keepers.append(x)
 else:
  pass

match.columns = ['index2','a_2','b_2','c_2','d_2','e_2','g_2','ALT2','REF2']
result = pd.concat([reads, match], axis=1)#, join_axes=[reads.index])
result = result.dropna()
result = result[result.index.isin(keepers)].reset_index(drop=True)

result['h']= 0
all_result1 = []
for x, row in result.iterrows():
 all_result1 = row.e+row.e_2
 result.h[x] = all_result1

result['i']= 0
all_result2 = []
for x, row in result.iterrows():
 all_result2 = row.g+","+row.g_2
 result.i[x] = str.split(all_result2, ",")

#search over list of result over all result and add positions to those that match to a list
#use that list to pull and concatenate genotypes, adding them to a dictionary
result['j'] = 0
for x, row in result.iterrows():
 z = list(set(result.g[x].split(","))|(set(result.g_2[x].split(","))))
 result.j[x] = z

#Get an index of matching result 
result['read_pairs'] = 0
z = []
for x, row in result.iterrows():
 z = []
 for y in result.j[x]: #goes through each unduplicated value of J
  k=[]
  for i, j in enumerate(result.i[x]):
   if j == y:
    k.append(i)
  if len(k) == 2:
   z.append(k)
 z1 = [tuple(l) for l in z]
 read_pairs =[]
 for number in z1:
  n = [result.h[x][i] for i in number]
  read_pairs.append(tuple(n))
  result.read_pairs[x] = list((','.join(w) for w in read_pairs))

result["ab"] = result["a"].map(str) + ":" +result["b"].map(str)
result["ab_2"] = result["a_2"].map(str) + ":" +result["b_2"].map(str)

double_match = []
for x, row in result.iterrows():
 alt = list(row.ALT.split(","))
 alt2 = list(row.ALT2.split(","))
 double_match.append(len(set(alt).intersection(alt2)))


single_match1 = []
for x, row in result.iterrows():
 alt = list(row.ALT.split(","))
 ref2 = list(row.REF2.split(","))
 single_match1.append(len(set(alt).intersection(ref2)))


single_match2 = []
for x, row in result.iterrows():
 ref = list(row.REF.split(","))
 alt2 = list(row.ALT2.split(","))
 single_match2.append(len(set(ref).intersection(alt2)))


ref_match = []
for x, row in result.iterrows():
 ref = list(row.REF.split(","))
 ref2 = list(row.REF2.split(","))
 ref_match.append(len(set(ref).intersection(ref2)))


pos = result[['ab','ab_2','d','d_2']]
###the df index might need to be copied over from the pos index
#put all the count columns together
df = pd.DataFrame( OrderedDict({'double_match': double_match, 'single_match1': single_match1, 'single_match2': single_match2,'ref_match': ref_match,}))

#this is the per SNV count for read co-occurrences
df2 = pd.concat([pos,df],axis=1)

#generate a binary- presence absence dataset
df2["bd"] = df2["double_match"]
df2["bs1"] = df2["single_match1"]
df2["bs2"] = df2["single_match2"]
df2 = df2.assign(bd=df2.bd.gt(0) * 1)
df2 = df2.assign(bs1=df2.bs1.gt(0) * 1)
df2 = df2.assign(bs2=df2.bs2.gt(0) * 1)

#add column binary count to assess presence absence in a single column
df2['presab'] = [''.join(str(x) for x in y) for y in map(tuple, df2[['bd', 'bs1', 'bs2']].values)]
#combine df2 and result dfs
df2['index'] = (result.index)
result2= pd.concat([df2.reset_index(drop=True),result.reset_index(drop=True)], axis=1)

result2['e']=result2.e.str.replace(',', '').str.upper()
result2['e']=result2['e'].replace('', 'NN')

for x, row in result2.iterrows():
 ##annotate with most prominant ALT
 ref = result2.c[x]
 alt = result2.loc[x,'e'].replace(ref, '')
 result2.loc[x,'e'] = str(Counter(alt).most_common(1))[3]


result2['c'] = list(map(lambda x: x.upper(), result2['c']))
result2["ab"] = result2["a"].map(str) + ":" +result2["b"].map(str) + "_" +result2["c"].map(str) + "/" +result2["e"].map(str)

result2['e_2']=result2.e_2.str.replace(',', '').str.upper()
result2['e_2']=result2['e_2'].replace('', 'NN')
##the problem with this bit is it inserts the REF allele in place same with approximately the same position above
for x, row in result2.iterrows():
 ##annotate with most prominant ALT
 ref = result2.c_2[x]
 alt = result2.loc[x,'e_2'].replace(ref, '')
 result2.loc[x,'e_2'] = str(Counter(alt).most_common(1))[3]

result2['c_2'] = list(map(lambda x: x.upper(), result2['c_2']))
result2["ab_2"] = result2["a_2"].map(str) + ":" +result2["b_2"].map(str)+ "_" +result2["c_2"].map(str) + "/" +result2["e_2"].map(str)

result2 = result2.loc[:,~result2.columns.duplicated()]


result2 = result2[(result2.read_pairs != 0)]
result2['diversity']= 0
result2['diversity_count']= 0
diversity = []
diversity_count = []
for x, row in result2.iterrows():
 diversity = Counter(result2.read_pairs[x])
 result2.diversity[x] = tuple(list(diversity.items()))
 diversity_count = len(diversity)
 result2.diversity_count[x] = diversity_count
##unhash this if you want to get rid of boring stuff
##result = result[(result.diversity_count >1)]

##want distance between markers needed.
result2['dist'] = result2.b_2 - result2.b

result3= result2[['ab', 'c', 'ab_2', 'c_2','dist','presab','diversity_count', 'diversity']]
from io import StringIO
output = StringIO()
result3.to_csv(output,sep="\t",header=False, index=False)
output.seek(0)
print(output.read())
