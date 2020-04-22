##/hps/nobackup2/flicek/user/cander21/LCE20190618/bin/allelic_combination_grader.py
##python 3
##script to count specific classes of proximal mutation combinations and grade the occurrence.
##marker grades explained:
#1001 type 1
#1011 or 1101 type 2
#10(>1)1 or 1(>1)01 type 3
#1111 (or zero zeros) type 4

from io import StringIO
import sys
import pandas as pd
import itertools
df = pd.read_csv(sys.stdin,sep=" ",header=None)
df.columns = ['REF1','REF2','COMBOS']
##read in 3 fields

#create 4 new columns for each new marker type and count the occurences.
df['a_c'] = 0
df['ra_c'] = 0
df['ar_c'] = 0
df['rr'] = 0 

for b,row in df.iterrows():
 nuc_list1 = ["A","C","G","T"]
 R1 = row.REF1
 nuc_list1.remove(R1)
 nuc_list2 = ["A","C","G","T"]
 R2 = row.REF2
 nuc_list2.remove(R2)
  ##then cycle through the two lists and produce a string representing each combination, adding them to a list
  ##then cycle through both lists and append them to a REF+ALT list for each position
  ##Then a last list of both ALTs and both REFs ... is that important to have both refs??
 alt_combos = list(itertools.product(nuc_list1, nuc_list2))
 a_c = []
 for x, y in alt_combos:
  z = ''.join([x,y])
  a_c.append(z)
 refalt_combos = list(itertools.product(R1, nuc_list2))
 ra_c = []
 for x, y in refalt_combos:
  z = ''.join([x,y])
  ra_c.append(z)
 altref_combos = list(itertools.product(nuc_list1, R2))
 ar_c = []
 for x, y in altref_combos:
  z = ''.join([x,y])
  ar_c.append(z)
 rr = []
 z = ''.join([R1,R2])
 rr.append(z)
 ##scan through COMBOS to find matches to any of the for lists
 A_C = []
 RA_C = []
 AR_C = []
 RR = []
 for a in row.COMBOS.split(','):
###change here so that you add specific combinations to the relevant empty list, then measure the length of the list, then flush once that has been appended to the df.
  if a in str(rr):
   RR.append(a)
  if a in str(ar_c):
   AR_C.append(a)
  if a in str(ra_c):
   RA_C.append(a)
  if a in str(a_c):
   A_C.append(a)
 df.at[b,'a_c'] = len(A_C)
 df.at[b,'rr'] = len(RR)
 df.at[b,'ra_c'] = len(RA_C)
 df.at[b,'ar_c'] = len(AR_C)



df['grade'] = 0
for x, row in df.iterrows():
 if (row.a_c ==1) & (row.rr ==1) & (row.ra_c ==0) & (row.ar_c ==0):
  df.at[x,'grade'] = 1
 elif (row.a_c ==1) & (row.rr ==1) & (row.ra_c ==1) & (row.ar_c ==0):
  df.at[x,'grade'] = 2
 elif (row.a_c ==1) & (row.rr ==1) & (row.ra_c ==0) & (row.ar_c ==1):
  df.at[x,'grade'] = 2
 elif (row.a_c >=1) & (row.rr ==1) & (row.ra_c ==0) & (row.ar_c >=1):
  df.at[x,'grade'] = 3
 elif (row.a_c >=1) & (row.rr ==1) & (row.ra_c >=1) & (row.ar_c ==0):
  df.at[x,'grade'] = 3
 elif (row.a_c ==1) & (row.rr ==1) & (row.ra_c ==1) & (row.ar_c ==1):
  df.at[x,'grade'] = 4
 elif (row.a_c >=1) & (row.rr ==1) & (row.ra_c >=1) & (row.ar_c >=1):
  df.at[x,'grade'] = 5
 elif (row.a_c >=1) & (row.rr ==1) & (row.ra_c == 0) & (row.ar_c == 0):
  df.at[x,'grade'] = 5

output = StringIO()
df['grade'].to_csv(output,sep="\n", header=False,index=False)
output.seek(0)
print(output.read())
