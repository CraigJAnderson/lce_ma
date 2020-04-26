###start script###
##make a bed file that is sensitive to SCE boundaries, using the fnod file.
#python 3 
#use as: python fnod_bed.py do12458_90793_N2.fnod
#input is fnod file from Martin
import pandas as pd
import sys
from collections import Counter
from io import StringIO

df = pd.read_csv(sys.argv[1],sep=",",header=0)

chr_list = list(df.chr)
counts = Counter(chr_list)
for x, y in counts.items():
 if y > 1:
  for z in range(1, y + 1):
   chr_list[chr_list.index(x)] = x + "_" + str(z)

df['new_chr'] = chr_list
#create a bed format list of locations
df = df[['chr','cStart','cEnd','new_chr']]

output = StringIO()
df.to_csv(output,index=False,sep="\t",header=False)
output.seek(0)
print(output.read())
###end script###
