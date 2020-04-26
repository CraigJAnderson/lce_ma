###script start###
##script to take mao file and calculate the numbers of T>N and C>N sites within each SCE boundary
## originally in /hps/nobackup2/flicek/user/cander21/LCE20190618/bin/sce_sensitive_ma_prop.py
##python 3
import sys
import pandas as pd
from io import StringIO
import math
pd.options.mode.chained_assignment = None

df = pd.read_csv(sys.argv[1],sep="\t",header=None)
df.columns = ["chr","start","end","REF","mut","A","C","G","T","sce"]
df = df.sort_values(['chr', 'sce', 'start'], ascending=[True, True, True])
df = df.reset_index()
df = df[['chr','start','sce','mut','REF','A','C','G','T']]

##this section calls sites as multiallelic or biallelic, regardless of coverage (so 1X can be MA)
df['type'] = "Nan"
for x, row in df.iterrows():
 nuc_list = ["A","C","G","T"]
 nuc_list.remove(df.REF[x])
 ##annotate multiallelic
 if ((df.ix[x,nuc_list] == 0) .sum()) < 2 :
  df.at[x,'type'] = "MA"
 else:
  df.at[x,'type'] = "BA"

##this part starts by folding half of the variants so there are only C>N and T>N. It then defines the number of each within each sister chromatid exchange region... ####important#### I need it to do this for MA and BA sites...
df = df[['chr','start','sce','mut','type']]
df = df.replace({'/':'_'}, regex=True)
df.columns = ["a","b","c","d","e"]
df = df.replace({'A_T':'T_A'}, regex=True)
df = df.replace({'A_G':'T_C'}, regex=True)
df = df.replace({'A_C':'T_G'}, regex=True)
df = df.replace({'G_T':'C_A'}, regex=True)
df = df.replace({'G_C':'C_G'}, regex=True)
df = df.replace({'G_A':'C_T'}, regex=True)

##split genome into chunks
ws = int(sys.argv[2]) ##e.g. for 10 Mb chuncks 10000000
df['g'] = 0
df['f'] = 0
for a,b in df.iterrows():
 df['g'].iloc[a] = int(math.floor(df['b'].iloc[a]/(ws))*ws) ##ws/2 originally
 df['g'] = df['g'].astype(int)
 df['f'].iloc[a] = str(df['c'].iloc[a])+":"+str(df['g'].iloc[a]) ##df.at[x,'ALT_MAX'] = alt_max alternative??

dfma = df[df['e'].str.contains("MA")]
dfma = dfma.reset_index()
allcgma=[]
allc2nma = []
allt2nma = []
C2Nma = []
T2Nma = []
for x,y in dfma.iterrows():
 if dfma['f'][x] == dfma['f'].shift(-1)[x]:
  if dfma['d'][x] == 'C_T' or dfma['d'][x] == 'C_G' or dfma['d'][x] == 'C_A':
   C2Nma.append(x)
  elif dfma['d'][x] == 'T_C' or dfma['d'][x] == 'T_G' or dfma['d'][x] == 'T_A':
   T2Nma.append(x)
 elif dfma['f'][x] != dfma['f'].shift(-1)[x]:
  if dfma['d'][x] == 'C_T' or dfma['d'][x] == 'C_G' or dfma['d'][x] == 'C_A':
   C2Nma.append(x)
  elif dfma['d'][x] == 'T_C' or dfma['d'][x] == 'T_G' or dfma['d'][x] == 'T_A':
   T2Nma.append(x)
  allt2nma.append(len(T2Nma))
  allc2nma.append(len(C2Nma))
  allcgma.append(dfma['f'][x])
  C2Nma = []
  T2Nma = []
 elif dfma['f'].shift(+1)[x] == 'nan':
  if dfma['d'][x] == 'C_T' or dfma['d'][x] == 'C_G' or dfma['d'][x] == 'C_A':
   C2Nma.append(x)
  elif dfma['d'][x] == 'T_C' or dfma['d'][x] == 'T_G' or dfma['d'][x] == 'T_A':
   T2Nma.append(x)
  allt2nma.append(len(T2Nma))
  allc2nma.append(len(C2Nma))
  allcgma.append(dfma['f'][x])
  exit


dfba = df[df['e'].str.contains("BA")]
dfba = dfba.reset_index()
allcgba=[]
allc2nba = []
allt2nba = []
C2Nba = []
T2Nba = []
for x,y in dfba.iterrows():
 if dfba['f'][x] == dfba['f'].shift(-1)[x]:
  if dfba['d'][x] == 'C_T' or dfba['d'][x] == 'C_G' or dfba['d'][x] == 'C_A':
   C2Nba.append(x)
  elif dfba['d'][x] == 'T_C' or dfba['d'][x] == 'T_G' or dfba['d'][x] == 'T_A':
   T2Nba.append(x)
 elif dfba['f'][x] != dfba['f'].shift(-1)[x]:
  if dfba['d'][x] == 'C_T' or dfba['d'][x] == 'C_G' or dfba['d'][x] == 'C_A':
   C2Nba.append(x)
  elif dfba['d'][x] == 'T_C' or dfba['d'][x] == 'T_G' or dfba['d'][x] == 'T_A':
   T2Nba.append(x)
  allt2nba.append(len(T2Nba))
  allc2nba.append(len(C2Nba))
  allcgba.append(dfba['f'][x])
  C2Nba = []
  T2Nba = []
 elif dfba['f'].shift(+1)[x] == 'nan':
  if dfba['d'][x] == 'C_T' or dfba['d'][x] == 'C_G' or dfba['d'][x] == 'C_A':
   C2Nba.append(x)
  elif dfba['d'][x] == 'T_C' or dfba['d'][x] == 'T_G' or dfba['d'][x] == 'T_A':
   T2Nba.append(x)
  allt2nba.append(len(T2Nba))
  allc2nba.append(len(C2Nba))
  allcgba.append(dfba['f'][x])
  exit


newdf = pd.DataFrame({'cgba': allcgba,'C2Nba': allc2nba,'T2Nba': allt2nba,'cgma': allcgma,'C2Nma': allc2nma,'T2Nma': allt2nma})
newdf['chr'] = newdf['cgba'].str.split('[:_]',expand=True)[0]

##the proportion of sites that are multiallelic in that region
newdf['mapc2n'] = 0.0
for x,y in newdf.iterrows():
 newdf['mapc2n'].iloc[x] = (newdf.C2Nma[x] /(newdf.C2Nba[x] +newdf.C2Nma[x]))

newdf['mapt2n'] = 0.0
for x,y in newdf.iterrows():
 newdf['mapt2n'].iloc[x] = (newdf.T2Nma[x] /(newdf.T2Nba[x] +newdf.T2Nma[x]))


##make new chr column to fill in missing chromosomes as zeros.

chr_list = list(range(1,20))+list('X')

chr_list = pd.DataFrame(chr_list)
chr_list.columns = ['chr']
chr_list["chr"] = chr_list.chr.apply(str)
results = newdf.append(chr_list,ignore_index=True).fillna(0, downcast='integer')
results = results.reset_index()

results.chr = pd.to_numeric(results.chr, errors='coerce')
results = results.sort_values(by=['chr'],ascending=True).fillna('X')

list_of_rows_to_remove = []
for x, row in results.iterrows():
 if results['cgma'][x] == 0 :
  if results['chr'][x] == results['chr'].shift(+1)[x]:
   list_of_rows_to_remove.append(x)


results = results[~results['index'].isin(list_of_rows_to_remove)]
results = results.reset_index()
results = results[['chr', 'C2Nba', 'C2Nma', 'T2Nba', 'T2Nma', 'mapc2n', 'mapt2n','cgba']]
results[['C2Nba','C2Nma','T2Nba','T2Nma']] = results[['C2Nba','C2Nma','T2Nba','T2Nma']].astype(int)
results.chr = pd.to_numeric(results.chr,errors='ignore')

from io import StringIO
output = StringIO()
results.to_csv(output,sep="\t",header=False, index=False)
output.seek(0)
print(output.read())
###endscript###
