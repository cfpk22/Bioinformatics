# CRE motif genome search for GC content and logo of flanking residues
# Created by Conor Kelly 29/4/20. Test Conor Kelly 29/4/20
# Libraries required to run code, require logomaker and newer version of matplotlib.pyplot library installed 
# system arguments are 1) text file fo DNA data to install, 2) number of flanking residues that you want to extract, 3) exported file string for exported files (usually: chr1_50bp for eg.)

import sys
import pandas as pd
import numpy as np
import regex as re
import matplotlib  
matplotlib.use('TkAgg')   
import matplotlib.pyplot as plt  
import logomaker as lm

with open(sys.argv[1], 'r') as file:
	dna_lower = file.read().replace('\r\n', '')
dna = dna_lower.upper()

D = int(sys.argv[2]) # number for the flanking residues
# Open list for start of each CRE motif
CRE_reading_Start = []
CRE_antisense_Start = []
run_start = []
run_start_2 = []
# Start position of each CRE motif
runs = re.finditer(r"TGACGTCA", dna)
runs_anti = re.finditer(r"ACTGCAGT", dna)
for match in runs: 
	run_start = match.start()
	CRE_reading_Start.append(run_start)
	
for match in runs_anti:
	run_start_2 = match.start()
	CRE_antisense_Start.append(run_start_2)

# print(len(CRE_reading_Start))	
# print(len(CRE_antisense_Start))

df_CRE_start = pd.DataFrame(CRE_reading_Start, columns=['Start'])
df_antisenseCRE_start = pd.DataFrame(CRE_antisense_Start, columns=['Start'])
# print(df_reading_start)
# print(df_antisense_start)
df_CRE_end = df_CRE_start + 8
df_antisenseCRE_end = df_antisenseCRE_start + 8

df_CREsense_3prime_flank = df_CRE_start + D
df_CREsense_5prime_flank = df_CRE_end - D

df_CREantisense_3prime_flank = df_antisenseCRE_start + D
df_CREantisense_5prime_flank = df_antisenseCRE_end - D

# Using the start and end of each region that we want to extract from above, code to do this and end up with a list of DNA sequences
# CRE sense orientation and 3 orime flanking residues
df_CREsense_3prime_flank_final = pd.concat([df_CRE_start, df_CREsense_3prime_flank.rename(columns={"Start":'End'})], axis=1)
list_CREsense_3prime_flank = []

for row in df_CREsense_3prime_flank_final.itertuples():
	start_string = row.Start
	end_string = row.End	
	splicestring = dna[int(start_string):(end_string)]
	splicestring.replace("\r\n", "")
	# print(len(splicestring))
# 	print(splicestring)
	list_CREsense_3prime_flank.append(splicestring)

list_CREsense_3prime_flank = [x.replace("\r\n", "") for x in list_CREsense_3prime_flank] # removes the artefact \r\n that appears in each of the scripts

# CRE sense orientation and 5 orime flanking residues
df_CREsense_5prime_flank_final = pd.concat([df_CREsense_5prime_flank, df_CRE_end.rename(columns={"Start":'End'})], axis=1)
list_CREsense_5prime_flank = []

for row in df_CREsense_5prime_flank_final.itertuples():
	start_string = row.Start
	end_string = row.End	
	splicestring = dna[int(start_string):(end_string)]
# 	print(len(splicestring))
	list_CREsense_5prime_flank.append(splicestring)

list_CREsense_5prime_flank = [x.replace("\r\n", "") for x in list_CREsense_5prime_flank]	
# print(list_CREsense_5prime_flank)
# CRE antisense orientation and 5 orime flanking residues
df_CREantisense_5prime_flank_final = pd.concat([df_CREantisense_5prime_flank, df_antisenseCRE_end.rename(columns={"Start":'End'})], axis=1)
list_CREantisense_5prime_flank = []

for row in df_CREantisense_5prime_flank_final.itertuples():
	start_string = row.Start
	end_string = row.End	
	splicestring = dna[int(start_string):(end_string)]
# 	print(len(splicestring))
	list_CREantisense_5prime_flank.append(splicestring)

list_CREantisense_5prime_flank_temp = [x.replace("\r\n", "") for x in list_CREantisense_5prime_flank]	
list_CREantisense_5prime_flank = [x[::-1] for x in list_CREantisense_5prime_flank_temp] # Reverses the string to enable this to align to the sense strands for logomaker
# CRE antisense orientation and 3 orime flanking residues
df_CREantisense_3prime_flank_final = pd.concat([df_antisenseCRE_start, df_CREantisense_3prime_flank.rename(columns={"Start":'End'})], axis=1)
list_CREantisense_3prime_flank = []

for row in df_CREantisense_3prime_flank_final.itertuples():
	start_string = row.Start
	end_string = row.End	
	splicestring = dna[int(start_string):(end_string)]
# 	print(len(splicestring))
	list_CREantisense_3prime_flank.append(splicestring)

list_CREantisense_3prime_flank_temp = [x.replace("\r\n", "") for x in list_CREantisense_3prime_flank]
list_CREantisense_3prime_flank = [x[::-1] for x in list_CREantisense_3prime_flank_temp]
# Copy of the code over from v1 - where the use of regular expressions unveiled the lack of overlapping matches
df_5prime = pd.DataFrame(list_CREsense_5prime_flank, columns=['DNA'])
df_3prime = pd.DataFrame(list_CREsense_3prime_flank, columns=['DNA'])
df_5prime_f = pd.DataFrame(list_CREantisense_5prime_flank, columns=['DNA'])
df_3prime_f = pd.DataFrame(list_CREantisense_3prime_flank, columns=['DNA'])
df_reading_1 = df_5prime.append(df_3prime)
df_reading_2 = df_reading_1.append(df_5prime_f)
df_reading = df_reading_1.append(df_3prime_f)

split_data = df_reading.DNA.str.split('(\w)', expand=True)
split_data2 = split_data.ix[:, 1::2]
# print(split_data2)

# Convert the dataframe iteratively to a list 
list_dna_data = split_data2.values.tolist()
# print(list_dna_data)

# Calculate GC content of each list item and create new list of results
length_list_dna_data = len(list_dna_data)

list_df_frac = []

for i in range(length_list_dna_data):
	GC_count = list_dna_data[i].count("C") + list_dna_data[i].count("G")
	GC_frac = float(GC_count)/len(list_dna_data[i])
	GC_frac_final = 100*GC_frac
	list_df_frac.append(GC_frac_final)

df_GC_frac = pd.DataFrame(list_df_frac, columns=['GC_content'])	
summary_df_GC_frac = df_GC_frac[["GC_content"]].describe()
summary_df_GC_frac.to_csv(sys.argv[3]+'_GCcontent_summary.csv', index=True, header=True)

# Make a histogram of the data and save this.
plt.figure(1)
plt.hist(df_GC_frac['GC_content'], bins = 20, color='orange', edgecolor='black')
plt.xlabel('GC content (%)')
plt.ylabel('Frequency')
plt.savefig(sys.argv[3]+'_GC_content_histogram.png')

# Make the DNA propensity logo of each of the sites.
# Will append each of the lists from the above to one list below (only taking into account each 5' and 3' flanking residues)
list_CREsense_5prime_flank.extend(list_CREantisense_3prime_flank)
list_CREsense_3prime_flank.extend(list_CREantisense_5prime_flank)
# print(list_CREsense_5prime_flank)


counts_5prime_flank = lm.alignment_to_matrix(list_CREsense_5prime_flank)
# print(counts_5prime_flank)
counts_3prime_flank = lm.alignment_to_matrix(list_CREsense_3prime_flank)
counts_5prime_flank.head()
counts_3prime_flank.head()
plt.figure(2)
logo_5prime = lm.Logo(counts_5prime_flank, stack_order='small_on_top')
logo_5prime.ax.set_xlabel('Position')
logo_5prime.ax.set_ylabel('Counts')
plt.savefig(sys.argv[3]+'_countslogo_5prime.png')
# plt.show()
plt.figure(3)
logo_3prime = lm.Logo(counts_3prime_flank, stack_order='small_on_top')
logo_3prime.ax.set_xlabel('Position')
logo_5prime.ax.set_ylabel('Counts')
plt.savefig(sys.argv[3]+'_logo_3prime.png')
# plt.show()