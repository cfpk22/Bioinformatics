#Script for translating CcpNmr TALOS exported data to SSP format for secondary structure predictions.
#Created Conor Kelly - 200304. Tested 200304.
import pandas as pd
import numpy as np
import sys
#need to add the column headers to the TALOS file. Res Name Atom Chemshift separated by a single space

#This command enables you to enter the .txt file with the chemical shifts in the command line as an argument.
file_name = sys.argv[1]
data = pd.read_csv(file_name, sep='\s+',header=0) #reads in the data to columns and rows

# Testing a for loop for doing the data extraction
# Creating a list of the different atoms in the table
uniqueatomlist = data.Atom.unique()

#Attempt at creating for loop to do all the work below slightly quicker? 
for i in uniqueatomlist:
	shift = data[data.Atom == i]
	shift_final = shift.drop(['Name', 'Atom'], axis=1)
	shift_final.to_csv(i+'_shift.txt', sep='\t', index=False, header=False)



