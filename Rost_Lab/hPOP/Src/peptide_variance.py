import csv
import numpy as np
import pandas as pd

with open('hpop2_aligned.mat.csv') as hpop2:
	hpop_reader = csv.reader(hpop2)
	data = list(hpop_reader)

# pv = {}	

# for row in reader:
# 	for i in 
# 		pv[row[0]] = 
	#turn into some array
	# convert to array of expression
	# sum up top 3 peptide
	# find variance across array

#clean up Na before = turn in 

peptide_name = []
protein_name = []
expression = []
lrow = len(data[0])-3

for row in data[:1]:
	peptide_name.append(row[0])
	protein_name.append(row[1])
	exp = []
	exp.append(row[2:lrow])
	expression.append(exp)

hpop_array = np.zeros(len(data)-1, dtype = {'names':('peptide', 'protein', 'expression'), 'formats':('U15', 'U15', (1,lrow-2))})

hpop_array['peptide'] = peptide_name
hpop_array['protein'] = protein_name
hpop_array['expression'] = expression
print(hpop_array[0])



