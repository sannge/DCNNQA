import os, sys
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

#comment out whole block of code w ctrl + /
# resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
# 	    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
# 	    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
# 	    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

inputFile = sys.argv[1]
outputFile=sys.argv[2]

if not os.path.isfile(inputFile):
    print "Input file does not exist"
    exit()

res = {
'ALA' : [],'CYS' : [],'ASP' : [],'GLU' : [],'PHE' : [],'GLY' : [],'HIS' : [],'ILE' : [],\
'LYS' : [],'LEU' : [],'MET' : [],'ASN' : [],'PRO' : [],'GLN' : [],'ARG' : [],'SER' : [],\
'THR' : [],'VAL': [],'TRP' : [], 'TYR' : [] }

with open(inputFile, "r+") as f:
    file = []
    for line in f:
        line = line.strip().split()
        file.append(line[1] + ' ' + line[2])
        #print line[1] + ' ' + line[2]
    f.close()


letters = []
for i in file:
    i = i.split()
    letter = i[0]

    if letter not in letters:
        letters.append(letter)
#print(letters)

combLine = []
for letter in letters:
    for line in file:
        if letter in line:
            combLine.append(line.split())
#print combLine

i = 0
for lists in combLine:
    if combLine[i][0] in res.keys():
        let = combLine[i][0]
        del combLine[i][0]
        res[let].append(combLine[i])
    i = i + 1
w=open(outputFile,'w')
for key, val in res.items():
    if not res[key]:
	res[key] = 0
        w.write(str(key)+' '+'0'+'\n')
    else:
        val = np.array(val, dtype=np.float64)
#	val_normed = (val - val.min(0)) / val.ptp(0)
        val_normed = list(np.array(val).sum(axis=0)/len(val))
	val_normed = np.nan_to_num(val_normed)
        w.write(str(key)+' '+' '.join(map(str, val_normed))+'\n')
