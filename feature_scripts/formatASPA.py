import os, sys
import numpy as np

inputFile = sys.argv[1]
outputFile=sys.argv[2]

if not os.path.isfile(inputFile):
    print "Input file does not exist"
    exit()

res = {
'1': [],
'2': [],
'3': [],
'4': [],
'5': [],
'6': [],
'7': [],
'8': [],
'9': [],
'10': [],
'11': [],
'12': [],
'13': [],
'14': [],
'15': [],
'16': [],
'17': [],
'18': [],
'19': [],
'20': [],
'21': [],
'22': [],
'23': [],
'24': [],
'25': [],
'26': [],
'27': [],
'28': [],
'29': [],
'30': [],
'31': [],
'32': [],
'33': [],
'34': [],
'35': [],
'36': [],
'37': [],
'38': [],
'39': [],
'40': [],
}

with open(inputFile, "r+") as f:
    aspaFile = []
    for line in f:
        #print line.strip()
        line = line.strip()
        line = line.split()
        del line[0]
        #print line
        aspaFile.append(line)
    f.close()


numbers = []
for i in aspaFile:
    number = i[0]
    if number not in numbers:
        numbers.append(number)

numbers = list(map(int, numbers))
numbers.sort()
numbers = list(map(str, numbers))
#print numbers

combLine = []
for number in numbers:
    for line in aspaFile:
        if number in line:
            combLine.append(line)

#print combLine

i = 0
for lists in combLine:
    if combLine[i][0] in res.keys():
        #print combLine[i][0]
        num = combLine[i][0]
        del combLine[i][0]

        res[num].append(combLine[i])
        i = i + 1

#print res.keys()
#print res.values()
w=open(outputFile,'w')

for key, val in res.items():
    if len(val) == 0:

        w.write('0'+'\n')

    else:

        val = np.array(val, dtype=np.float64)
        val = list(np.array(val).sum(axis=0)/len(val))
        #print val

        w.write(' '.join(map(str, val))+'\n')
w.close()
#print(len(res))
