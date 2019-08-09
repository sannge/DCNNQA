import os

#The script has a method that will combine all the output files
#the called feature scripts generated.
def combine(pdb,dir_output):
    files=os.listdir(dir_output)
    w=open(dir_output+'/'+'allresults.txt','w')
    title='model'
    scores=[]


    for file in files:
        #get the feature name by splitting the file names
        #make sure you give the output file name right.
        if not file=='allresults.txt':
            if not os.path.isdir(dir_output+file):
                if file.endswith('.txt'):
                    featureName=file.split('_')[0]
                    title=title+' '+featureName
                    r=open(dir_output+'/'+file,'r')
                    all=r.readlines()
                    for one in all:
                        one=one.rstrip('\n')
                        names=one.split()
                        if len(names)==1:
                            score=names[0]
                            scores.append(score)
                        elif len(names)>1:
                            name=names[0]
                            score=names[1]
                            scores.append(score)




    w.write(title+'\n')
    nums=' '.join(scores)
    w.write(pdb+' '+ nums)
    print(len(scores))
