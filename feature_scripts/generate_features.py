import os
import sys
import numpy
from combined_output import combine
if len(sys.argv)!=4:
    print("The program should be excecuted as the following!")
    print("python generate_features.py <dir_pdb> <new_dir_output> <dir_root>")
    print("<dir_pdb> - path to single pdb file")
    print("<new_dir_output> - directory of new output folder")
    print("<dir_root> - path to softwareforprediction")
    sys.exit()
first_dir_pdb=sys.argv[1]
first_dir_pdb=os.path.abspath(first_dir_pdb)

dir_output=sys.argv[2]
dir_output=os.path.abspath(dir_output)


#asking the root of the software
general_path=sys.argv[3]
#getting argv[1]
if os.path.exists(dir_output):
    sys.exit('The directory already exists!')

os.mkdir(dir_output)
os.mkdir(dir_output+'/input_dir/')
os.system('cp '+first_dir_pdb+' '+dir_output+'/input_dir/')
dir_pdb=dir_output+'/input_dir/'
#This program will take pdb as an input, and write the features in the output file\\
#in the folder, there must be only one pdb
abs_dir_pdb=''
pdbFileName=''
for r,d,f in os.walk(dir_pdb):
    for file in f:
        if not file.endswith('.pdb'):
            file=file+'.pdb'
        pdbFileName=file
        abs_dir_pdb=dir_pdb+'/'+file



EuclideanOutputFile=dir_output+'/'+'Euclidean_output.txt'
EuclideanCommand="perl Euclidean_compact.pl "+abs_dir_pdb+" "+EuclideanOutputFile


#FOR dssp disperse

#getting general path to software system

general_path=general_path+'/'



DSSPOutputFolder=dir_output+'/'+'DSSP_output'
mkDSSP_path=general_path+'/tools/mkdssp'
datasetDSSP_path=general_path+'/tools/dssp2dataset.pl'
#DSSP only accepts the pdb location of the folder, not the import pdb;
#pdb_location will get the folder of import pdb
#perl dssp2secondary.pl /data/store1/summer2019/DSSP/ /data/store1/summer2019/DSSP/mkdssp
#/data/store1/summer2019/DSSP/dssp2dataset.pl
#/data/store1/summer2019/DSSP/output_dispersd
DSSPCommand="perl dssp2secondary.pl "+dir_pdb+" "+mkDSSP_path+" "+datasetDSSP_path+" "+DSSPOutputFolder
os.system(DSSPCommand)

#For surface output file
SurfaceOutputFile=dir_output+'/'+'Surface_output.txt'
#getting surface score
SurfaceCommand='perl surface_score.pl '+DSSPOutputFolder+' '+SurfaceOutputFile
#For weighted_exposed output file
WExposedOutputFile=dir_output+'/'+'WeightedExposedArea_output.txt'
#getting wexposedarea score
WExposedCommand='perl weight_exposed1.pl '+DSSPOutputFolder+' '+WExposedOutputFile
#getting RWplus unnormalized Output file
unnormalizedFolder=dir_output+'/'+'unnormalized/'
if not os.path.exists(unnormalizedFolder):
    os.mkdir(unnormalizedFolder)
#getting the name of pdb for unnormalized RWplus score
db=pdbFileName.replace(".pdb","")
print("DB is ",db)
RWplusOutputFile=dir_output+'/'+'unnormalized/'+db
#getting the RWplus score
RWplusCommand='perl RWPlus_score.pl '+general_path+'tools/ '+dir_pdb+' '+RWplusOutputFile
#getting pdb 2 fasta output file
pdb2fastaOutputFile=dir_output+'/'+'fastaforpdb'
#getting pdb2fasta
pdb2fastaCommand=general_path+'tools/pdb2fasta '+abs_dir_pdb+' > '+pdb2fastaOutputFile

#getting unnormalized output file for RWPlus
normalizedRWplusOutputFile=dir_output+'/'+'normalizedRWplus_output.txt'
#Normalized RWplus RWplusCommand
normalizedRWplusCommand='perl norm_RWplus.pl '+dir_output+'/fastaforpdb '+dir_output+'/unnormalized/'+db+' '+normalizedRWplusOutputFile
#setting up for STRIDE
#Output file path for STRIDE
if not os.path.exists(dir_output+'/'+'STRIDE/'):
    os.mkdir(dir_output+'/'+'STRIDE/')
StrideOutputFile=dir_output+'/'+'STRIDE/'+'stride_output.txt'
#Command for STRIDE
StrideCommand=general_path+'/tools/STRIDE/stride '+abs_dir_pdb+' > '+StrideOutputFile
#MASS potential outputs
if not os.path.exists(dir_output+'/'+'MASS/'):
    os.mkdir(dir_output+'/'+'MASS/')
MassOutputFiles=dir_output+'/'+'MASS/'+db

MassCommand='perl MASS_potentials.pl . '+abs_dir_pdb+' '+StrideOutputFile+' '+MassOutputFiles
#ASPA output file(unnormalized original version)
ASPAOutputFile=dir_output+'/ASPA_output.txt'
#ASPA command file(unnormalized original version)
ASPACommand1='python2.7 formatASPA.py '+MassOutputFiles+'.ASPA '+ASPAOutputFile
#ASPACommand2='python2.7 formatASPA2.py '+MassOutputFiles+'.ASPA '+ASPAOutputFile
#SBROD score generation and write the output to a file
SBRODOutputFile=dir_output+'/SBROD_output.txt'
SBRODCommand=general_path+'/tools/SBROD/sbrod --scale '+abs_dir_pdb+' > '+SBRODOutputFile
#ASPA output file(unnormalized original version)
DDP_CAOutputFile=dir_output+'/DDPCA_output.txt'
#DDP_CA command file(unnormalized original version)
DDPCACommand1='python2.7 formatCDP.py '+MassOutputFiles+'.DDP_CA '+DDP_CAOutputFile
#DDPCACommand2='python2.7 formatCDP2.py '+MassOutputFiles+'.DDP_CA '+DDP_CAOutputFile
#FOr DDP_CB
DDP_CBOutputFile=dir_output+'/DDPCB_output.txt'
DDPCBCommand1='python2.7 formatCDP.py '+MassOutputFiles+'.DDP_CB '+DDP_CBOutputFile
#For SSDP_CA
SSDP_CAOutputFile=dir_output+'/SSDPCA_output.txt'
SSDPCACommand1='python2.7 formatCDP.py '+MassOutputFiles+'.SSDP_CA '+SSDP_CAOutputFile
#normalized here below commented
#SSDPCACommand2='python2.7 formatCDP2.py '+MassOutputFiles+'.SSDP_CA '+SSDP_CAOutputFile

#For SSDP_CB
SSDP_CBOutputFile=dir_output+'/SSDPCB_output.txt'

SSDPCBCommand1='python2.7 formatCDP.py '+MassOutputFiles+'.SSDP_CB '+SSDP_CBOutputFile
#normalized here below commented
#SSDPCBCommand2='python2.7 formatCDP2.py '+MassOutputFiles+'.SSDP_CB '+SSDP_CBOutputFile
ASPRCAOutputFile=dir_output+'/ASPRCA_output.txt'
ASPRCACommand1='python2.7 formatRSAP.py '+MassOutputFiles+'.ASPR_CA '+ASPRCAOutputFile

ASPRCBOutputFile=dir_output+'/ASPRCB_output.txt'
ASPRCBCommand1='python2.7 formatRSAP.py '+MassOutputFiles+'.ASPR_CB '+ASPRCBOutputFile

CDPCAOutputFile=dir_output+'/CDPCA_output.txt'
CDPCACommand1='python2.7 formatCDP.py '+MassOutputFiles+'.CDP_CA '+CDPCAOutputFile

CDPCBOutputFile=dir_output+'/CDPCB_output.txt'
CDPCBCommand1='python2.7 formatCDP.py '+MassOutputFiles+'.CDP_CB '+CDPCBOutputFile

CSPCAOutputFile=dir_output+'/CSPCA_output.txt'
CSPCACommand1='python2.7 formatRSAP.py '+MassOutputFiles+'.CSP_CA '+CSPCAOutputFile

CSPCBOutputFile=dir_output+'/CSPCB_output.txt'
CSPCBCommand1='python2.7 formatRSAP.py '+MassOutputFiles+'.CSP_CB '+CSPCBOutputFile

PAPOutputFile=dir_output+'/PAP_output.txt'
PAPCommand1='python2.7 formatTAP.py '+MassOutputFiles+'.PAP '+PAPOutputFile

RSAPOutputFile=dir_output+'/RSAP_output.txt'
RSAPCommand1='python2.7 formatRSAP.py '+MassOutputFiles+'.RSAP '+RSAPOutputFile

TAPOutputFile=dir_output+'/TAP_output.txt'
TAPCommand1='python2.7 formatTAP.py '+MassOutputFiles+'.TAP '+TAPOutputFile

VDPCAOutputFile=dir_output+'/VDPCA_output.txt'
VDPCACommand1='python2.7 formatRSAP.py '+MassOutputFiles+'.VDP_CA '+VDPCAOutputFile
#write command using os.system
os.system(EuclideanCommand)
os.system(SurfaceCommand)
os.system(WExposedCommand)
os.system(RWplusCommand)
os.system(pdb2fastaCommand)
os.system(normalizedRWplusCommand)
os.system(StrideCommand)
os.system(MassCommand)
os.system(ASPACommand1)
os.system(SBRODCommand)
os.system(DDPCACommand1)
os.system(DDPCBCommand1)
os.system(SSDPCACommand1)
os.system(SSDPCBCommand1)
os.system(ASPRCACommand1)
os.system(ASPRCBCommand1)
os.system(CDPCACommand1)
os.system(CDPCBCommand1)
os.system(CSPCACommand1)
os.system(CSPCBCommand1)
os.system(PAPCommand1)
os.system(RSAPCommand1)
os.system(TAPCommand1)
os.system(VDPCACommand1)



#getting the name of the pdb to display in the file
combine(pdbFileName,dir_output)

#setting up for prediction
allScores=[]
r=open(dir_output+'/allresults.txt','r')
titleLine=r.readline()
scoreLine=r.readline()
scoreLine=scoreLine.rstrip('\n')
names=scoreLine.split()
pdbName=names[0]
for i in range(1,len(names)-1):
    allScores.append(names[i])


for i in range(2606,3138):
    allScores.append('0')

#Reshaping into 54*54
predict_data=numpy.reshape(allScores,(56,56))
