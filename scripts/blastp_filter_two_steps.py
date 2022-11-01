#Author: Jiang Qian
#Data: 2022.10.05
#modified: 2022.10.18
import argparse,os

parser = argparse.ArgumentParser(description='filter blastp results with identity, length & coverage')
parser.add_argument('-i1', '--parentalDB', required=True, help="parentalDB queried results")
parser.add_argument('-i2', '--targetDB', required=True, help="fused_genesDB queried results")
parser.add_argument('-s', '--speciesOrder', required=True, help="order file of species")
parser.add_argument('-o', '--outputFile', required=True, help="key words of output files")
parser.add_argument('-nl', '--nlidentity', type=int, default=80, help="NL identity")
parser.add_argument('-tx', '--txidentity', type=int, default=80, help="TX identity")
parser.add_argument('-nl2', '--nllength', type=int, default=206, help="NL length")
parser.add_argument('-tx2', '--txlength', type=int, default=156, help="TX length")
parser.add_argument('-tnl1', '--tnlidentity', type=float, default=36, help="TNL identity")
parser.add_argument('-tnl2', '--tnlcov', type=int, default=80, help="TNL coverage")

args = parser.parse_args()

inFile1 = os.path.abspath(args.parentalDB) #parentalDB queried results
inFile2 = os.path.abspath(args.targetDB)
speciesOrder = os.path.abspath(args.speciesOrder)
outFile = os.path.abspath(args.outputFile)
outFile1 = outFile + '.filter1' #filter TNL with NL & TX cutoff
outFile2 = outFile + '.filter2' #filter non-new fused TNL
outFile3 = outFile + '.common' #filter in-group common genes
nlPer = args.nlidentity
txPer = args.txidentity
nlLen = args.nllength
txLen = args.txlength
tnlPer = args.tnlidentity
tnlCov = args.tnlcov

#filter nl & tx candidates separately
def processLine(lineData, parentalDict):
    query = lineData[0]
    subject = lineData[1]
    parentalType = subject.split('|')[2]
    identity = float(lineData[2])
    length = float(lineData[5])
    if (parentalType == "NL"):
        if (identity >= nlPer and length >= nlLen):
            parentalDict['NL'].append(subject)
    elif (parentalType == "TX"):
        if (identity >= txPer and length >= txLen):
            parentalDict['TX'].append(subject)
    return parentalDict

#filter NL & TX candidates separately, and print all hits under that NL & TX ID
df1 = open(inFile1,'r') #parentalDB queried results
parentalCandidates = {}
for line in df1:
    line = line.strip()
    dataInLine = line.split('\t')
    if len(dataInLine) < 14:
        continue

    query = dataInLine[0]
    if query not in parentalCandidates.keys():
        parentalCandidates[query] = {'NL':[],'TX':[]}
    parentalCandidates[query] = processLine(dataInLine,parentalCandidates[query])
df1.close()

#filter TNL with NL & TX cutoff
df1 = open(inFile1,'r') #parentalDB queried results
of1 = open(outFile1,'w')
for line in df1:
    line = line.strip()
    if line.startswith('#'):
        continue

    dataInLine = line.split('\t')
    id1 = dataInLine[0]
    if parentalCandidates[id1]['NL'] != [] and parentalCandidates[id1]['TX'] != []:
        id2 = dataInLine[1]
        if id2 in parentalCandidates[id1]['NL'] or id2 in parentalCandidates[id1]['TX']:
            print(line,file=of1)
of1.close()
df1.close()

#get order of species
sOrder = {}
sf = open(speciesOrder,'r')
for line in sf:
    line = line.strip()
    dataInLine = line.split()
    species = dataInLine[0]
    order = float(dataInLine[1])
    sOrder[species] = order
sf.close()

#filter non-new fused TNL
df2 = open(inFile2,'r')
nonNew = []
paralogy = {}
commonGenes = {}
for line in df2:
    line = line.strip()
    dataInLine = line.split('\t')
    if len(dataInLine) < 14:
        continue

    query = dataInLine[0]
    querySpecies = query.split('|')[0]
    subject = dataInLine[1]
    subjectSpecies = subject.split('|')[0]
    identity = float(dataInLine[2])
    coverage = float(dataInLine[3])

    if coverage < tnlCov: #filter genes with cov < tnlCov
        continue
    if sOrder[querySpecies] == sOrder[subjectSpecies] and querySpecies != subjectSpecies:
        if query not in paralogy.keys():
            paralogy[query]={}
            paralogy[query]={'sub':subject,'ide':identity,'cov':coverage}
        else:
            if identity > paralogy[query]['ide']:
                paralogy[query]={'sub':subject,'ide':identity,'cov':coverage}
               
df2.seek(0)
for line in df2:
    line = line.strip()
    dataInLine = line.split('\t')
    if len(dataInLine) < 14:
        continue

    query = dataInLine[0]
    querySpecies = query.split('|')[0]
    subject = dataInLine[1]
    subjectSpecies = subject.split('|')[0]
    identity = float(dataInLine[2])
    coverage = float(dataInLine[3])

    if coverage < tnlCov: #filter genes with cov < tnlCov
        continue
    if sOrder[querySpecies] > sOrder[subjectSpecies]: #query genes in acestor
        if identity < tnlPer: #filter genes with low identity
            continue
        if query in paralogy.keys():
            if identity >= paralogy[query]['ide']:
                nonNew.append(query)
        else:
            nonNew.append(query)
    elif sOrder[querySpecies] == sOrder[subjectSpecies]: #query genes in in-group
        if querySpecies != subjectSpecies and identity >= 80 and coverage >= 80: #query genes in related species with homology
            if query not in commonGenes.keys():
                commonGenes[query]={}
            if subject not in commonGenes[query].keys():
                commonGenes[query][subject]={'ide':identity,'cov':coverage}
df2.close()

#final output
of1 = open(outFile1,'r')
of2 = open(outFile2,'w')
fusedGenes=[]
for line in of1:
    line = line.strip()
    dataInLine = line.split('\t')
    tnl = dataInLine[0]
    if tnl not in nonNew:
        print(line,file=of2)
        if tnl not in fusedGenes:
            fusedGenes.append(tnl)
    '''print TNLs in out-group
    else:
        if tnl not in paralogy.keys():
            print(tnl)
        else:
            print(tnl+'\t'+str(paralogy[tnl]['sub'])+'\t'+ str(paralogy[tnl]['ide'])+'\t'+str(paralogy[tnl]['cov']))
    '''
of1.close()
of2.close()
#filter common genes
#find common fused genes in-group
commonfusedGenes=[]
for g1 in fusedGenes:
    if g1 in commonGenes.keys():
        commonfusedGenes.append(g1)
        
#output
of2 = open(outFile2,'r')
of3 = open(outFile3,'w')
for line in of2:
    line = line.strip()
    dataInLine = line.split('\t')
    tnl = dataInLine[0]
    if tnl in commonfusedGenes:
        print(line,file=of3)
of3.close()


