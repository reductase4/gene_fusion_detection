#!/user/bin/python3

#Description: extract hmm names of target genes for target species
#Author: Jiang Qian
#Data: 2021.11.23
#V2: add hmm stat(counts of TIR, NBS, LRR, others)
import argparse
import os
import re

parser = argparse.ArgumentParser(description='extract hmm names of target genes for target species')
parser.add_argument('-i', '--input', required=True, help="id type list")
parser.add_argument('-p', '--pfam', required=True, help="pfam search out")
parser.add_argument('-s', '--short', required=True, help="sepcies short name")
parser.add_argument('-m', '--match', required=True, help="match word")
parser.add_argument('-o', '--output', required=True, help="output dir")
parser.add_argument('-k', '--key', required=True, help="output keyword (e.g., species full name)")
args = parser.parse_args()

species = args.short
#matchWord = args.match
matchWord = (args.match).split(',')
output1 = args.output + os.sep + args.key + "_pfam_out." + "_".join(matchWord)
output2 = args.output + os.sep + args.key + "_pfam_stat." + "_".join(matchWord)

Input = open(args.input,'r')
IDs = {}
for line in Input:
    line = line.strip()
    geneID = line.split('\t')[0]
    geneType = line.split('\t')[1]
    #if re.match(matchWord,geneType):
    if geneType in matchWord:
        newID = species + '|' + geneID + '|' + geneType
        IDs[geneID] = newID
Input.close()

Pfam = open(args.pfam,'r')
Out1 = open(output1,'w')
Out2 = open(output2,'w')
hmms = {}
for line in Pfam:
    line = line.strip()
    gene = line.split()[0]
    hmm_name = line.split()[6]
    if gene in IDs.keys():
        if gene not in hmms.keys():
            hmms[gene]=[]
        hmms[gene].append(hmm_name)
for gene in hmms.keys():
    print(IDs[gene],'\t'.join(hmms[gene]),file=Out1, sep='\t')
    n_TIR,n_NBS,n_LRR,n_others = 0,0,0,0
    for domain in hmms[gene]:
        if re.match('TIR',domain):
            n_TIR += 1
        elif re.match('NB-ARC',domain):
            n_NBS += 1
        elif re.match('LRR',domain):
            n_LRR += 1
        else:
            n_others += 1
    print(IDs[gene],n_TIR,n_NBS,n_LRR,n_others,file=Out2, sep='\t')
Pfam.close()
Out1.close()
Out2.close()

