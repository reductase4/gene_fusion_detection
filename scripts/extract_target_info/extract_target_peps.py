#!/user/bin/python3

#Description: extract target peps with new IDs for target sepecies
#Author: Jiang Qian
#Data: 2021.11.18

import argparse
import os
import re

parser = argparse.ArgumentParser(description='extract target peps with new IDs for target sepecies')
parser.add_argument('-i', '--input', required=True, help="id type list")
parser.add_argument('-f', '--fasta', required=True, help="fas file of sequences")
parser.add_argument('-s', '--short', required=True, help="sepcies short name")
parser.add_argument('-m', '--match', required=True, help="match word")
parser.add_argument('-o', '--output', required=True, help="output dir")
parser.add_argument('-k', '--key', required=True, help="output keyword(e.g., species full name)")
args = parser.parse_args()

species = args.short
#matchWord = args.match
matchWord = (args.match).split(',')
output = args.output + os.sep + args.key + "_target_pep." + "_".join(matchWord)

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

F = open(args.fasta,'r')
Output = open(output,'w')
seq = {}
gene = ''
for line in F:
    line = line.strip()
    if line.startswith('>'):
        if gene in IDs.keys() and seq[gene] != '':
            print('>'+IDs[gene],file=Output)
            print(seq[gene],file=Output)
        gene=line.replace('>','')
        seq[gene]=''
    else:
        seq[gene]+=line
if gene in IDs.keys():
    print('>'+IDs[gene],file=Output)
    print(seq[gene],file=Output)

F.close()
Output.close()
