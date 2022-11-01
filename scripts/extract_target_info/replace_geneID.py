#!/user/bin/python3

#Description: replace gene IDs
#Author: Jiang Qian
#Data: 2021.11.15

import argparse
import os

parser = argparse.ArgumentParser(description='replace gene IDs')
parser.add_argument('-i', '--input', required=True, help="id type list")
parser.add_argument('-s', '--sequence', required=True, help="fas file of sequences")
parser.add_argument('-o', '--output', required=True, help="output dir")
args = parser.parse_args()

output = args.output + os.sep + "new.fas"

Input = open(args.input,'r')
IDs = {}
for line in Input:
    line = line.strip()
    geneID = line.split('\t')[0]
    type = line.split('\t')[1]
    newID = geneID + '|' + type
    IDs[geneID] = newID
Input.close()

F = open(args.sequence,'r')
Output = open(output,'w')
for line in F:
    line = line.strip()
    if line.startswith('>'):
        gene=line.replace('>','')
        print('>'+IDs[gene],file=Output)
    else:
        print(line,file=Output)
F.close()
Output.close()
