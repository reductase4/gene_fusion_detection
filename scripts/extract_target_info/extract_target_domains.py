#!/user/bin/python3

#Author: Jiang Qian
#Data: 2021.11.23
#extract all domains in domain pos file

import argparse
import os
import re

parser = argparse.ArgumentParser(description='extract domains of target genes for target species')
parser.add_argument('-i', '--input', required=True, help="target id/type list")
parser.add_argument('-d', '--domain', required=True, help="domain pos file")
parser.add_argument('-s', '--short', required=True, help="sepcies short name")
parser.add_argument('-m', '--match', required=True, help="match word")
parser.add_argument('-o', '--output', required=True, help="output dir")
parser.add_argument('-k', '--key', required=True, help="output keyword (e.g., species full name)")
args = parser.parse_args()

species = args.short
matchType = args.match
key = args.key
domainPos = args.domain
seqfile = os.path.dirname(os.path.abspath(domainPos)) + os.sep + key + '.formated.protein.input.fas'
output1 = args.output + os.sep + args.key + "_domain_name." + matchType
output2 = args.output + os.sep + args.key + "_domain_seq." + matchType

Input = open(args.input,'r')
IDs = {}
for line in Input:
    line = line.strip()
    geneID = line.split('\t')[0]
    geneType = line.split('\t')[1]
    if re.match(matchType,geneType):
        newID = species + '|' + geneID + '|' + geneType
        IDs[geneID] = newID
Input.close()

PosFile = open(domainPos,'r')
posdict = {}
for line in PosFile:
    line = line.strip()
    if line.startswith('id\tLen'):
        continue
    gene = line.split()[0]
    if gene not in IDs.keys():
        continue
    domains = line.split()[2:]
    posdict[gene] = []
    for domain in domains:
        if domain == '.':
            continue
        string,pos = domain.split('|') #domain_NBS|133-370 or IPR_domain_LRR|633-890
        name = string.split('domain_')[1]
        #start,end = pos.split('-')
        start,end = [int(i) for i in pos.split('-')]
        posdict[gene].append({'name':name,'start':start,'end':end})
'''
for gene in posdict.keys():
    print(posdict[gene])
    for i in range(len(posdict[gene])):
        print(posdict[gene][i])
        exit()
'''
#sort posdict
posdict_sort = {}
for gene in posdict.keys():
    posdict_sort[gene] = []
    #sort
    sortList = sorted(posdict[gene],key=lambda i: i['start'])
    #remove duplications
    noDupList = []
    for item in sortList:
        if not item in noDupList:
            noDupList.append(item)
    #merge overlaps
    for i in range(len(noDupList)):
        name2 = noDupList[i]['name']
        start2 = noDupList[i]['start']
        end2 =noDupList[i]['end']
        #print(name2,start2,end2)
        if i == 0:
            name1,start1,end1 = name2,start2,end2
            continue
        if end1 < start2:
            posdict_sort[gene].append({'name':name1,'start':start1,'end':end1})
            name1,start1,end1 = name2,start2,end2
        elif end1 >= start2: #overlap
            if end1 > end2: #record 1 is larger
                start = start1
                end = end1
                name = name1 + '_' + name2
            else: #record 1 + record 2
                start = start1
                end = end2
                name = name1 + '_' + name2
            name1,start1,end1 = name,start,end
        #print(name,start,end)
        #posdict_sort[gene].append({'name':name,'start':start,'end':end})
    posdict_sort[gene].append({'name':name1,'start':start1,'end':end1})

#read sequence file
SeqFile = open(seqfile,'r')
Out1 = open(output1,'w')
Out2 = open(output2,'w')
seq = {}
gene = ''
for line in SeqFile:
    line = line.strip()
    if line.startswith('>'):
        gene = line.replace('>','')
        seq[gene] = ''
    else:
        seq[gene] += line

for gene in posdict_sort.keys():
    domain_name = []
    domain_seq = ''
    for i in range(len(posdict_sort[gene])):
        info = posdict_sort[gene][i]
        domain_name.append(info['name'])
        domain_seq += seq[gene][info['start']-1:info['end']]
    print(IDs[gene],'\t'.join(domain_name),file=Out1,sep='\t')
    print('>'+IDs[gene],file=Out2)
    print(domain_seq,file=Out2)

SeqFile.close()
Out1.close()
Out2.close()

