#!/user/bin/python3

#Author: Jiang Qian
#Data: 2021.11.19
#usage: python extract_target_info.py -d /data01/jiangqian/RGAugury/all_results/RGAugury_results_2/all_results_465 -s species.list -m RNL -t pep/hmm/domainSeq -o output
#modified: 2022.8.30
#v2.0: LRR filter by pfam results

import argparse
import os,sys
import re

def get_LRR(pfam_results):
    LRR_dict = {}
    Infile = open(pfam_results,'r')
    for line in Infile:
        line = line.strip()
        gene = line.split()[0]
        domain = line.split()[6]
        evalue = line.split()[12]
        if "LRR" in domain and float(evalue) < 0.1: #pfam_scam 0.001, NCBI CD-search 0.01
            if gene not in LRR_dict.keys():
                LRR_dict[gene]=[]
            LRR_dict[gene].append(domain)
    Infile.close()
    return(LRR_dict)

def LRR_filter(gene_list,LRR_dict):
    Infile = open(gene_list,'r')
    Outfile = open(gene_list+'.new','w')
    for line in Infile:
        line = line.strip()
        gene, gene_type = line.split()
        if "L" in gene_type and gene not in LRR_dict.keys():
            print(line)
            gene_type = gene_type.replace("L","")
            if gene_type == "N":
                gene_type = "NBS"
            #print(gene+'\t'+ gene_type)
            print(gene,gene_type,sep='\t',file=Outfile)
            continue
        print(line,file=Outfile)
    Infile.close()
    Outfile.close()
    return(Outfile)

parser = argparse.ArgumentParser(description='extract target info from RGAugury results')
parser.add_argument('-d', '--dir', required=True, help="RGAugury results dir")
parser.add_argument('-s', '--sepcies', required=True, help="sepcies list")
parser.add_argument('-m', '--match', required=True, help="match word(e.g., RNL/CNL/TNL)")
parser.add_argument('-n', '--name', default='all', help="domain name(e.g., all/NBS/NBS,LRR)")
parser.add_argument('-t', '--info', required=True, help="extract info type(e.g.,pep/hmm/domainSeq)")
parser.add_argument('-o', '--output', required=True, help="output dir")
args = parser.parse_args()

RGA_results = os.path.abspath(args.dir)
sepcies = args.sepcies
matchType = args.match
infoType = args.info
domain_name = args.name
output = os.path.abspath(args.output)
bin_file = os.path.dirname(sys.argv[0])
extract_pep_script= bin_file + os.sep + "extract_target_peps.py"
extract_hmm_script= bin_file + os.sep +"extract_target_hmms_v2.py"
extract_domain_script= bin_file + os.sep + "extract_target_domains_sep.py"

if os.path.exists(output):
    os.system('rm -r ' + output)
os.makedirs(output)

F = open(sepcies,'r')
for line in F:
    full_name = line.strip()
    names = full_name.split('_')
    short_name = full_name
    #target gene list
    listFile = RGA_results + os.sep + full_name + os.sep + full_name + '.NBS.candidates.lst'
    #filter LRR domain by pfamfile
    pfamFile = RGA_results + os.sep + full_name + os.sep + full_name + '.pfam.local.search.out'
    newlist = listFile + '.new'
    if not os.path.isfile(newlist):
        genes_with_LRR = get_LRR(pfamFile)
        LRR_filter(listFile, genes_with_LRR)
    #extract peps
    if infoType == "pep":
        fasFile = RGA_results + os.sep + full_name + os.sep + full_name + '.NBS.candidates.fas'
        if not os.path.exists(fasFile):
            print(fasFile + 'does not exist!')
            continue
        cmd = 'python ' + extract_pep_script + ' -i ' + newlist + ' -f ' + fasFile + ' -s ' + short_name + ' -m ' + matchType + ' -k ' + full_name + ' -o ' + output
        postfix = '.pep'
    #extract hmm names
    if infoType == "hmm":
        pfamFilter = RGA_results + os.sep + full_name + os.sep + full_name + '.pfam.local.search.out.filter'
        cmd = 'python ' + extract_hmm_script + ' -i ' + newlist + ' -p ' + pfamFilter + ' -s ' + short_name + ' -m ' + matchType + ' -k ' + full_name + ' -o ' + output
        postfix = '.hmm'
    #extract domain name & seq
    if infoType == "domainSeq":
        domainFile = RGA_results + os.sep + full_name + os.sep + full_name + '.NBS.merged.domains.txt'
        cmd = 'python ' + extract_domain_script + ' -i ' + newlist + ' -d ' + domainFile + ' -s ' + short_name + ' -m ' + matchType + ' -n ' + domain_name + ' -k ' + full_name + ' -o ' + output
        postfix = '.domain'
    #run command
    os.system(cmd)
#Out.close
#merge files
if infoType == "domainSeq":
    cmd = 'cat ' + output + os.sep + '*' + '_domain_name.'+ matchType + ' > ' + output + os.sep + 'all_domain_name_' + matchType + postfix
    print(cmd)
    os.system(cmd)
    cmd = 'cat ' + output + os.sep + '*' + '_domain_seq.'+ matchType + ' > ' + output + os.sep + 'all_domain_seq_' + matchType + postfix
    print(cmd)
    os.system(cmd)
elif infoType == "hmm":
    cmd = 'cat ' + output + os.sep + '*' + '_pfam_out.'+ matchType + ' > ' + output + os.sep + 'all_pfam_out_' + matchType + postfix
    print(cmd)
    os.system(cmd)
    cmd = 'cat ' + output + os.sep + '*' + '_pfam_stat.'+ matchType + ' > ' + output + os.sep + 'all_pfam_stat_' + matchType + postfix
    print(cmd)
    os.system(cmd)
else:
    cmd = 'cat ' + output + os.sep + '*' + matchType + ' > ' + output + os.sep + 'all_' + matchType + postfix
    print(cmd)
    os.system(cmd)
F.close()
