#!/user/bin/python3

#Author: Jiang Qian
#Data: 2021.11.19
#usage: python extract_target_info.py -d /data01/jiangqian/RGAugury/all_results/RGAugury_results_2/all_results_465 -s species.list -m RNL -t pep/hmm/domainSeq -o output

import argparse
import os,sys
import re

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
extract_hmm_script= bin_file + os.sep + "extract_target_hmms_v2.py"
extract_domain_script= bin_file + os.sep + "extract_target_domains_sep.py"

if os.path.exists(output):
    os.system('rm -r ' + output)
os.makedirs(output)

F = open(sepcies,'r')
Out = open(output+os.sep+'all_short_name.list','w')
for line in F:
    full_name = line.strip()
    names = full_name.split('_')
    #short_name = names[0][0]+names[1][0:2]
    short_name = full_name
    #produce full name echo short name file
    print(full_name,short_name,file=Out,sep='\t')
    #target gene list
    listFile = RGA_results + os.sep + full_name + os.sep + full_name + '.NBS.candidates.lst'
    #extract peps
    if infoType == "pep":
        fasFile = RGA_results + os.sep + full_name + os.sep + full_name + '.NBS.candidates.fas'
        if not os.path.exists(fasFile):
            print(fasFile + 'does not exist!')
            continue
        cmd = 'python ' + extract_pep_script + ' -i ' + listFile + ' -f ' + fasFile + ' -s ' + short_name + ' -m ' + matchType + ' -k ' + full_name + ' -o ' + output
        postfix = '.pep'
    #extract hmm names
    if infoType == "hmm":
        pfamFile = RGA_results + os.sep + full_name + os.sep + full_name + '.pfam.local.search.out'
        cmd = 'python ' + extract_hmm_script + ' -i ' + listFile + ' -p ' + pfamFile + ' -s ' + short_name + ' -m ' + matchType + ' -k ' + full_name + ' -o ' + output
        postfix = '.hmm'
    #extract domain name & seq
    if infoType == "domainSeq":
        domainFile = RGA_results + os.sep + full_name + os.sep + full_name + '.NBS.merged.domains.txt'
        cmd = 'python ' + extract_domain_script + ' -i ' + listFile + ' -d ' + domainFile + ' -s ' + short_name + ' -m ' + matchType + ' -n ' + domain_name + ' -k ' + full_name + ' -o ' + output
        postfix = '.domain'
    #run command
    os.system(cmd)
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
