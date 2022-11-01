#!/user/bin/python3
#Author: Jiang Qian
#Data: 2022.10.31

import argparse
import os,sys
import re

parser = argparse.ArgumentParser(description='gene fusion detection')
parser.add_argument('-d', '--dir', required=True, help="genome files")
parser.add_argument('-s', '--species', required=True, help="target species")
parser.add_argument('-o', '--order', required=True, help="species order")
parser.add_argument('-e', '--expression', required=False, help="expression dir")
args = parser.parse_args()

bin_file = os.path.dirname(sys.argv[0])
Genome_File = os.path.abspath(args.dir)
species = os.path.abspath(args.species)
species_order = os.path.abspath(args.order)
step4 = False
if args.expression:
    exp_dir = os.path.abspath(args.expression)
    step4 = True

#scripts
run_RGAugury_script = bin_file + os.sep + "run_RGAugury.sh"
extract_RGAugury_results_script = bin_file + os.sep + "extract_RGAugury_results.sh"
extract_target_script = bin_file + os.sep + "extract_target_info/extract_target_info_v2.0.py"
filter_script = bin_file + os.sep + "blastp_filter_two_steps.py"
extract_gene_expression_script = bin_file + os.sep + "extract_gene_expression.py"

#dirs
pre_data_file = os.path.dirname(bin_file) + os.sep + 'pre_data'
Results_File = os.getcwd() + os.sep + 'results'

if os.path.exists(Results_File):
    print(Results_File + ' exists!')
else:
    os.makedirs(Results_File)

work_steps = Results_File + os.sep + 'work_steps'
RGAugury_results = Results_File + os.sep + 'RGAugury_results'
pep_results = Results_File + os .sep + 'pep'
fusion_results = Results_File + os .sep + 'fusion_detection_results'
expression_results = Results_File + os .sep + 'expression_results'

#work_steps
if os.path.exists(work_steps):
    os.system('rm -r ' + work_steps)
os.makedirs(work_steps)

#step1 RGAugury
print('#####################')
print('Step1.RGAugury...')
SH = work_steps + os.sep + 'step1.RGAugury.sh'
sh = open(SH,'w')
#run RGAugury
cmd = 'sh ' + run_RGAugury_script + ' ' + Genome_File + ' ' + RGAugury_results
print(cmd,file=sh)
#extract results
cmd = 'sh ' + extract_RGAugury_results_script + ' ' + RGAugury_results + ' RGAugury'
print(cmd,file=sh)
sh.close()
sys.stdout.flush()
#run command
os.system('sh ' + SH)

#step2 extract peps
#extract_target_info.py for all; extract_target_info_v2.0.py for LRRfilter version
#extract TNL/NL/TX/TN peps
print('#####################')
print('Step2 extract peps...')
SH = work_steps + os.sep + 'step2.extract_pep.sh'
sh = open(SH,'w')
target_type = ['TNL','NL','TX']
for t in target_type:
    cmd = 'python ' + extract_target_script + ' -d ' + RGAugury_results + ' -s ' + species + ' -m ' + t +' -t pep -o ' + pep_results + os.sep + t
    print(cmd,file=sh)
sh.close()
sys.stdout.flush()
#run command
os.system('sh ' + SH)

#step3 fused gene detection
print('#####################')
print('Step3 fused gene detection...')
SH = work_steps + os.sep + 'step3.fused_gene_detection.sh'
sh = open(SH,'w')
pre_TNL = pre_data_file + os.sep + 'all_TNL_45_species_pep.txt'
pre_NL = pre_data_file + os.sep + 'all_NL_45_species_pep.txt'
pre_TX = pre_data_file + os.sep + 'all_TX_45_species_pep.tx'
new_TNL = pep_results + os.sep + 'TNL' +os.sep + 'all_TNL.pep'
new_NL = pep_results + os.sep + 'NL' +os.sep + 'all_NL.pep'
new_TX = pep_results + os.sep + 'TX' +os.sep + 'all_TX.pep'

fused_candidates = fusion_results + os.sep + 'fused_gene_candidates.pep'
parental_candidates = fusion_results + os.sep + 'parental_gene_candidates.pep'

fusedDB = fusion_results + os.sep + 'fusedDB.pep'
parentalDB = fusion_results + os.sep + 'parentalDB.pep'

fusion_parents_blast = fusion_results + os.sep + 'fusion_parents_fmt7.blast'
fusion_orthology_blast = fusion_results + os.sep + 'fusion_orthology_fmt7.blast'

keywords = 'fused_gene_candidates'
detection_results = fusion_results + os.sep + keywords + '.filter2'
detection_results_list = detection_results + '.list'
common_results = fusion_results + os.sep + keywords + '.common'
common_results_list = common_results + '.list'

if os.path.exists(fusion_results):
    print(fusion_results + ' exists!')
else:
    os.makedirs(fusion_results)

if os.path.exists(fusion_parents_blast) and os.path.exists(fusion_orthology_blast):
    print(fusion_parents_blast + ' exists!')
    print(fusion_orthology_blast + ' exists!')
    print('Skip blastp steps: blast result files exist.')
else:
    #fusedDB
    cmd = 'cat ' + new_TNL + ' ' + pre_TNL + ' > ' + fused_candidates
    print('#fused_candidates',file=sh)
    print(cmd,file=sh)

    #parentalDB
    cmd = 'cat ' + new_NL + ' ' + pre_NL + ' ' + new_TX + ' ' + pre_TX + ' > ' + parental_candidates
    print('#parental_candidates',file=sh)
    print(cmd,file=sh)

    #makedb
    print('#makedb',file=sh)
    cmd = 'makeblastdb -in ' + fused_candidates + ' -dbtype prot -out ' + fusedDB
    print(cmd,file=sh)

    cmd = 'makeblastdb -in ' + parental_candidates + ' -dbtype prot -out ' + parentalDB
    print(cmd,file=sh)

    #run blastp
    print('#run blastp',file=sh)
    string = '-outfmt \"7 qseqid sseqid pident qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore\" -evalue 1e-5 -num_threads 10'
    cmd = 'blastp -query ' + fused_candidates + ' -out ' + fusion_orthology_blast + ' -db ' + fusedDB + ' ' + string
    print(cmd,file=sh)

    cmd = 'blastp -query ' + fused_candidates + ' -out ' + fusion_parents_blast + ' -db ' + parentalDB + ' ' + string
    print(cmd,file=sh)

#run filter
print('#run filter',file=sh)
cmd = 'python ' + filter_script + ' -i1 ' + fusion_parents_blast + ' -i2 ' + fusion_orthology_blast + ' -s ' + species_order + ' -o ' + fusion_results + os.sep + keywords
print(cmd,file=sh)

cmd = 'cat ' + detection_results + '|awk -F \'[\\t]\' \'{print $1}\'|sort|uniq > ' + detection_results_list
print(cmd,file=sh)

cmd = 'cat ' + common_results + '|awk -F \'[\\t]\' \'{print $1}\'|sort|uniq > ' + common_results_list
print(cmd,file=sh)
sh.close()
sys.stdout.flush()
#run command
os.system('sh ' + SH)

#step4 gene expression data extraction
print('#####################')
if step4:
    print('Step4 extract gene expression data of fused & parental genes...')
    SH = work_steps + os.sep + 'step4.extract_gene_expression.sh'
    sh = open(SH,'w')
    
    if os.path.exists(expression_results):
        os.system('rm -r ' + expression_results)
    os.makedirs(expression_results)
    
    f = open(species,'r')
    target_species = []
    for line in f:
        line = line.strip()
        target_species.append(line)
    f.close()
    
    for sp in target_species:
        exp_data = exp_dir + os.sep + sp + '_total_gene_expression.csv'
        target_results = expression_results + os.sep + sp + '_' + keywords + '.txt'
        if os.path.exists(exp_data): 
            cmd = 'grep ' + sp + ' ' + detection_results + ' > ' + target_results
            print(cmd,file=sh)
            cmd = 'python ' + extract_gene_expression_script + ' -i ' + target_results + ' -e ' + exp_data
            print(cmd,file=sh)
    sh.close()
    sys.stdout.flush()
    #run command
    os.system('sh ' + SH)
    print('Done')
else:
    print('Skip step4: no expression data input.')

