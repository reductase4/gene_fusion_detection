#!/user/bin/python3
#usage: python extract_gene_expression.py -i Camellia_sinensis_fusion_gene.txt -e Total_gene_fpkm.csv
import sys, os
import argparse

parser = argparse.ArgumentParser(description='extract gene expression data')
parser.add_argument('-i', '--blast', required=True, help="blast file of fused and parental genes derived from gene fusion detection pipeline")
parser.add_argument('-e', '--expressionData', required=True, help="gene expression data, like pfkm etc.")
parser.add_argument('-s', '--sepSymbol', default= ',', help="sep symbol used in gene expression data file, like ',', '\t'")
args = parser.parse_args()

blastFile = os.path.abspath(args.blast)
fpkmFile = os.path.abspath(args.expressionData)
sepSymbol = args.sepSymbol

blastResult = open(blastFile,'r')
blastDict = {}
for line in blastResult:
	line = line.strip()
	info = line.split('|')
	query = info[1]
	hit = info[3]
	hitType = info[4].split('\t')[0]
	if query not in blastDict.keys():
		blastDict[query] = {}
	elif hit in blastDict[query].keys():
		continue
	blastDict[query][hit]=hitType
blastResult.close()

totalFpkm = open(fpkmFile,'r')
fpkmDict={}
for line in totalFpkm:
	line = line.strip()
	if line.startswith('#'):
		colnames = line.lstrip('#')
		continue
	gene, fpkm = line.split(sepSymbol,1)
	fpkmDict[gene] = fpkm
totalFpkm.close()

Output = open(blastFile.rsplit('.',1)[0]+'_'+ fpkmFile.rsplit('_',1)[1],'w')
print('geneType',colnames,sep=',',file=Output)
for query in blastDict.keys():
	print('TNL',query,fpkmDict[query],sep=',',file=Output)
	for hit in blastDict[query].keys():
		print(blastDict[query][hit],hit,fpkmDict[hit],sep=',',file=Output)
Output.close()


