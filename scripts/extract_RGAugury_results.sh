#!/bin/bash
#Creator: Jiang Qian
#Date: 2021.11.02
#Description: extract stat info from results of RGAugury
#usage: sh extract_RGAugury_results_v2.sh /data01/jiangqian/RGAugury/all_results/ensemble/ ensemble

#Version: v2.0
#Version modified: for more domain class.
#Version: v2.1 2021.12.13
#Version modified: add counts of LRR and CC domain.

DATA_DIR=$1 #results dir
Key=$2

dOut=`dirname $DATA_DIR`
outfile=$dOut'/results_'$Key'.csv'

echo `date`

echo -e "Name\tfileSize\tseqN\tNBS\tTX\tRX\tLRR\tCC\tCL\tCN\tTN\tRN\tNL\tTC\tRC\tTCN\tRCN\tCNL\tTNL\tRNL\tTCNL\tRCNL\tRTCNL\tOthers\tRLP\tRLK\tTM-CC" > $outfile

for file in $DATA_DIR/*
do
    ((num++))
    echo "################ $num #######################"
    species=`basename $file`        ##get the file name
	echo "Processing No.$num: $species"
	fileSize=`ls -lh $file/*formated.protein.input.fas|awk '{print $5}'`        ##get protein file size
    seqN=`grep -c '^>' $file/*formated.protein.input.fas`        ##get sequence number

    statFile=$file/*.summaries.txt
    if test -e $statFile
    then
        statN=`sed -n "2p" $statFile|sed 's/[ \t]*$//g'`        ##get RGAs stat numbers
    fi
    echo -e "$species\t$fileSize\t$seqN\t$statN" >> $outfile
    runTime=''
    statN=''
done

echo "Done"
echo `date`

