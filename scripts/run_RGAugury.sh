#!/bin/bash

#usage: nohup sh run_RGAugury_v3.sh /data01/jiangqian/RGAugury/all_data/ensemble /data01/jiangqian/RGAugury/all_results/ensemble & 

DATA_DIR=$1 #protein data dir
Result_DIR=$2 #results dir
taskn=5
key_pep="*.longest_transcript.pep.fa"

#max_bg_procs: used to check jobs & limit jobs number
#run x tasks all the time
function max_bg_procs {
    if [[ $# -eq 0 ]] ; then
            echo "Usage: max_bg_procs NUM_PROCS.  Will wait until the number of background (&)"
            echo "           bash processes (as determined by 'jobs -pr') falls below NUM_PROCS"
            return
    fi
    local max_number=$((0 + ${1:-0}))
    while true; do
            local current_number=$(jobs -pr | wc -l)
            if [[ $current_number -lt $max_number ]]; then
                    break
            fi
            sleep 10
    done
}

function runRGAugury {
    echo "$1 start time $(date +"%Y-%m-%d %T")"
    mkdir -p ${Result_DIR}/${1}     ##mkdir result dir
    cd ${Result_DIR}/${1}    ##cd to the result dir
    path=`pwd`
    echo "RGAugury_modified_v2.2.pl -p $2 -pfx $1 -c 4 -d SUPERFAMILY,SMART,panther,gene3d,Pfam"
    RGAugury_modified_v2.2.pl -p $2 -pfx $1 -c 4 -d SUPERFAMILY,SMART,panther,gene3d,Pfam
    echo "$1 finish time $(date +"%Y-%m-%d %T")"
}

for file in $DATA_DIR/*
do
    species=`basename $file`    ##get the file name
    if [ ! -d "$Result_DIR/$species" ]; then
        max_bg_procs $taskn
        ((count++))
        runRGAugury $species $file/$key_pep&
    fi
done

wait
echo "All $count files done"



