#!/bin/sh

#  split_and_create_cns_batches.sh
#
#
#  Created by Chengxi Ye on 9/16/15.
#


#!/usr/bin/env bash
###
# USAGE: ./split_and_run_sparc.sh [BACKBONE_FASTA] [CONSENSUS_FASTA] [READS_FASTA] [OUTPUT_DIR] [ITERATIONS] [BATCHES]
###

backbone_fasta=$1
consensus_fasta=$2
reads_fasta=$3
split_dir=$4
iterations=$5
batches=$6

echo "batches: "${batches}


prefix="cns_batch"


#clean the directory first
rm ./${prefix}*
rm ${split_dir}/backbone-*

./split_reads_by_backbone.py -b ${backbone_fasta} -o ${split_dir} -r ${reads_fasta} -c ${consensus_fasta}


total_files=$(ls ${split_dir}/*.reads.fasta | wc -l)
echo ${total_files}

files_per_batch=$(( ${total_files} / ${batches} ))

echo ${files_per_batch}

count=0
batch=0
cmd="source ~/.bashrc; "

for file in $(ls ${split_dir}/*.reads.fasta); do
chunk=`basename $file .reads.fasta`

cmd+="if [ -f ${split_dir}/${chunk}.reads.fasta ]; then "

for iter in `seq 1 ${iterations}`; do

#echo $iter
cmd+="blasr -nproc 64 ${split_dir}/${chunk}.reads.fasta ${split_dir}/${chunk}.fasta -bestn 1 -m 5 -minMatch 19 -out ${split_dir}/${chunk}.mapped.m5; "

cmd+="./Sparc m ${split_dir}/${chunk}.mapped.m5 b ${split_dir}/${chunk}.fasta k 1 c 2 g 1 HQ_Prefix Contig boost 5 t 0.2 o ${split_dir}/${chunk}; "

if [ ${iter} -lt ${iterations} ]
then
#rename
cmd+="mv ${split_dir}/${chunk}.consensus.fasta ${split_dir}/${chunk}.fasta;"
fi

done

#to save space
cmd+="rm ${split_dir}/${chunk}.mapped.m5; "
cmd+="rm ${split_dir}/${chunk}.reads.fasta; "
cmd+=" fi;"

count=$(($count + 1))

if [ $count -gt ${files_per_batch} ]
then
batch=$(($batch + 1))
echo $cmd >"${prefix}${batch}.txt"
#eval $cmd
cmd="source ~/.bashrc; "
count=0
fi

done


# Print left over sequences.
if [ $count -gt 0 ]
then
batch=$(($batch + 1))
echo $cmd >"${prefix}${batch}.txt"
fi


