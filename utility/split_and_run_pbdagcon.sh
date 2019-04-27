#!/usr/bin/env bash

###
# USAGE: ./split_and_run_pbdagcon.sh [BACKBONE_FASTA] [CONSENSUS_FASTA] [READS_FASTA] [OUTPUT_DIR]
###

backbone_fasta=$1
consensus_fasta=$2
reads_fasta=$3
split_dir=$4

#clear the directory first
rm -rf ${split_dir}/backbone-*

./split_reads_by_backbone.py -b ${backbone_fasta} -o ${split_dir} -r ${reads_fasta} -c ${consensus_fasta} 

for file in $(find ${split_dir} -name "*.reads.fasta"); do
    chunk=`basename $file .reads.fasta`
    #echo $chunk
    cmd="blasr -nproc 64 ${split_dir}/${chunk}.reads.fasta ${split_dir}/${chunk}.fasta -bestn 1 -m 5 -minMatch 19 -out ${split_dir}/${chunk}.mapped.m5"
    echo $cmd
    eval $cmd

    cmd="pbdagcon ${split_dir}/${chunk}.mapped.m5 -j 1 -c 1 -m 200 > ${split_dir}/${chunk}.consensus.fasta"
    echo $cmd
    eval $cmd

done

echo="merging consensus reads to ${split_dir}/final_assembly.fasta"
for f in ${split_dir}/*.consensus.fasta; do cat "$f" >> ${split_dir}/final_assembly.fasta; done
