#!/usr/bin/env bash

###
# USAGE: ./split_and_run_sparc.sh [BACKBONE_FASTA] [CONSENSUS_FASTA] [READS_FASTA] [OUTPUT_DIR] [ITERATIONS] ###

backbone_fasta=$1
consensus_fasta=$2
reads_fasta=$3
split_dir=$4
iterations=$5



#clean the directory first
rm ${split_dir}/backbone-*

./split_reads_by_backbone.py -b ${backbone_fasta} -o ${split_dir} -r ${reads_fasta} -c ${consensus_fasta} 

for file in $(ls ${split_dir}/*.reads.fasta); do
    chunk=`basename $file .reads.fasta`

    cmd=""
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

    echo $cmd
    eval $cmd


    #to save space
    cmd="rm ${split_dir}/${chunk}.mapped.m5"
    echo $cmd
    eval $cmd
    cmd="rm ${split_dir}/${chunk}.reads.fasta"
    echo $cmd
    eval $cmd

done

cmd="cat ${split_dir}/*.consensus.fasta > ${split_dir}/final_assembly.fasta"
echo=$cmd
eval $cmd