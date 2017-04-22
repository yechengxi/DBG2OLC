#!/usr/bin/env bash

if [ $# -eq 0 ]; then
 echo "
 ###
 # USAGE: split_and_run_pbdagcon.sh [BACKBONE_FASTA] [CONSENSUS_FASTA] [READS_FASTA] [OUTPUT_DIR][split_reads_by_backbone_version 1,2,3 - 3 by default]
 ###
 ## split_reads_by_backbone_version
 ## 1 = split_reads_by_backbone.py -- can cause "Too many files open" error. Good for small genomes and/or assemblies with few contigs.
 ## 2 = split_reads_by_backbone_readdict.py -- much faster, but requires enough memory to load all reads into memory. Will not cause "too many open files". Operation is a little different so order of reads.fastas can be different.
 ## 3 = split_reads_by_backbone_openclose.py -- a little slower than #1, but will not cause "too many open files". Will produce fastas in same order as #1. Make sure outdir is empty as it will append to existing files if they have the same name.
 ###
 "
 exit
fi

backbone_fasta=$1
consensus_fasta=$2
reads_fasta=$3
split_dir=$4

if [ $# -eq 6 ]; 
 then splitversion=$6
else
 splitversion=3
fi


# if dir not there, make it; else clean the dir out
if [ ! -d $split_dir ]; then mkdir $splitdir ; else  rm ${split_dir}/backbone-* ; fi


if [ $splitversion -eq 1 ]; then
 split_reads_by_backbone.py -b ${backbone_fasta} -o ${split_dir} -r ${reads_fasta} -c ${consensus_fasta} 
elif [ $splitversion -eq 2 ]; then
 split_reads_by_backbone_readdict.py -b ${backbone_fasta} -o ${split_dir} -r ${reads_fasta} -c ${consensus_fasta} 
elif [ $splitversion -eq 3 ]; then
 split_reads_by_backbone_openclose.py -b ${backbone_fasta} -o ${split_dir} -r ${reads_fasta} -c ${consensus_fasta} 
else
 echo "splitversion (argument #6) needs to be 1, 2, or 3. If left blank, default is 1."
 exit
fi


for file in $(ls ${split_dir}/*.reads.fasta); do
    chunk=`basename $file .reads.fasta`
    #echo $chunk
    cmd="blasr -nproc 64 ${split_dir}/${chunk}.reads.fasta ${split_dir}/${chunk}.fasta -bestn 1 -m 5 -minMatch 19 -out ${split_dir}/${chunk}.mapped.m5"
    echo $cmd
    eval $cmd

    cmd="pbdagcon ${split_dir}/${chunk}.mapped.m5 -j 1 -c 1 -m 200 > ${split_dir}/${chunk}.consensus.fasta"
    echo $cmd
    eval $cmd

done

cmd="cat ${split_dir}/*.consensus.fasta > ${split_dir}/final_assembly.fasta"
echo=$cmd
eval $cmd
