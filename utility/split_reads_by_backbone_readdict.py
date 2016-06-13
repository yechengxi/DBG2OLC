#!/usr/bin/env python

from SeqIO import ParseFasta
from collections import defaultdict
from optparse import OptionParser
from sets import Set
import os, sys
import re


def setup_options():
    parser = OptionParser()
    parser.add_option("-b", "--backbone", dest="backbone_filename", \
            help="Backbone fasta file containing contigs with gaps.", metavar="FILE")
    parser.add_option("-c", "--consensus", dest="consensus_filename", \
            help="Fasta file where header corresponds to backbone id and entries correspond to read ids. Typically 'DBG2OLC_Consensus_info.txt' in pipeline.")
    parser.add_option("-r", "--reads", dest="reads_filename", \
            help="Fasta file of the reads. Typically 'ctg_pb.fasta' in pipeline. ")
    parser.add_option("-o", "--output-dir", dest="output_dir", \
            help="Output directory of split fasta files.")
    parser.add_option("-p", "--prefix", dest="prefix", default="backbone", \
            help="Default prefix for split fasta.")


    (options, args) = parser.parse_args()

    return (options,args)


def split_backbone(options):
    """ 
    Split backbone fasta file into chunks. 
    Returns dictionary of backbone -> id.
    """

    backbone_to_id = {}
    id_counter = 0

    # Write all backbone files to their own fasta file.
    pf = ParseFasta(options.backbone_filename)
    tuple = pf.getRecord() # [seqname, sequence]
    while tuple is not None:
        print tuple[0]

        split_backbone = open(options.output_dir + '/' + options.prefix + '-' + str(id_counter) + '.fasta', 'w')
        split_backbone.write('>' + tuple[0] + '\n' + tuple[1])
        split_backbone.close()

        # backbone_to_id[seqname] = fileprefix-number
        backbone_to_id[tuple[0]] = options.prefix + '-' + str(id_counter)

        id_counter += 1
        tuple = pf.getRecord()

    return backbone_to_id


def build_reads_to_backbone_dict(options):
    """
    Return dictionary of reads to the corresponding backbone.
    """

    """
    Usually in DBG2OLC_Consensus_info.txt
    >Backbone_1
    Contig_59850
    Contig_61226
    Contig_87853
    Contig_42901
    Contig_308247
    Contig_180195
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/167/0_4524
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/1566/0_9679
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/1888/9474_18333
    """
    reads_to_backbone = defaultdict(list)

    curr_backbone = None
    with open(options.consensus_filename, 'r') as cf:
        for line in cf:
            line = line.strip()
            if line.startswith('>'):
                curr_backbone = line[1:]

            else:
                if curr_backbone:
                    print line + '\t' + curr_backbone
                    ## reads_to_backbone[read_or_contig_name].append(backbone_name)
                    reads_to_backbone[line].append(curr_backbone)

    return reads_to_backbone

def build_backbone_to_reads_dict(options):
    """
    Return dictionary of backbones to the corresponding reads.
    """

    """
    Usually in DBG2OLC_Consensus_info.txt
    >Backbone_1
    Contig_59850
    Contig_61226
    Contig_87853
    Contig_42901
    Contig_308247
    Contig_180195
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/167/0_4524
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/1566/0_9679
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/1888/9474_18333
    """
    backbone_to_reads = defaultdict(list)

    curr_backbone = None
    with open(options.consensus_filename, 'r') as cf:
        for line in cf:
            line = line.strip()
            if line.startswith('>'):
                curr_backbone = line[1:]

            else:
                if curr_backbone:
##                    print line + '\t' + curr_backbone
                    ## backbone_to_reads[backbone_name].append(read_or_contig_name)
                    backbone_to_reads[curr_backbone].append(line)

    return backbone_to_reads


def build_readseqs_dict(options):
    readseqs = {}
    pf = ParseFasta(options.reads_filename) ## typically 'ctg_pb.fasta' in pipeline
    tuple = pf.getRecord()
    while tuple is not None:
        readseqs[tuple[0]] = tuple[1]
        tuple = pf.getRecord()
    return readseqs
        

def ensure_dir(f):
    print f
    d = os.path.dirname(f)
    print d
    if not os.path.exists(d):
        os.makedirs(d)


def main():
    (options, args) = setup_options()

    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)

    # Split backbone fasta file into chunks.
    backbone_to_id = split_backbone(options)
    # use: backbone_to_id[seqname] = fileprefix-number

    # Find the reads that correspond to a given backbone.
##    reads_to_backbone = build_reads_to_backbone_dict(options)
    ## use: reads_to_backbone[read_or_contig_name].append(backbone_name)
    
    
    ## New approach: faster than original, requires more memory
    ## create readname:seq dict -> readseqs
    ## have a back_bone_to_reads dict instead
    ## do:
    ## for backbone in backbone_to_reads:
    ##    id = backbone_to_id[backbone]
    ##    new_fp = open(options.output_dir + '/' + str(id) + '.reads.fasta', 'w')
    ##    for read in that backbone_to_reads[backbone]:
    ##        new_fp.write('>' + reads + '\n' + readseqs[read] + '\n')

    backbone_to_reads = build_backbone_to_reads_dict(options)
    ## use: backbone_to_reads[backbone_name].append(read_or_contig_name)

    readseqs = build_readseqs_dict(options)

    for backbone in backbone_to_reads:
        id = backbone_to_id[backbone]
        new_fp = open(options.output_dir + '/' + str(id) + '.reads.fasta', 'w')
        for read in backbone_to_reads[backbone]:
            new_fp.write('>' + read + '\n' + readseqs[read] + '\n')
        new_fp.close()

if __name__ == '__main__':
    main()
