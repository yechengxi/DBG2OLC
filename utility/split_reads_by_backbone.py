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
            help="Fasta file where header corresponds to backbone id and entries correspond to read ids.")
    parser.add_option("-r", "--reads", dest="reads_filename", \
            help="Fasta file of the reads.")
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
    tuple = pf.getRecord()
    while tuple is not None:
        print tuple[0]

        split_backbone = open(options.output_dir + '/' + options.prefix + '-' + str(id_counter) + '.fasta', 'w')
        split_backbone.write('>' + tuple[0] + '\n' + tuple[1])
        split_backbone.close()

        backbone_to_id[tuple[0]] = options.prefix + '-' + str(id_counter)

        id_counter += 1
        tuple = pf.getRecord()

    return backbone_to_id


def build_reads_to_backbone_dict(options):
    """
    Return dictionary of reads to the corresponding backbone.
    """

    """
    >Backbone_1
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/167/0_4524
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/1566/0_9679
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/1888/9474_18333
    """
    reads_to_backbone = defaultdict(list)

    curr_backbone = None
    for line in open(options.consensus_filename, 'r'):
        line = line.strip()
        if line.startswith('>'):
            curr_backbone = line[1:]

        else:
            if curr_backbone:
                print line + '\t' + curr_backbone
                reads_to_backbone[line].append(curr_backbone)

    return reads_to_backbone


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

    # Find the reads that correspond to a given backbone.
    reads_to_backbone = build_reads_to_backbone_dict(options)

    # Split the reads based on their corresponding backbone file.
    file_pointers_dict = {}

    pf = ParseFasta(options.reads_filename)
    tuple = pf.getRecord()
    id = None
    while tuple is not None:

        if len(reads_to_backbone[tuple[0]]) > 0:
            for backbone in reads_to_backbone[tuple[0]]:
                id = backbone_to_id[backbone]

                #print tuple[0] + '\t-->\t' + options.output_dir + '/' + str(id) + '.reads.fasta'

                if id in file_pointers_dict:
                    file_pointers_dict[id].write('>' + tuple[0] + '\n' + tuple[1] + '\n')

                else:
                    new_fp = open(options.output_dir + '/' + str(id) + '.reads.fasta', 'w')
                    file_pointers_dict[id] = new_fp
                    file_pointers_dict[id].write('>' + tuple[0] + '\n' + tuple[1] + '\n')
        #else:
        #    print 'MISSING\t' + tuple[0]

        tuple = pf.getRecord()

if __name__ == '__main__':
    main()
