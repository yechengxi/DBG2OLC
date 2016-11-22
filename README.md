#DBG2OLC:Efficient Assembly of Large Genomes Using Long Erroneous Reads of the Third Generation Sequencing Technologies
=======

Our work is published in Scientific Reports: 

Ye, C. et al. DBG2OLC: Efficient Assembly of Large Genomes Using Long Erroneous Reads of the Third Generation Sequencing Technologies. Sci. Rep. 6, 31900; doi: 10.1038/srep31900 (2016).

http://www.nature.com/articles/srep31900

The manual can be downloaded from:

https://github.com/yechengxi/DBG2OLC/raw/master/Manual.docx

To use precompiled versions,please go to:

https://github.com/yechengxi/DBG2OLC/tree/master/compiled

To compile, download the source code of each program into separate folders and use: 

g++ -O3 -o SparseAssebmler *.cpp

g++ -O3 -o DBG2OLC *.cpp

g++ -O3 -o Sparc *.cpp


#Commands for Hybrid Assembly:

It is important to follow the steps before you invent your own pipeline. As combining DBG2OLC with other software tools may produce worse results. The reason is that many existing assembly/error correction tools introduce errors.  

##Step0. [Optional] Preparations:

We have provided code to help you select a subset of the reads:

https://github.com/yechengxi/AssemblyUtility

After compilation, you can use the following command to select a subset of reads from fasta/fastq files. Note that longest 0 is used here, if you set it to 1 it will select the longest reads.

./SelectLongestReads sum 600000000 longest 0 o Illumina_50x.fastq f Illumina_500bp_2x300_R1.fastq

./SelectLongestReads sum 260000000 longest 0 o Pacbio_20x.fasta f Pacbio.fasta

And you can use the following command to evaluate an assembly.

./AssemblyStatistics contigs YourAssembly.fasta

The program will generate two txt files containing essential statistics about your assembly.

Step1. Use an accurate DBG-assembler to construct short but accurate contigs. Please make sure they are the raw DBG contigs without using repeat resolving techniques such as gap closing or scaffolding. Otherwise you may have poor final results due to the errors introduced by the heuristics used in short read assembly pipelines.

SparseAssembler command format:

./SparseAssembler GS GENOME_SIZE NodeCovTh FALSE_KMER_THRESHOLD EdgeCovTh FALSE_EDGE_THRESHOLD k KMER_SIZE g SKIP_SIZE f YOUR_FASTA_OR_FASTQ_FILE1 f YOUR_FASTA_OR_FASTQ_FILE2 f YOUR_FASTA_OR_FASTQ_FILE3_ETC

A complete example on the S.cer w303 dataset:

Download the Illumina reads from

ftp://qb.cshl.edu/schatz/ectools/w303/Illumina_500bp_2x300_R1.fastq.gz

or here:

http://pan.baidu.com/s/1sk9bEnN

Normally with ~50x coverage, NodeCovTh 1 EdgeCovTh 0 can produce nice results.

./SparseAssembler LD 0 k 51 g 15 NodeCovTh 1 EdgeCovTh 0 GS 12000000 f ../Illumina_data/Illumina_50x.fastq

In this test run, the N50 is 29 kbp. As we have selected the beginning part of the sequencing file, which can be of lower quality, the next step may help to improve the assembly quality.

##[Miscellaneous]
For other more complex genomes or a different coverage, the first run may not generate the best result. The previous computations can be loaded and two parameters can be fine-tuned to construct a cleaner de Bruijn/ k-mer graph:

./SparseAssembler LD 1 NodeCovTh 2 EdgeCovTh 1 k 51 g 15 GS 12000000 f ../Illumina_data/Illumina_50x.fastq

The N50 is improved to 32kbp in my run.

The output Contigs.txt will be used by DBG2OLC.

##Step2. Overlap and layout.

Feed DBG2OLC with the contig file in fasta format from the previous step (Contigs.txt in this example). 

Download the PacBio reads from:

ftp://qb.cshl.edu/schatz/ectools/w303/Pacbio.fasta.gz

or here:

http://pan.baidu.com/s/1sk9bEnN

The basic command format of DBG2OLC is:

./DBG2OLC k KmerSize AdaptiveTh THRESH_VALUE1 KmerCovTh THRESH_VALUE2 MinOverlap THRESH_VALUE3 Contigs NGS_CONTIG_FILE f LONG_READS.FASTA RemoveChimera 1 

In the following example, the first 20x PacBio reads are extracted from the abovementioned file and we can assemble with:

./DBG2OLC k 17 AdaptiveTh 0.0001 KmerCovTh 2 MinOverlap 20 RemoveChimera 1 Contigs Contigs.txt f ../Pacbio_data/Pacbio _20x.fasta 

In our test run, the N50 is 583kbp.

There are three major parameters that affect the assembly quality:

M = matched k-mers between a contig and a long read.

AdaptiveTh: adaptive k-mer matching threshold. If M < AdaptiveTh* Contig_Length, this contig cannot be used as an anchor to the long read.

KmerCovTh: fixed k-mer matching threshold. If M < KmerCovTh, this contig cannot be used as an anchor to the long read.

MinOverlap: minimum overlap score between a pair of long reads.
For each pair of long reads, an overlap score is calculated by aligning the compressed reads and score with the matching k-mers.
 
##[Miscellaneous]

At this point, the parameters may be fine-tuned to get better performance. As with SparseAssembler, LD 1 can be used to load the compressed reads/anchored reads. 

Suggested tuning range is provided here:

For 10x/20x PacBio data: KmerCovTh 2-5, MinOverlap 10-30, AdaptiveTh 0.001~0.01.

For 50x-100x PacBio data: KmerCovTh 2-10, MinOverlap 50-150, AdaptiveTh 0.01-0.02. 

Some other less flexible or less important parameters:

k: k-mer size, 17 works well.

Contigs: the fasta contigs file from existing assembly.

MinLen: minimum read length. 

RemoveChimera: remove chimeric reads in the dataset, suggest 1 if you have >10x coverage. 

For high coverage data (100x), there are two other parameters:

ChimeraTh: default: 1, set to 2 if coverage is ~100x.

ContigTh: default: 1, set to 2 if coverage is ~100x.

These two are used in multiple alignment to remove problematic reads and false contig anchors. When we have high coverage, some more stringent conditions shall be applied as with the suggested parameters.

##Step 3. Call consensus. 

Install blasr and the consensus module (sparc/pbdagcon). Make sure they are in your path variable. 
The input files for consensus are: 

(1) backbone_raw.fasta by DBG2OLC

(2) DBG2OLC_Consensus_info.txt by DBG2OLC

(3) DBG contigs (in fasta format)

(4) PacBio reads (in fasta format)

You can check the N50 of (1) to see if you are satisfied, otherwise keep tuning and don¡¯t proceed.

// this is to concatenate the contigs and the raw reads for consensus

cat Contigs.txt pb_reads.fasta > ctg_pb.fasta

// we need to open a lot of files to distribute the above file into lots of smaller files

ulimit -n unlimited

//run the consensus scripts

sh ./split_and_run_sparc.sh backbone_raw.fasta DBG2OLC_Consensus_info.txt ctg_pb.fasta ./consensus_dir 2 >cns_log.txt

#Commands used to assemble other genomes:

The A. thaliana Ler-0 dataset:

20x PacBio reads:

./DBG2OLC KmerCovTh 2 AdaptiveTh 0.005 MinOverlap 20 RemoveChimera 1 Contigs Contigs.txt k 17 f ../PacBio/20x.fasta

40x PacBio reads:

./DBG2OLC  KmerCovTh 2 AdaptiveTh 0.01 MinOverlap 20 RemoveChimera 1 Contigs Contigs.txt k 17 f ../PacBio/40x.fasta

The H. sapiens dataset:

Longest 30x PacBio reads:

./DBG2OLC k 17 KmerCovTh 2 MinOverlap 20 AdaptiveTh 0.01 RemoveChimera 1 Contigs Contigs.txt f 30x.fasta >DBG2OLC_LOG.txt


#Commands for Non-hybrid NGS Assembly:


The program command is slightly different for purly Illumina reads assembly.  

Example command: 
./DBG2OLC LD 0 MinOverlap 70 PathCovTh 3 Contigs Contigs.txt k 31 KmerCovTh 0 f ReadsFile1.fa f ReadsFile2.fq f MoreFiles.xxx

There are four critical parameters:

k: k-mer length (max size: 31).

KmerCovTh: explained above, suggest 0-1.

MinOverlap: explained before.

PathCovTh: the minimum occurrence for a compressed read for a compressed read to be used, suggest  1-3.

Assembly is reported as DBG2OLC_Consensus.fasta.

The command we used for E. coli Illumina Miseq dataset:

./DBG2OLC k 31 PathCovTh 2 MinLen 50 MinOverlap 31 Contigs Contigs.txt KmerCovTh 0 f reads.fasta 


