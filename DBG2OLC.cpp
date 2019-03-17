#include "iostream"
#include "stdio.h"
#include "string"
#include "vector"
#include "cstdlib"
#include "bitset"
#include <map>
#include <math.h>
#include "memory"
#include <algorithm>
#include "fstream"
#include "sstream"
#include "list"
#include "stdlib.h"
#include "time.h"
#include "stdint.h"
#include "BasicDataStructure.h"
#include "GraphConstruction.h"
#include "GraphSimplification.h"
#include "GraphSearch.h"
#include "BuildContigs.h"
using namespace std;
/*
The goal here is to take a set of illumina contigs and re-index them using small k-mers.
done
For each long read, we find the contigs that are covered and then we build another contig-long_read index
partially done
Based on this new index we analyze the overlap relation and glue together the long reads.
idea: sort the a-stats of the ctgs and start from the best ctg to glue the best pacbio reads
*/

/*
bug to fix:
tips removal,bubble merging to be improved
alignment to be improved
graph structure to be improved
*/

int main(int argc, char* argv[])
{
	/*
	test
	string str0 = "";
	stringstream ss(str0);
	int coord = -1;
	ss >> coord;
	*/

	int K_size = 17, KmerCovTh = 0, ContigCovTh = 1, ChimeraTh = 1;
	double Redundancy = 1.1;
	bool CloseGaps = 0;
	double AdaptiveTh = -0.005;
	bool SparseAlign = 1;
	double mm_ratio = 0.5, mm_penalty = 0.3, mm_FR = 3.0, min_ext=1.0;
	int max_readlen = 1000000000;
	int max_dist = 5;
	int max_depth = 5;
	int ConsensusBestN = 10;
	int  MinOverlap = 50;
	int MaxOffset = 10000;
	int MaxMismatchFR = 5000;
	int MatchMethod = 2;
	int ScoreMethod = 3;
	string ctg_filename = "Contigs.txt", ref_filename = "RefContigIndex_info.txt", qry_filename, qry_filename2, contig_graph_filename = "RefinedContigGraph.txt";
	string Prefix = "LongRead";
	vector<string> qry_filename_vec;
	vector<string> qry_filename_vec2;
	bool LOAD_CONTIG_KMER_INDEX = 0, LOAD_LONGREAD_CONTIG_INDEX = 0, LOAD_OVERLAP_GRAPH = 0;
	bool ConstructBackbonesOnly=0;
	bool ASSEMBLY = 1, ALIGN = 0;
	bool CORRECT_WITH_CONTIG_GRAPH = 1;
	bool RemoveChimera = 0, RecoverFalseNegatives = 0;
	bool OUTPUT_NON_CONTAINED_READS = 0;
	bool OutputGraph = 1;
	bool FixOrientation = 0;
	bool Debug = 0;
	bool MSA = 0;
	int CleanRound = 1;
	int TopOverlaps = 50;
	int PathCovTh = 0;
	int64_t TotalReads = -1;
	bool InwardPairs = 0, OutwardPairs = 0;
	string ctg_kmer_idx_name = "ContigKmerIndex_";
	vector<string> LongReadFilenameVec,LongReadsInfoFiles,NonContainedReadsFiles;
	bool HELP = 1;
	bool Consensus = 1;

	
	int MinLen = 1;
	int MaxLen = 150;// max compressed read length
	for (int i = 1; i < argc; ++i)
	{


		if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "h") == 0)
		{
			i++;
			HELP = 1;
			break;
		}

		if (strcmp(argv[i], "Contigs") == 0)
		{
			i++;
			ctg_filename = (argv[i]);

			continue;
		}
		if (strcmp(argv[i], "ContigGraph") == 0)
		{
			i++;
			contig_graph_filename = (argv[i]);

			continue;
		}
		if (strcmp(argv[i], "SparseAlign") == 0)
		{
			i++;
			SparseAlign = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "r") == 0)
		{
			i++;
			ref_filename = (argv[i]);

			continue;
		}
		if (strcmp(argv[i], "R") == 0)
		{
			i++;
			TotalReads = atoi(argv[i]);

			continue;
		}
		
		if (strcmp(argv[i], "CloseGaps") == 0)
		{
			i++;
			CloseGaps = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "BB") == 0)
		{
			i++;
			ConstructBackbonesOnly = atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "FixOrientation") == 0)
		{
			i++;
			FixOrientation = atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "ConsensusBestN") == 0)
		{
			i++;
			ConsensusBestN = atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "q") == 0)
		{
			i++;
			qry_filename = (argv[i]);
			qry_filename_vec.push_back(qry_filename);
			continue;
		}
		if (strcmp(argv[i], "qr") == 0)
		{
			i++;
			qry_filename2 = (argv[i]);
			qry_filename_vec2.push_back(qry_filename2);
			continue;
		}
		if (strcmp(argv[i], "f") == 0)
		{
			i++;
			LongReadFilenameVec.push_back(argv[i]);
			string name1= argv[i];
			int n = 0;
			for (n = name1.size() - 1; n >= 0; --n)
			{
				if (name1[n] == '\\'||name1[n]=='/')
				{
					break;
				}
			}
			name1 = name1.substr(n + 1, name1.size());

			name1="ReadsInfoFrom_" + name1;
			LongReadsInfoFiles.push_back(name1);
			name1 = argv[i];
			n = 0;
			for (n = name1.size() - 1; n >= 0; --n)
			{
				if (name1[n] == '\\' || name1[n] == '/')
				{
					break;
				}
			}
			name1 = name1.substr(n + 1, name1.size());

			name1 = "NonContainedReadsFrom_" + name1;
			NonContainedReadsFiles.push_back(name1);
			continue;
		}

		if (strcmp(argv[i], "k") == 0)
		{
			i++;
			K_size = atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "Match") == 0)
		{
			i++;
			MatchMethod = atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "Score") == 0)
		{
			i++;
			ScoreMethod = atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "ip") == 0)
		{
			i++;
			InwardPairs = 1;
			qry_filename_vec.push_back(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "ipr") == 0)
		{
			i++;
			qry_filename_vec2.push_back(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "NC") == 0)
		{
			i++;
			OUTPUT_NON_CONTAINED_READS = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "Consensus") == 0)
		{
			i++;
			Consensus = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "TopOverlaps") == 0)
		{
			i++;
			TopOverlaps = atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "PathCovTh") == 0)
		{
			i++;
			PathCovTh = atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "AdaptiveTh") == 0)
		{
			i++;
			AdaptiveTh = atof(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "mm_penalty") == 0)
		{
			i++;
			mm_penalty = atof(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "mm_FR") == 0)
		{
			i++;
			mm_FR = atof(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "min_ext") == 0)
		{
			i++;
			min_ext = atof(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "mm_ratio") == 0)
		{
			i++;
			mm_ratio = atof(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "max_depth") == 0)
		{
			i++;
			max_depth = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "max_dist") == 0)
		{
			i++;
			max_dist = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "MSA") == 0)
		{
			i++;
			MSA = atof(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "CleanRound") == 0)
		{
			i++;
			CleanRound = atoi(argv[i]);

			continue;
		}
	
		if (strcmp(argv[i], "Debug") == 0 || strcmp(argv[i], "DEBUG") == 0)
		{
			i++;
			Debug = atoi(argv[i]);

			continue;
		}

		if (strcmp(argv[i], "MinLen") == 0)
		{
			i++;
			MinLen = atoi(argv[i]);

			continue;
		}

		if (strcmp(argv[i], "MaxLen") == 0)
		{
			i++;
			MaxLen = atoi(argv[i]);

			continue;
		}

		if (strcmp(argv[i], "MinOverlap") == 0)
		{
			i++;
			MinOverlap = atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "MaxOffset") == 0)
		{
			i++;
			MaxOffset = atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "MaxMismatchFR") == 0)
		{
			i++;
			MaxMismatchFR = atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "LD0") == 0)
		{
			i++;
			LOAD_CONTIG_KMER_INDEX = (bool)atoi(argv[i]);
			//LOAD_LONGREAD_CONTIG_INDEX = LOAD_CONTIG_KMER_INDEX;
			continue;
		}
		if (strcmp(argv[i], "LD1") == 0 || strcmp(argv[i], "LD") == 0)
		{
			i++;
			LOAD_LONGREAD_CONTIG_INDEX = (bool)atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "LD2") == 0)
		{
			i++;
			LOAD_OVERLAP_GRAPH = (bool)atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "RemoveChimera") == 0)
		{
			i++;
			RemoveChimera = (bool)atoi(argv[i]);
			
			continue;
		}

		if (strcmp(argv[i], "RecoverFalseNegatives") == 0)
		{
			i++;
			RecoverFalseNegatives = (bool)atoi(argv[i]);

			continue;
		}
		

		if (strcmp(argv[i], "Assembly") == 0)
		{
			i++;
			ASSEMBLY = (bool)atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "Prefix") == 0)
		{
			i++;
			Prefix = (argv[i]);

			continue;
		}
		if (strcmp(argv[i], "Align") == 0)
		{
			i++;
			ALIGN = (bool)atoi(argv[i]);

			continue;
		}
		
		if (strcmp(argv[i], "KmerCovTh") == 0)
		{
			i++;
			KmerCovTh = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "Redundancy") == 0)
		{
			i++;
			Redundancy = atof(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "ContigCovTh") == 0)
		{
			i++;
			ContigCovTh = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "ChimeraTh") == 0)
		{
			i++;
			ChimeraTh = atoi(argv[i]);
			continue;
		}
	
	}

	if (HELP)
	{
		cout << " Example command: " << endl;
		
		cout << "For third-gen sequencing: DBG2OLC LD1 0 Contigs contig.fa k 17 KmerCovTh 2 MinOverlap 20 AdaptiveTh 0.005 f reads_file1.fq/fa f reads_file2.fq/fa" << endl;
		cout << "For sec-gen sequencing: DBG2OLC LD1 0 Contigs contig.fa k 31 KmerCovTh 0 MinOverlap 50 PathCovTh 1 f reads_file1.fq/fa f reads_file2.fq/fa" << endl;

		cout << "Parameters:" << endl;
		cout << "MinLen: min read length for a read to be used." << endl;
		cout << "Contigs:  contig file to be used." << endl; 
		cout << "k: k-mer size." << endl;
		cout << "LD: load compressed reads information. You can set to 1 if you have run the algorithm for one round and just want to fine tune the following parameters." << endl;
		cout << "PARAMETERS THAT ARE CRITICAL FOR THE PERFORMANCE:" << endl;
		cout << "If you have high coverage, set large values to these parameters." << endl;
		cout << "KmerCovTh: k-mer matching threshold for each solid contig. (suggest 2-10)" << endl;
		cout << "MinOverlap: min matching k-mers for each two reads. (suggest 10-150)" << endl;
		cout << "AdaptiveTh: [Specific for third-gen sequencing] adaptive k-mer threshold for each solid contig. (suggest 0.001-0.02)" << endl;		
		cout << "PathCovTh: [Specific for Illumina sequencing] occurence threshold for a compressed read. (suggest 1-3)" << endl;
		cout << "Author: Chengxi Ye cxy@umd.edu." << endl;
		cout << "last update: Jun 11, 2015." << endl;
		
	}
	
	//
	if (RemoveChimera)
	{
		MSA = 1;
	}
	contigs_info contigs_info;
	contigs_info.KmerCovTh = KmerCovTh;
	contigs_info.ContigsFileName = ctg_filename;
	if (AdaptiveTh > 0.0)
	{
		contigs_info.CallConsensus = Consensus;
	}
	else
	{
		contigs_info.CallConsensus = Consensus;
	}



	contigs_info.ContigCovTh = ContigCovTh;
	contigs_info.ContigGraphName = contig_graph_filename;
	ifstream ctgs_in(ctg_filename.c_str());
	contigs_info.contig_tag.clear();
	contigs_info.contig_tag.push_back("");
	string str, contig_str;
	size_t hashTableSZ = 0;
	struct hashtable ht;
	ht.ht_sz = 0;


//	string LongReadContigIndex_name = Prefix + "ContigIndex_info.txt";
	if (ASSEMBLY)
	if (LOAD_LONGREAD_CONTIG_INDEX == 0)
	if (LOAD_CONTIG_KMER_INDEX)
	{
		cout << "Loading kmer index..." << endl;
		LoadContigKmerIndex(&ht, ctg_kmer_idx_name);
		cout << "Done." << endl;
	}
	else
	{
		cout << "Loading contigs." << endl;
		int max_ctg_len = -1;
		int tot_ctgs = 0;
		string tag, n_tag;
		
		while (get_a_fasta_read(ctgs_in, tag, contig_str, n_tag))
		{
			contigs_info.contig_tag.push_back(tag);
			int sz = contig_str.size();
			hashTableSZ += sz - K_size + 1;
			if (max_ctg_len < sz)
			{
				max_ctg_len = sz;
			}
			tot_ctgs++;

		}

		ctgs_in.close();
		contigs_info.contig_sz_vt.clear();
		contigs_info.contig_sz_vt.resize(tot_ctgs + 1);

		ref_read_t ref;

		ref.read_bits = (uint64_t *)myMalloc((size_t)(max_ctg_len / 4) + 100);
		ref.alloc_sz = (size_t)(max_ctg_len / 4 + 100);

		Init_HT(&ht, hashTableSZ);
		
		for (int round = 1; round <= 2; ++round)
		{

		
			ht.round = round;
			int bucket_count = 0;
			ctgs_in.open(ctg_filename.c_str());
			tot_ctgs = 0;
			while (get_a_fasta_read(ctgs_in, tag, contig_str, n_tag))
			{
				int sz = contig_str.size();
				tot_ctgs++;
				ref.contig_no = tot_ctgs;
				contigs_info.contig_sz_vt[tot_ctgs] = sz;
				Init_Ref_Read(contig_str, ref);
				//if (contig_str.size() > 100)
				if (1)
				{
					Contig_Kmer_Index(&ref, &ht, K_size, &bucket_count);

				}
			}
			ctgs_in.close();
			if (ht.round == 1)
			{
				SwitchBuckets(&ht);
			}
			contigs_info.total_contigs = tot_ctgs;
			cout << bucket_count << " k-mers in round " << round << "." << endl;
		}
		SaveContigKmerIndex(&ht, ctg_kmer_idx_name);
	}


	
	contigs_info.contig_sz_vt.clear();
	contigs_info.contig_sz_vt.push_back(0);
	int tot_ctgs = 0, ContigLen = 0;
	string tag, n_tag;
	ctgs_in.close();

	ctgs_in.open(ctg_filename.c_str());
	contigs_info.contig_tag.clear();

	contigs_info.contig_tag.push_back("");
	while (get_a_fasta_read(ctgs_in, tag, contig_str, n_tag))
	{
		int ContigLen = contig_str.size();
		tot_ctgs++;
		contigs_info.contig_tag.push_back(tag);
		contigs_info.contig_sz_vt.push_back(ContigLen);
	}
	contigs_info.total_contigs = tot_ctgs;
	ctgs_in.close();

	



	/*
	LoadRefinedContigGraph(&contigs_info,contig_graph_filename);


	str.clear();
	ctgs_in.open(ctg_filename.c_str());
	contigs_info.contigs_str.push_back(str);

	while (getline(ctgs_in, str) && str.size() > 0)
	{
	getline(ctgs_in, contig_str);
	while (contig_str[contig_str.size() - 1] == '\n' || contig_str[contig_str.size() - 1] == '\r')
	{
	int sz = contig_str.size();
	contig_str[sz - 1] = '\0';
	contig_str.resize(sz - 1);
	}
	int sz = contig_str.size();
	contigs_info.contigs_str.push_back(contig_str);
	tot_ctgs++;

	}
	ctgs_in.close();
	int cov = contigs_info.contig_adjacency_right[abs(17303)][-26604].cov;
	OutputFastaContigs(&contigs_info, "ContigsByLongReadsNodes_info.txt");
	*/



	if (!LOAD_LONGREAD_CONTIG_INDEX&&LongReadFilenameVec.size())
	{
		ref_read_t long_read;

		if (LongReadFilenameVec.size() == 0)
		{
			cout << "No input files." << endl;
			return 0;
		}
		
		string Histogram_name = Prefix + "ContigIndex_Histogram.txt";
		ofstream LongReadContigIndex_hist(Histogram_name.c_str());
		ofstream LongReadContigIndex_log("LongReadContigIndex_log.txt");
		ofstream out_refined_graph;
			
		//ofstream out_positions_debug("CenterEstimates.txt");
		ofstream NonContainedReads_o;
		if (AdaptiveTh < 0.0)
		{
			OUTPUT_NON_CONTAINED_READS = 1;
			out_refined_graph.open("RefinedContigGraph.txt");
		}
		
		long_read.read_bits = (uint64_t *)myMalloc((size_t)(max_readlen / 4) + 100);
		long_read.alloc_sz = (size_t)(max_readlen / 4 + 100);
		cout << "Analyzing reads..." << endl;
		time_t beg_time, end_time;
		time(&beg_time);
		map<int, int> contigs_match_hist;
		uint64_t total_kmers = 0, MatchingUniqueKmers = 0;
		int read_idx = 0;
		for (int i = 0; i < LongReadFilenameVec.size(); ++i)
		{
			cout << "File" << i + 1 << ": " << LongReadFilenameVec[i] << endl;
			ofstream LongReadContigIndex_info(LongReadsInfoFiles[i].c_str());
			if (OUTPUT_NON_CONTAINED_READS)
			{
				NonContainedReads_o.close();

				NonContainedReads_o.open(NonContainedReadsFiles[i].c_str());
			}
			bool fq_flag = 0;
			string tag, n_tag, read_str, quality_score;

			ifstream LongReads_in(LongReadFilenameVec[i].c_str());

			getline(LongReads_in, tag);
			if (tag[0] == '@')
			{
				fq_flag = 1;
			}
			LongReads_in.close();

			LongReads_in.open(LongReadFilenameVec[i].c_str());

			bool read_success = 1;

			while (read_success)
			{
				if (fq_flag == 0)
				{
					read_success = get_a_fasta_read(LongReads_in, tag, read_str, n_tag);
				}
				else
				{
					read_success = get_a_fastq_read(LongReads_in, tag, read_str, quality_score);
				}
				if (read_success == 0)
				{
					break;
				}


				//deal with 'N's in a stupid way
				string temp_str, raw_read = read_str;
				temp_str.resize(read_str.size());
				size_t nn = 0;
				for (int n = 0; n < read_str.size(); ++n)
				{
					if (read_str[n] != 'N')
					{
						temp_str[nn]=read_str[n];
						nn++;
					}
				}
				temp_str.resize(nn);
				read_str = temp_str;
				if (read_str.size() > max_readlen)
				{
					read_str.resize(max_readlen - 1);
				}
				if (read_str.size() < MinLen)
				{
					continue;
				}

				read_idx++;
				//cout << tag << endl;
				Init_Ref_Read(read_str, long_read);
				long_read.read_idx = read_idx;
				struct LongReadContigIndex LongReadContigIndex;
				LongReadContigIndex.KmerCovTh = KmerCovTh;
				LongReadContigIndex.nMatches = 0;
				LongReadContigIndex.BlockSize = 1000;
				LongReadContigIndex.FastMap = 1;

				ReadCompression(&long_read, &ht, &contigs_info, K_size, &LongReadContigIndex);

				MatchingUniqueKmers += LongReadContigIndex.nMatches;
				total_kmers += long_read.readLen - K_size + 1;


				map<int, ContigInRead>::iterator tmp_it;
				int contig_matches = 0;

				for (tmp_it = LongReadContigIndex.layout.begin(); tmp_it != LongReadContigIndex.layout.end(); ++tmp_it)
				{

					if (tmp_it->second.cov > KmerCovTh)
					{
						contig_matches++;

					}
				}

				if (contig_matches > 1||AdaptiveTh>0.0)
				{
					LongReadContigIndex_info << tag << endl;
					LongReadContigIndex_info << read_str.size() << endl;
					for (tmp_it = LongReadContigIndex.layout.begin(); tmp_it != LongReadContigIndex.layout.end(); ++tmp_it)
					{
						if (tmp_it->second.cov > KmerCovTh)
						{
							LongReadContigIndex_info << tmp_it->first << ", " << tmp_it->second.ctg_no << ", " << tmp_it->second.cov << ", " << tmp_it->second.coord2 << endl;
						}
					}
					contigs_match_hist[contig_matches]++;
					if (OUTPUT_NON_CONTAINED_READS)
					{
						NonContainedReads_o << tag << endl;
						NonContainedReads_o << raw_read << endl;
					}
				}
			


				//newly added...
				if (OUTPUT_NON_CONTAINED_READS)
				{
					map<int, KmerInContig>::iterator it1, it2;

					for (it1 = LongReadContigIndex.LR2CTG.begin(); it1 != LongReadContigIndex.LR2CTG.end(); ++it1)
					{

						it2 = it1;
						it2++;
						if (it2 == LongReadContigIndex.LR2CTG.end())
						{
							break;
						}
						if (it2->second.contig_no != it1->second.contig_no)
						{
							bool flip1, flip2;
							flip1 = it1->second.flip;
							flip2 = it2->second.flip;
							int pos1 = it1->second.pos;
							int pos2 = it2->second.pos;
							int extra_bases1, extra_bases2;
							int contig1 = it1->second.contig_no;
							int contig2 = it2->second.contig_no;

							if (flip1)
							{
								contig1 = -contig1;
							}
							if (flip2)
							{
								contig2 = -contig2;
							}
							int dist = it2->first - it1->first;
							if (flip1 == 0)
							{
								extra_bases1 = contigs_info.contig_sz_vt[abs(contig1)] - pos1;
							}
							else
							{
								extra_bases1 = pos1;
							}

							if (flip2 == 0)
							{
								extra_bases2 = pos2;
							}
							else
							{
								extra_bases2 = contigs_info.contig_sz_vt[abs(contig2)] - pos2;
							}
							int extra_bases = extra_bases1 + extra_bases2;
							dist = dist - extra_bases;
							if (contig1 != 0 && contig2 != 0)
							{
								if (contig1 > 0)
								{
									out_refined_graph << contig1 << " " << contig2 << " " << dist << endl;
								}
								else
								{
									out_refined_graph << contig1 << " " << -contig2 << " " << dist << endl;
								}
							}

						}

					}

					//newly added...
				}



			}

		}
		cout << "Long reads indexed. " << endl;
		for (map<int, int>::iterator hist_it = contigs_match_hist.begin(); hist_it != contigs_match_hist.end(); ++hist_it)
		{
			LongReadContigIndex_hist << hist_it->first << ", " << hist_it->second << endl;
		}
		cout << "Total Kmers: " << total_kmers << endl;

		cout << "Matching Unique Kmers: " << MatchingUniqueKmers << endl;
		LongReadContigIndex_log << "Total Kmers: " << total_kmers << endl;
		LongReadContigIndex_log << "Matching Unique Kmers: " << MatchingUniqueKmers << endl;

		time(&end_time);
		cout <<"Compression time: "<< difftime(end_time, beg_time) << " secs." << endl;

	}

	if (ht.ht_sz > 0)
	{
		FreeSparseKmerGraph(&ht);

	}

	reads_info reads_info;

	align_info align_info;
	align_info.SparseAlign = SparseAlign;
	align_info.Debug = Debug;
	align_info.flip = 0;
	align_info.force_flip = 0;
	align_info.fix_orientation = 1;

	align_info.max_mismatch_FR_score = MinOverlap * mm_FR;
	if (align_info.max_mismatch_FR_score > 300)
	{
		//align_info.max_mismatch_FR = 300;
	}
	align_info.min_overlap = MinOverlap;
	align_info.min_extension_FR = MinOverlap * min_ext;
	align_info.mm_ratio = mm_ratio;// 0.4;//0.4
	align_info.mm_penalty = mm_penalty;// 0.
	align_info.max_len = 250;//modified on 20190315 to allow aligning ultra long reads
	align_info.band_width = 20;
	reads_info.K_size = K_size;
	reads_info.CloseGaps = CloseGaps;
	reads_info.ConsensusBestN = ConsensusBestN;
	align_info.MaxOffset = MaxOffset;
	align_info.MaxMismatchFR = MaxMismatchFR;
	reads_info.Clean = 1;
	if (AdaptiveTh < 0.0)
	{
		MatchMethod = 1;//not ready for sec gen yet
	}

	align_info.matching_method = MatchMethod;
	cout << "Scoring method: " << ScoreMethod << endl;
	cout << "Match method: " << MatchMethod << endl;
	align_info.scoring_method = ScoreMethod;

	reads_info.AdaptiveTh = AdaptiveTh;
	reads_info.input_files = LongReadFilenameVec;
	reads_info.TopOverlaps = TopOverlaps;
	reads_info.PathCovTh = PathCovTh;
	reads_info.ChimeraTh = ChimeraTh;
	reads_info.RemoveChimera = RemoveChimera;
	//cout << "Recover false negatives: " << RecoverFalseNegatives << endl;
	reads_info.RecoverFalseNegatives = RecoverFalseNegatives;
	reads_info.MSA = MSA;
	reads_info.Debug = Debug;
	reads_info.MinLen = MinLen;
	reads_info.max_reads = 500;
	reads_info.TotalReads = TotalReads;
	


	reads_info.KmerCovTh = KmerCovTh;
	reads_info.Redundancy = Redundancy;
	if (ALIGN)
	{
		ofstream o_align("align2ref.txt");
		ofstream o_align_rc("align2rc.txt");
		ofstream o_align_rc_tags("align2rc_tags.txt");
		ofstream o_align_chimera("chimera.txt");
		ofstream o_align_chimera_tags("chimera_tags.txt");
		ofstream o_corr_len("ContigLenEstCorr.txt");

		align_info.fix_orientation = 1;
		//align_info.InnerProduct = 0;


		cout << "Loading reference." << endl;

		//string RefContigIndex_name = "RefContigIndex_info.txt";
		LoadingRefIndex(ref_filename, &reads_info, &contigs_info);
		cout << "Reference loaded. Size: " << reads_info.RefIndexVec.size() - 1 << endl;
		LoadLongReadIndex(qry_filename_vec, &reads_info, &contigs_info);
		cout << "Queries loaded. Size: " << reads_info.LongReadIndexVec.size() - 1 << endl;

		for (int i = 1; i != reads_info.LongReadIndexVec.size(); ++i)
		{

			vector<Coord_CTG_Cov> qry = reads_info.LongReadIndexVec[i];
			align_maps align_maps;

			vector<struct align_info> align_info_vec;
			align_info_vec.push_back(align_info);

			align2ref(&contigs_info, &reads_info, qry, &align_maps, align_info_vec);
			//cout << align_info_vec.size() << endl;
			for (int v = 0; v < align_info_vec.size(); ++v)//
			{
				struct align_info align_info_tmp = align_info_vec[v];
				o_align << ">Read_" << i << "_ref_" << v + 1 << "(" << align_info_tmp.ref_idx << ")_score:" << align_info_tmp.max_score << endl;
				for (int j = 0; j != align_info_tmp.ref_aligned.size(); ++j)
				{
					o_align << align_info_tmp.ref_aligned[j] << ", ";
				}
				o_align << endl;
				for (int j = 0; j != align_info_tmp.qry_aligned.size(); ++j)
				{
					o_align << align_info_tmp.qry_aligned[j] << ", ";
				}
				o_align << endl;

				o_corr_len << align_info_tmp.max_score << endl;
			}
			if (align_info_vec.size() > 1)
			{
				bool rc_match = 0;
				if (align_info_vec[0].ref_idx == -align_info_vec[1].ref_idx)
				{

					//if ((align_info_vec[0].max_match_qry < align_info_vec[1].min_match_qry))// || (align_info_vec[1].max_match_qry < align_info_vec[0].min_match_qry))//non overlap match
					if ((align_info_vec[0].max_match_ref < align_info_vec[1].min_match_ref))// || (align_info_vec[1].max_match_qry < align_info_vec[0].min_match_qry))//non overlap match

					{
						map<int, int> local_index_qry;

						for (int cc = 0; cc <= align_info_vec[0].max_match_ref; ++cc)
						{
							local_index_qry[qry[cc].contig_no] = 1;
						}
						for (int cc = align_info_vec[1].min_match_ref; cc <= align_info_vec[1].max_match_ref; ++cc)
						{
							if (local_index_qry.count(-qry[cc].contig_no))
							{
								rc_match = 1;
							}
						}

					}
					//					if ((align_info_vec[1].max_match_qry < align_info_vec[0].min_match_qry))//non overlap match
					if ((align_info_vec[1].max_match_ref < align_info_vec[0].min_match_ref))//non overlap match
					{
						map<int, int> local_index_qry;

						for (int cc = 0; cc <= align_info_vec[1].max_match_ref; ++cc)
						{
							local_index_qry[qry[cc].contig_no] = 1;
						}
						for (int cc = align_info_vec[0].min_match_ref; cc <= align_info_vec[0].max_match_ref; ++cc)
						{
							if (local_index_qry.count(-qry[cc].contig_no))
							{
								rc_match = 1;
							}
						}

					}
					if (rc_match)
					{
						for (int v = 0; v < align_info_vec.size(); ++v)//
						{
							struct align_info align_info_tmp = align_info_vec[v];
							o_align_rc << ">Read_" << i << "_ref_" << v + 1 << "(" << align_info_tmp.ref_idx << ")_score:" << align_info_tmp.max_score << endl;
							for (int j = 0; j != align_info_tmp.ref_aligned.size(); ++j)
							{
								o_align_rc << align_info_tmp.ref_aligned[j] << ", ";
							}
							o_align_rc << endl;
							for (int j = 0; j != align_info_tmp.qry_aligned.size(); ++j)
							{
								o_align_rc << align_info_tmp.qry_aligned[j] << ", ";
							}
							o_align_rc << endl;


						}

						if (reads_info.tag_vec.size() > 0)
						{
							o_align_rc_tags << reads_info.tag_vec[i] << endl;
						}
					}
				}

				if (rc_match == 0)
				{
					bool chimera_match = 0;
					if ((align_info_vec[0].max_match_qry < align_info_vec[1].min_match_qry) || (align_info_vec[1].max_match_qry < align_info_vec[0].min_match_qry))//non overlap match
					{
						chimera_match = 1;
						for (int v = 0; v < align_info_vec.size(); ++v)//
						{
							struct align_info align_info_tmp = align_info_vec[v];
							o_align_chimera << ">Read_" << i << "_ref_" << v + 1 << "(" << align_info_tmp.ref_idx << ")_score:" << align_info_tmp.max_score << endl;
							for (int j = 0; j != align_info_tmp.ref_aligned.size(); ++j)
							{
								o_align_chimera << align_info_tmp.ref_aligned[j] << ", ";
							}
							o_align_chimera << endl;
							for (int j = 0; j != align_info_tmp.qry_aligned.size(); ++j)
							{
								o_align_chimera << align_info_tmp.qry_aligned[j] << ", ";
							}
							o_align_chimera << endl;
							if (reads_info.tag_vec.size() > 0)
							{
								o_align_chimera_tags << reads_info.tag_vec[i] << endl;
							}
						}
					}
				}
			}

		}

	}

	if (ASSEMBLY)
	{


	

		
	
		if (ConstructBackbonesOnly)
		if (AdaptiveTh > 0.0)
		{

			cout << "debug." << endl;
			

			/*

			string final_graph_name = "RawOverlapGraph", raw_bog_name = "RawBestOverlapGraph", final_bog_name2 = "BiDirectedBestOverlapGraph.dot";

			vector<string> file_vt;
			file_vt.push_back("CleanedLongReads.txt");
			reads_info.Clean = 0;
			LoadLongReadIndex(file_vt, &reads_info, &contigs_info);//for long reads
			
			LoadOverlapGraph(&reads_info, final_graph_name);
			LoadBestOverlapGraph(&reads_info, raw_bog_name);
			LoadReadsInfo(&reads_info);
			


			//aggressive_cleaning(&reads_info);
			for (int clean_round = 1; clean_round <= CleanRound; ++clean_round)
			{
				reads_info.mode = 1;//1
				ConstructUndirectedBestOverlapGraph(&reads_info, final_bog_name2);//
				//return 0;
				reads_info.LeftOverlaps = reads_info.LeftBestOverlapsTemp;
				reads_info.RightOverlaps = reads_info.RightBestOverlapsTemp;
				//OutputOverlapGraph(&reads_info, final_graph_name);
				reads_info.n_deleted_edges = 0;
				if (clean_round > 0)
				{
					//aggressive_cleaning(&reads_info);


					bool RemoveTipsOnly = 0;
					cout << "Graph simplification." << endl;

					for (int iter = 0; iter < 4; ++iter)
					{
						int nBranches = 0, n_linear = 0;

						cout << "Iteration: " << iter << endl;
						if (iter % 2 == 0)
						{
							RemoveTipsOnly = 1;
						}
						else
						{
							RemoveTipsOnly = 0;
						}

						for (int i = 1; i < reads_info.LeftOverlaps.size(); ++i)
						{
							int beg_read = i;
							if (i == 9)
							{
								cout << "";
							}

							if (reads_info.RightOverlaps[i].size() == 1 && reads_info.LeftOverlaps[i].size() == 1)
							{
								n_linear++;
							}
							if (reads_info.RightOverlaps[i].size() > 1)
							{
								nBranches++;
								BFSearchBubbleRemoval_read(&reads_info, beg_read, max_depth, max_dist, RemoveTipsOnly);
							}

							if (reads_info.LeftOverlaps[i].size() > 1)
							{
								nBranches++;
								BFSearchBubbleRemoval_read(&reads_info, -beg_read, max_depth, max_dist, RemoveTipsOnly);
							}

						}

						cout << nBranches << " branching positions." << endl;
						cout << n_linear << " linear nodes." << endl;
						//aggressive_cleaning(&reads_info);
					}
					cout << reads_info.n_deleted_edges << " edges deleted." << endl;
					if (OutputGraph)
					{
						string clean_graph_filename = "CleanedGraph.dot";
						//OutputCleanedOverlapGraph(&reads_info, clean_graph_filename);

						OutputOverlapGraph(&reads_info, clean_graph_filename);
					}

					ConstructCleanedBestOverlapGraph(&reads_info, &contigs_info);

					aggressive_cleaning(&reads_info);
					if (OutputGraph)
					{
						string final_bog_name = "CleanedBestOverlapGraph";
						OutputBestOverlapGraph(&reads_info, final_bog_name);
					}

				}

			}



			*/














			
			
			string consensus_name = "backbone.fasta";
			
			string lr_name = "DBG2OLC_";;
			//LoadSparseOverlapGraph(&contigs_info, &reads_info, lr_name);

			
			cout << "Loading contigs." << endl;
			ctgs_in.close();
			ctgs_in.open(ctg_filename.c_str());
			int tot_ctgs = 0;
			contigs_info.contigs_str.clear();
			str.clear();
			contigs_info.contigs_str.push_back(str);
			contigs_info.contig_tag.clear();
			contigs_info.contig_tag.push_back(str);
			bool read_success = 1;
			string ctg_tag, ctg_tag_n;
			while (read_success)
			{

				read_success = get_a_fasta_read(ctgs_in, ctg_tag, contig_str, ctg_tag_n);
				contigs_info.contig_tag.push_back(ctg_tag);

				while (contig_str[contig_str.size() - 1] == '\n' || contig_str[contig_str.size() - 1] == '\r')
				{
					int sz = contig_str.size();
					contig_str[sz - 1] = '\0';
					contig_str.resize(sz - 1);
				}
				int sz = contig_str.size();
				contigs_info.contigs_str.push_back(contig_str);
				tot_ctgs++;

			}
			ctgs_in.close();
			
			cout << "Construct backbone sequence." << endl;

			ConstructBackbone(&reads_info, &contigs_info, &align_info);

			return 0;
		}



		align_info.fix_orientation = 0;

		if (qry_filename_vec.size() == 0)
		{
			for (int q = 0; q<LongReadsInfoFiles.size(); ++q)
			{
				qry_filename_vec.push_back(LongReadsInfoFiles[q]);
			}//default name

		}
		cout << "Loading long read index" << endl;
		if (AdaptiveTh > 0.0)
		{
			if (!LOAD_OVERLAP_GRAPH)
			{
				reads_info.Clean = 1;
				LoadLongReadIndex(qry_filename_vec, &reads_info, &contigs_info);//for long reads
			}
			else
			{
				vector<string> file_vt;
				file_vt.push_back("CleanedLongReads.txt");
				reads_info.Clean = 0;
				LoadLongReadIndex(file_vt, &reads_info, &contigs_info);//for long reads
			}
			align_info.FAST = 1;
		}
		else
		{
			LoadLongReadIndexWithMerging(qry_filename_vec, &reads_info, &contigs_info);//for short reads
			//LoadPairedReadIndexWithMerging(qry_filename_vec, &reads_info, &contigs_info);//for short reads

			align_info.FAST = 0;
		}



		cout << "Loaded." << endl;
		string final_graph_name = "RawOverlapGraph", raw_bog_name = "RawBestOverlapGraph", final_bog_name2 = "BiDirectedBestOverlapGraph.dot";

		if (!LOAD_OVERLAP_GRAPH)
		{
			align_info.fix_orientation = FixOrientation;//latest test

		
			BuildOverlapGraph(&reads_info, &contigs_info, &align_info);
			
			if (OutputGraph)
			{
				OutputOverlapGraph(&reads_info, final_graph_name);
				OutputBestOverlapGraph(&reads_info, raw_bog_name);
				OutputReadsInfo(&reads_info);
				LoadOverlapGraph(&reads_info, final_graph_name);
				LoadBestOverlapGraph(&reads_info, raw_bog_name);
				LoadReadsInfo(&reads_info);
			}


			//aggressive_cleaning(&reads_info);
			for (int clean_round = 1; clean_round <= CleanRound; ++clean_round)
			{
				reads_info.mode = 1;//1
				ConstructUndirectedBestOverlapGraph(&reads_info, final_bog_name2);//
				//return 0;
				reads_info.LeftOverlaps = reads_info.LeftBestOverlapsTemp;
				reads_info.RightOverlaps = reads_info.RightBestOverlapsTemp;
				//OutputOverlapGraph(&reads_info, final_graph_name);
				reads_info.n_deleted_edges = 0;
				if (clean_round > 0)
				{
					//aggressive_cleaning(&reads_info);


					bool RemoveTipsOnly = 0;
					cout << "Graph simplification." << endl;

					for (int iter = 0; iter < 4; ++iter)
					{
						int nBranches = 0, n_linear = 0;

						cout << "Iteration: " << iter << endl;
						if (iter % 2 == 0)
						{
							RemoveTipsOnly = 1;
						}
						else
						{
							RemoveTipsOnly = 0;
						}

						for (int i = 1; i < reads_info.LeftOverlaps.size(); ++i)
						{
							int beg_read = i;

							
							if (reads_info.RightOverlaps[i].size() == 1 && reads_info.LeftOverlaps[i].size() == 1)
							{
								n_linear++;
							}
							if (reads_info.RightOverlaps[i].size() > 1)
							{
								nBranches++;
								BFSearchBubbleRemoval_read(&reads_info, beg_read, max_depth, max_dist, RemoveTipsOnly);
							}

							if (reads_info.LeftOverlaps[i].size() > 1)
							{
								nBranches++;
								BFSearchBubbleRemoval_read(&reads_info, -beg_read, max_depth, max_dist, RemoveTipsOnly);
							}

						}

						cout << nBranches << " branching positions." << endl;
						cout << n_linear << " linear nodes." << endl;
						//aggressive_cleaning(&reads_info);
					}
					cout << reads_info.n_deleted_edges << " edges deleted." << endl;
					if (OutputGraph)
					{
						string clean_graph_filename = "CleanedGraph.dot";
						//OutputCleanedOverlapGraph(&reads_info, clean_graph_filename);

						OutputOverlapGraph(&reads_info, clean_graph_filename);
					}

					ConstructCleanedBestOverlapGraph(&reads_info, &contigs_info);

					aggressive_cleaning(&reads_info);
					if (OutputGraph)
					{
						string final_bog_name = "CleanedBestOverlapGraph";
						OutputBestOverlapGraph(&reads_info, final_bog_name);
					}

				}

			}




			string out_ctg_filename = "ContigsByLongReads.txt";
			reads_info.mode = 2;//2
			final_bog_name2 = "FinalBiDirectedBestOverlapGraph.dot";
			ConstructUndirectedBestOverlapGraph(&reads_info, final_bog_name2);//
			out_ctg_filename = "ContigsByLongReads";
			reads_info.LeftOverlaps = reads_info.LeftBestOverlapsTemp;
			reads_info.RightOverlaps = reads_info.RightBestOverlapsTemp;


			if (CloseGaps&&reads_info.AdaptiveTh > 0.0)
			{
				string graph_name = "OverlapGraphBeforeGapClosing";
				string bog_name = "BestOverlapGraphBeforeGapClosing";
				OutputOverlapGraph(&reads_info, graph_name);
				OutputBestOverlapGraph(&reads_info, bog_name);
				OutputReadsInfo(&reads_info);

				CloseGapsWithAllReads(&reads_info, &contigs_info, &align_info);

				reads_info.mode = 2;
				final_bog_name2 = "GapClosedUnDirectedBestOverlapGraph.dot";
				ConstructUndirectedBestOverlapGraph(&reads_info, final_bog_name2);//
				/*
				for (int i = 0; i < reads_info.LeftBestOverlapsTemp.size(); ++i)
				{
				if (reads_info.LeftBestOverlapsTemp[i].size() == 0 && reads_info.RightBestOverlapsTemp[i].size() == 0)
				{
				reads_info.contained[i] = 1;
				}
				}
				*/
				reads_info.LeftOverlaps = reads_info.LeftBestOverlapsTemp;
				reads_info.RightOverlaps = reads_info.RightBestOverlapsTemp;
				ConstructCleanedBestOverlapGraph(&reads_info, &contigs_info);

				/*
				ConstructUndirectedBestOverlapGraph(&reads_info, final_bog_name2);
				reads_info.LeftOverlaps = reads_info.LeftBestOverlapsTemp;
				reads_info.RightOverlaps = reads_info.RightBestOverlapsTemp;
				ConstructCleanedBestOverlapGraph(&reads_info, &contigs_info);
				*/
				if (OutputGraph)
				{
					string graph_name = "OverlapGraph_GapClosed";
					string bog_name = "BestOverlapGraph_GapClosed";
					OutputOverlapGraph(&reads_info, graph_name);
					OutputBestOverlapGraph(&reads_info, bog_name);
					OutputReadsInfo(&reads_info);


					LoadOverlapGraph(&reads_info, graph_name);
					LoadBestOverlapGraph(&reads_info, bog_name);
					LoadReadsInfo(&reads_info);
				}

			}


			if (0)
			if (AdaptiveTh > 0.0)
			{
				ofstream o_selected_reads("NonContainedReadNames2.txt");
				int numReads = reads_info.LeftBestOverlapsTemp.size();
				for (int i = 1; i < numReads; ++i)
				{
					if ((reads_info.LeftBestOverlapsTemp[i].size() != 0 || reads_info.RightBestOverlapsTemp[i].size() != 0))
					{
						o_selected_reads << reads_info.tag_vec[i] << endl;
					}
				}
			}

			

		}
		else
		{


			/*
			string graph_name = "OverlapGraphBeforeGapClosing";
			string bog_name = "BestOverlapGraphBeforeGapClosing";
			LoadOverlapGraph(&reads_info, graph_name);
			LoadBestOverlapGraph(&reads_info, bog_name);
			LoadReadsInfo(&reads_info);

			CloseGapsWithAllReads(&reads_info, &contigs_info, &align_info);
			*/


			string graph_name = "OverlapGraph_GapClosed";
			string bog_name = "BestOverlapGraph_GapClosed";
			LoadOverlapGraph(&reads_info, graph_name);
			LoadBestOverlapGraph(&reads_info, bog_name);
			LoadReadsInfo(&reads_info);
		}
		

		if (contigs_info.CallConsensus || reads_info.AdaptiveTh>0.0)
		{
			cout << "Loading contigs." << endl;

			ctgs_in.close();
			ctgs_in.open(ctg_filename.c_str());
			int tot_ctgs = 0;
			//contigs_info.contig_tag.clear();
			contigs_info.contigs_str.clear();
			str.clear();
			contigs_info.contigs_str.push_back(str);

			bool read_success = 1;
			string ctg_tag, ctg_tag_n;


			while (read_success)
			{

				read_success = get_a_fasta_read(ctgs_in, ctg_tag, contig_str, ctg_tag_n);
				//contigs_info.contig_tag.push_back(ctg_tag);

				while (contig_str[contig_str.size() - 1] == '\n' || contig_str[contig_str.size() - 1] == '\r')
				{
					int sz = contig_str.size();
					contig_str[sz - 1] = '\0';
					contig_str.resize(sz - 1);
				}
				int sz = contig_str.size();
				contigs_info.contigs_str.push_back(contig_str);
				tot_ctgs++;

			}
			ctgs_in.close();

		}


		if (!LOAD_OVERLAP_GRAPH)
		{

			if (reads_info.AdaptiveTh > 0)
			{
				CollectConsensusInfo(&reads_info, &contigs_info, &align_info);
			}
		}
		if (contigs_info.CallConsensus&&AdaptiveTh < 0.0)
		{
			qry_filename_vec2 = NonContainedReadsFiles;
			LoadingRawReads(qry_filename_vec, qry_filename_vec2, &reads_info, &contigs_info);

		}
		else
		{
			LoadNonContainedRawReads(&reads_info);
		}

		RemoveConflictingJoins(&contigs_info, &reads_info, &align_info);
		string out_ctg_filename = "DBG2OLC_";
		ConstructCompressedOverlapGraph(&contigs_info, &reads_info, &align_info, out_ctg_filename);
		if (AdaptiveTh > 0.0)
		{
			ConstructBackbone(&reads_info, &contigs_info, &align_info);
		}
		cout << "Assembly finished." << endl;
		

	}
	
	return 0;
}
