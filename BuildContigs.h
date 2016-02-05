#ifndef __BUILD_CONTIGS_H
#define __BUILD_CONTIGS_H


#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
#include "time.h"

#include "BasicDataStructure.h"
#include "GraphConstruction.h"

using namespace std;


struct consensus_info
{
	vector<int> ReadList;
	string consensus;
	int compare_len;
	int consensus_idx;
	double mm_penalty;
};




void extend_alignment(consensus_info *consensus_info, string &append_str)
{
	int consensus_sz = consensus_info->consensus.size();
	int beg = consensus_sz - consensus_info->compare_len;
	if (beg < 0)
	{
		beg = 0;
	}

	string seq1 = consensus_info->consensus.substr(beg, consensus_info->compare_len);
	string seq2 = append_str.substr(0, consensus_info->compare_len);
	int seq1_sz = seq1.size(), seq2_sz = seq2.size();
	//cout << "1";
	//cout << seq1 << endl;
	int max_score = 0, offset = 0;


	for (int j = 0; j < consensus_info->compare_len; ++j)
	{
	
		
		int len = min(seq1_sz - j, seq2_sz);
		int match_cnt = 0;
		int mm_cnt = 0;
		for (int k = 0; k < len; ++k)
		{
			if (seq1[j + k] == seq2[k])
			{
				match_cnt++;
			}
			else
			{
				mm_cnt++;
			}
		}
		int score = (match_cnt - consensus_info->mm_penalty*(double)mm_cnt);
		if (score>max_score)
		{
			//cout << seq1.substr(j,seq1.size()) << endl;
			//cout << seq2 << endl;
			max_score = score;
			offset = len;
		}

	}
	int append_sz = append_str.size();
	int append_beg = offset;
	if (append_beg < 0)
	{
		return;
	}
	/*
	cout << "2";
	for (int k = 0; k < seq1.size() - append_beg; ++k)
	{
		cout << " ";
	}
	cout << append_str << endl;
	*/
	for (int k = append_beg; k < append_sz; ++k)
	{
		consensus_info->consensus.push_back(append_str[k]);
	}
}

void extend_alignment_fast(consensus_info *consensus_info, string &append_str)
{
	int consensus_sz = consensus_info->consensus.size();
	int beg = consensus_sz - consensus_info->compare_len;
	if (beg < 0)
	{
		beg = 0;
	}

	string seq1 = consensus_info->consensus.substr(beg, consensus_info->compare_len);
	string seq2 = append_str.substr(0, consensus_info->compare_len);
	int seq1_sz = seq1.size(), seq2_sz = seq2.size();
	//cout << "1";
	//cout << seq1 << endl;
	int max_score = 0, offset1 = 0,offset2=0;
	map<string, int> key_pos;
	map<int, int> offset_count;
	int key_len = 17;
	
	//method1:exact match
	if (seq1_sz > 0)
	{
		//cout << "q" ;
		for (int i = 0; i < seq1_sz - key_len; ++i)
		{
			string key = seq1.substr(i, key_len);
			if (key_pos.count(key) == 0)
			{
				key_pos[key] = i;
			}
			else
			{
				key_pos[key] = -1;//repeat
			}
		}

		for (int i = 0; i < seq2_sz - key_len; ++i)
		{
			string key = seq2.substr(i, key_len);

			if (key_pos.count(key) > 0 )
			{
				offset1 = key_pos[key] - i;
				if (offset1 > 0)
				{
					offset_count[offset1]++;
				}
			}

		}
		offset1 = -1;
		if (offset_count.size() == 1)
		{
			offset1 = offset_count.begin()->first;
			if (offset1 > 0)
			{
				offset1 = min(seq1_sz - offset1, seq2_sz);
			}
			//cout << "q";
		}
		else
		{
			/*
			cout << endl;
			for (map<int, int>::iterator it = offset_count.begin(); it != offset_count.end(); ++it)
			{
				cout << it->first << " " << it->second << endl;
			}
			*/
		}

		//cout << offset1 << ", ";
	}
	else
	{
		offset1 = 0;
	}




	//method2
	
	if (offset1 < 0)//slow
	{
		//cout << "s";
		for (int j = 0; j < consensus_info->compare_len; ++j)
		{


			int len = min(seq1_sz - j, seq2_sz);
			int match_cnt = 0;
			int mm_cnt = 0;
			for (int k = 0; k < len; ++k)
			{
				if (seq1[j + k] == seq2[k])
				{
					match_cnt++;
				}
				else
				{
					mm_cnt++;
				}
			}
			int score = (match_cnt - consensus_info->mm_penalty*(double)mm_cnt);
			if (score>max_score)
			{
				//cout << seq1.substr(j,seq1.size()) << endl;
				//cout << seq2 << endl;
				max_score = score;
				offset1= len;
			}

		}
	}
	
	int append_sz = append_str.size();
	int append_beg = offset1;
	/*
	if (offset1 != offset2)
	{
		

		cout << offset1 << endl << offset2 << endl;
		cout <<"r"<< seq1<< endl;
		append_beg = offset1;

		cout << "q";
		for (int k = 0; k < seq1.size() - append_beg; ++k)
		{
			cout << " ";
		}
		cout << append_str.substr(0,150) << endl;
		append_beg = offset2;
		cout << "q";
		for (int k = 0; k < seq1.size() - append_beg; ++k)
		{
			cout << " ";
		}
		cout << append_str.substr(0, 150) << endl;
	}
	*/
	if (append_beg < 0)
	{
		return;
	}
	/*
	cout << "2";
	for (int k = 0; k < seq1.size() - append_beg; ++k)
	{
	cout << " ";
	}
	cout << append_str << endl;
	*/
	for (int k = append_beg; k < append_sz; ++k)
	{
		consensus_info->consensus.push_back(append_str[k]);
	}
}


void SimpleConsensus(contigs_info *contigs_info, reads_info *reads_info, consensus_info *consensus_info)
{
	// consensus with no indel
	consensus_info->consensus.clear();
	ofstream o_debug;
	if (reads_info->Debug)
	{
		o_debug.open("consensus_debug.txt");
	}
	bool append_ctg = 1;
	int current_len = 0;
	int last_ctg = 0, first_ctg = 0;
	
	for (int r = 0; r < consensus_info->ReadList.size(); ++r)
	{
		
		string append_str;

		vector<Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[abs(consensus_info->ReadList[r])];
		vector<int> read_vec;

		for (int cc = 0; cc < current_read_layout.size(); ++cc)
		{
			read_vec.push_back(current_read_layout[cc].contig_no);
		}

		//o_contigs << ReadsList[0] << endl;
		if (consensus_info->ReadList[r] < 0)
		{
			reverse(current_read_layout.begin(), current_read_layout.end());
			for (int cc = 0; cc < current_read_layout.size(); ++cc)
			{
				current_read_layout[cc].contig_no = -current_read_layout[cc].contig_no;
			}

		}

		first_ctg = current_read_layout[0].contig_no;

		if ((last_ctg == first_ctg) || (r == 0))
		{
			append_ctg = 1;
		}
		else
		{
			append_ctg = 0;
		}
		last_ctg = current_read_layout[current_read_layout.size() - 1].contig_no;
		if (append_ctg)
		{

			int ctg_no = current_read_layout[0].contig_no;

			append_str = contigs_info->contigs_str[abs(ctg_no)];

			if (ctg_no < 0)
			{

				reverse_complement_str(append_str);
			}

			if (reads_info->Debug)
			{
				cout << ">ctg_" << ctg_no << endl;
				
				o_debug << ">ctg_" << ctg_no << endl;
				o_debug << append_str << endl;
			}

			extend_alignment_fast(consensus_info, append_str);

		}


		if (reads_info->Selected_Overlaps[read_vec].size() > 0)
		{
			append_str = reads_info->Selected_Overlaps[read_vec][0].raw_seq;
			if (consensus_info->ReadList[r] < 0)
			{
				reverse_complement_str(append_str);
			}
			if (reads_info->Debug)
			{
				cout << ">read_" << consensus_info->ReadList[r] << endl;

				o_debug << ">read_" << consensus_info->ReadList[r] << endl;
				o_debug << append_str << endl;
			}
			extend_alignment(consensus_info, append_str);
		}
		else
		{
			cout << "Missing read: " << consensus_info->ReadList[r] << endl;

		}

		if (r + 1 == consensus_info->ReadList.size())
		{
			append_str = contigs_info->contigs_str[abs(last_ctg)];

			if (last_ctg < 0)
			{
				reverse(append_str.begin(), append_str.end());
				complement_str(append_str);
			}

			if (reads_info->Debug)
			{
				cout << ">ctg_" << last_ctg << endl;
				o_debug << ">ctg_" << last_ctg << endl;
				o_debug << append_str << endl << endl;
			}
			extend_alignment(consensus_info, append_str);
		}

	}

}

void LoadRefinedContigGraph(contigs_info *contigs_info, string filename)
{
	ifstream cg_in;
	cg_in.open(filename.c_str());

	//cg_in >> contigs_info->total_contigs;
	int contig_idx;
	contigs_info->contig_adjacency_left.clear();
	contigs_info->contig_adjacency_left.resize(contigs_info->total_contigs + 1);
	contigs_info->contig_adjacency_right.clear();
	contigs_info->contig_adjacency_right.resize(contigs_info->total_contigs + 1);

	string str;
	while (getline(cg_in, str))
	{
		stringstream ss(str);
		int contig1, contig2, dist, cov;
		string bridge;

		ss >> contig1 >> contig2 >> dist >> cov >> bridge;
		if (cov > 100)
		{
			cov = 100;
		}
		if (contig1 > 0)
		{
			if (contigs_info->contig_adjacency_right[contig1][contig2].cov<cov)
			{
				contigs_info->contig_adjacency_right[contig1][contig2].dist_sum = dist*cov;
				contigs_info->contig_adjacency_right[contig1][contig2].cov = cov;
				contigs_info->contig_adjacency_right[contig1][contig2].bridge = bridge;
			}


			if (contig2>0)
			{
				if (contigs_info->contig_adjacency_left[contig2][contig1].cov < cov)
				{
					contigs_info->contig_adjacency_left[contig2][contig1].dist_sum = dist*cov;
					contigs_info->contig_adjacency_left[contig2][contig1].cov = cov;
					contigs_info->contig_adjacency_left[contig2][contig1].bridge = bridge;
				}
			}
			else
			{
				if (contigs_info->contig_adjacency_right[-contig2][-contig1].cov < cov)
				{
					contigs_info->contig_adjacency_right[-contig2][-contig1].dist_sum = dist*cov;
					contigs_info->contig_adjacency_right[-contig2][-contig1].cov = cov;

					if (bridge.size()>0)
					{
						reverse(bridge.begin(), bridge.end());
						complement_str(bridge);
					}
					contigs_info->contig_adjacency_right[-contig2][-contig1].bridge = bridge;

				}
			}
		}
		else
		{
			if (contigs_info->contig_adjacency_left[-contig1][contig2].cov<cov)
			{
				contigs_info->contig_adjacency_left[-contig1][contig2].dist_sum = dist*cov;
				contigs_info->contig_adjacency_left[-contig1][contig2].cov = cov;
				contigs_info->contig_adjacency_left[-contig1][contig2].bridge = bridge;
			}

			if (contig2>0)
			{
				if (contigs_info->contig_adjacency_right[contig2][-contig1].cov < cov)
				{
					contigs_info->contig_adjacency_right[contig2][-contig1].dist_sum = dist*cov;
					contigs_info->contig_adjacency_right[contig2][-contig1].cov = cov;
					contigs_info->contig_adjacency_right[contig2][-contig1].bridge = bridge;
				}

			}
			else
			{
				if (contigs_info->contig_adjacency_left[-contig2][contig1].cov < cov)
				{
					contigs_info->contig_adjacency_left[-contig2][contig1].dist_sum = dist*cov;
					contigs_info->contig_adjacency_left[-contig2][contig1].cov = cov;
					if (bridge.size()>0)
					{
						reverse(bridge.begin(), bridge.end());
						complement_str(bridge);
					}
					contigs_info->contig_adjacency_left[-contig2][contig1].bridge = bridge;
				}

			}
		}
	}


}

void OutputFastaContigs(contigs_info *contigs_info, string filename)
{
	ifstream in_contig_paths(filename.c_str());
	bool read_success = 1;
	string tag, n_tag;
	int Len = 0;
	int n_gaps = 0;
	vector<Coord_CTG_Cov> TempIndex;
	ofstream out_fasta_contigs("ContigsByDBG2OLC.txt");
	int n_contigs = 0;
	while (read_success)
	{
		read_success = get_a_contig_path(in_contig_paths, tag, TempIndex, Len, 0, n_tag);
		n_contigs++;
		string joint_contig_str;
		for (int i = 0; i < TempIndex.size(); ++i)
		{
			int current_contig = TempIndex[i].contig_no;
			string contig_str = contigs_info->contigs_str[abs(current_contig)];

			if (current_contig<0)
			{
				reverse(contig_str.begin(), contig_str.end());
				complement_str(contig_str);
			}
			if (i < TempIndex.size() - 1)
			{
				int next_contig = TempIndex[i + 1].contig_no;
				int dist = 0;
				string bridge;
				if (current_contig < 0)
				{
					if (contigs_info->contig_adjacency_left[abs(current_contig)].count(-next_contig))
					{
						dist = contigs_info->contig_adjacency_left[abs(current_contig)][-next_contig].dist_sum / contigs_info->contig_adjacency_left[abs(current_contig)][-next_contig].cov;
						bridge = contigs_info->contig_adjacency_left[abs(current_contig)][-next_contig].bridge;
						reverse_complement_str(bridge);
					}
					else
					{
						n_gaps++;

					}


				}
				else
				{
					if (contigs_info->contig_adjacency_right[abs(current_contig)].count(next_contig))
					{
						int cov = contigs_info->contig_adjacency_right[abs(current_contig)][next_contig].cov;
						int dist_sum = contigs_info->contig_adjacency_right[abs(current_contig)][next_contig].dist_sum;
						dist = contigs_info->contig_adjacency_right[abs(current_contig)][next_contig].dist_sum / contigs_info->contig_adjacency_right[abs(current_contig)][next_contig].cov;
						bridge = contigs_info->contig_adjacency_right[abs(current_contig)][next_contig].bridge;
					}
					else
					{
						n_gaps++;
					}
				}

				if (dist < 0)
				{
					int contig_sz = contig_str.size();
					contig_sz += dist;
					if (contig_sz < 0)
					{
						contig_sz = 0;
					}
					contig_str = contig_str.substr(0, contig_sz);
				}
				else
				{
					contig_str += bridge;
				}

			}
			joint_contig_str += contig_str;
		}
		out_fasta_contigs << ">DBG2OLC_" << n_contigs << endl;
		out_fasta_contigs << joint_contig_str << endl;
	}
	cout << n_gaps << " unexpected gaps." << endl;

}




void ConstructCompressedOverlapGraph(contigs_info *contigs_info, reads_info *reads_info, align_info *align_info, string filename)
{
	string ContigName = filename + "Nodes.txt";
	ofstream o_contigs(ContigName.c_str());
	string ContigPathsName = filename + "ContigPaths.txt";
	ofstream o_ContigPaths(ContigPathsName.c_str());
	string SuperContigsName = filename + "Nodes_info.txt";
	ofstream o_supercontig(SuperContigsName.c_str());
	string ContigLenName = filename + "NodeSizeEst.txt";
	ofstream o_contigs_len(ContigLenName.c_str());
	string Consensus_info = filename + "Consensus_info.txt";
	ofstream o_consensus_info;
	if (reads_info->AdaptiveTh > 0.0)
	{
		o_consensus_info.open(Consensus_info.c_str());
	}

	string GraphName = filename + "CompressedOverlapGraph.dot";
	ofstream o_graph(GraphName.c_str());
	GraphName = filename + "CompressedOverlapGraph.txt";
	ofstream o_graph2(GraphName.c_str());
	ofstream o_consensus;
	ofstream o_layout;
	ofstream o_layout2;

	ofstream debug_warnings;
	ofstream ext_log;

	if (reads_info->AdaptiveTh > 0.0)
	{
		o_layout.open("layout_reads.txt");
		o_layout2.open("layout_compressed_reads.txt");

	}

	if (reads_info->Debug)
	{
		debug_warnings.open("layout_warnings.txt");
	}

	if (contigs_info->CallConsensus&&reads_info->AdaptiveTh<0.0)
	{
		string consensus_name = filename + "Consensus.fasta";
		o_consensus.open(consensus_name.c_str());
	}

	//vector<int> LeftOverlapsBest, RightOverlapsBest;
	vector< map<int, int> > LeftBestOverlapsTemp, RightBestOverlapsTemp;

	map<int, int> read2node;
	vector<vector<int> > NodeInReads;
	vector<bool> used, contigs_used_vt;
	int total_ctgs = 0;
	align_matrices align_matrices;
	int frag_sum = 0;
	int offset_sum = 0;
	int ctg_cnt = 0;

	LeftBestOverlapsTemp = reads_info->LeftOverlaps;
	RightBestOverlapsTemp = reads_info->RightOverlaps;
	//either is ok
	//LeftBestOverlapsTemp = reads_info->LeftBestOverlapsTemp;
	//RightBestOverlapsTemp = reads_info->RightBestOverlapsTemp;
	contigs_used_vt.resize(contigs_info->total_contigs + 1);
	int numReads = reads_info->LenVec.size();
	used.resize(numReads);
	consensus_info consensus_info;
	consensus_info.consensus.reserve(1000000);
	int n_layout = 0;
	for (int i = 1; i < numReads; ++i)
	{
		
		if (used[i] == 0 && (LeftBestOverlapsTemp[i].size() != 0 || RightBestOverlapsTemp[i].size() != 0))
		{
			vector<int> ReadsOnTheRight, ReadsOnTheLeft, ReadsList;

			used[i] = 1;
			bool Right = 0;
			total_ctgs++;

			ReadsList.clear();

			for (int round = 1; round <= 2; ++round)
			{

				int current_read = i;
				int next_read = 0;
				if (round == 1)
				{
					Right = 0;
				}
				else
				{
					Right = 1;
				}

				while (1)
				{
					if (Right)
					{
						if (RightBestOverlapsTemp[current_read].size() != 1)
						{
							break;
						}
						else
						{
							next_read = RightBestOverlapsTemp[current_read].begin()->first;
							if (used[abs(next_read)] == 1)
							{
								break;
							}

							if (next_read< 0)
							{
								Right = 0;
								if (RightBestOverlapsTemp[abs(next_read)].size() > 1)
								{
									break;
								}
							}
							else
							{
								Right = 1;
								if (LeftBestOverlapsTemp[abs(next_read)].size() > 1)
								{
									break;
								}
							}
							ReadsList.push_back(next_read);
							current_read = abs(next_read);

						}
					}
					else
					{
						if (!Right)
						{
							if (LeftBestOverlapsTemp[current_read].size() != 1)
							{
								break;
							}
							else
							{
								next_read = -LeftBestOverlapsTemp[current_read].begin()->first;
								if (used[abs(next_read)] == 1)
								{
									break;
								}

								if (next_read< 0)
								{
									Right = 0;
									if (RightBestOverlapsTemp[abs(next_read)].size() > 1)
									{
										break;
									}
								}
								else
								{
									Right = 1;
									if (LeftBestOverlapsTemp[abs(next_read)].size() > 1)
									{
										break;
									}
								}
								ReadsList.push_back(next_read);
								current_read = abs(next_read);
							}
						}
					}


					used[current_read] = 1;

				}

				if (round == 1)
				{
					reverse(ReadsList.begin(), ReadsList.end());
					for (int r = 0; r < ReadsList.size(); ++r)
					{
						ReadsList[r] = -ReadsList[r];

					}
					ReadsList.push_back(i);
				}
			}

			

			
			bool skip_contig = 1;
			for (int r = 0; r < ReadsList.size(); ++r)
			{
				if (reads_info->chimeric[abs(ReadsList[r])] == 0)
				{
					skip_contig = 0;
					break;
				}
			}
			if (ReadsList.size()>=3)
			{
				skip_contig = 0;
			}
			if (skip_contig)
			{
				//cout << "skip" << endl;
				ReadsList.clear();
				total_ctgs--;
				continue;

			}


			int len_est = 0;
			o_contigs << ">" << total_ctgs << endl;
			o_supercontig << ">" << total_ctgs << endl;
			o_ContigPaths << ">" << total_ctgs << endl;
			consensus_info.consensus_idx = total_ctgs;
			string consensus;
			NodeInReads.push_back(ReadsList);

			for (int c = 0; c < ReadsList.size(); ++c)
			{
				o_contigs << ReadsList[c] << ", ";
				o_ContigPaths << ReadsList[c] << ", ";
				if (ReadsList[c]>0)
				{
					read2node[ReadsList[c]] = total_ctgs;
				}
				else
				{
					read2node[-ReadsList[c]] = -total_ctgs;
				}
			}
			o_contigs << endl;
			o_ContigPaths << endl;


			vector<Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[abs(ReadsList[0])];
			if (ReadsList[0] < 0)
			{
				reverse(current_read_layout.begin(), current_read_layout.end());
				for (int cc = 0; cc < current_read_layout.size(); ++cc)
				{
					current_read_layout[cc].contig_no = -current_read_layout[cc].contig_no;

				}

			}
			int coord = 0;
			list<int> coord_lst, contig_lst, cov_lst;
			for (int cc = 0; cc < current_read_layout.size(); ++cc)
			{

				len_est += contigs_info->contig_sz_vt[abs(current_read_layout[cc].contig_no)];
			
				if (cc>0)
				{
					coord += abs(current_read_layout[cc].coord2 - current_read_layout[cc - 1].coord2);
				}

				coord_lst.push_back(coord);
				contig_lst.push_back(current_read_layout[cc].contig_no);
				cov_lst.push_back(contigs_info->contig_sz_vt[abs(current_read_layout[cc].contig_no)]);
				contigs_used_vt[abs(current_read_layout[cc].contig_no)] = 1;

			}


			if (contigs_info->CallConsensus&&reads_info->AdaptiveTh<0.0)
			{


				consensus_info.compare_len = reads_info->MaxReadLen;
				consensus_info.mm_penalty = align_info->mm_penalty;
				consensus_info.ReadList = ReadsList;
				if (reads_info->AdaptiveTh < 0.0)
				{
					o_consensus << ">" << total_ctgs << endl;

					SimpleConsensus(contigs_info, reads_info, &consensus_info);
					o_consensus << consensus_info.consensus << endl;
				}
			
			}

			vector<int> offset_vt;
			vector<string> raw_reads;
			offset_vt.push_back(0);
			for (int c = 0; c < ReadsList.size() - 1; ++c)
			{
				int current_read = ReadsList[c], next_read = ReadsList[c + 1];
				vector<Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[abs(current_read)];
				vector<Coord_CTG_Cov> next_read_layout = reads_info->LongReadIndexVec[abs(next_read)];
				if (current_read < 0)
				{
					int RLen = reads_info->LenVec[abs(current_read)];
					if (RLen == 0)
					{
						RLen = current_read_layout[current_read_layout.size() - 1].coord2;
					}
					reverse(current_read_layout.begin(), current_read_layout.end());
					for (int cc = 0; cc < current_read_layout.size(); ++cc)
					{
						current_read_layout[cc].contig_no = -current_read_layout[cc].contig_no;
						current_read_layout[cc].coord2 = RLen - current_read_layout[cc].coord2;
					}

				}
				if (next_read < 0)
				{
					int RLen = reads_info->LenVec[abs(next_read)];
					if (RLen == 0)
					{
						RLen = next_read_layout[next_read_layout.size() - 1].coord2;
					}
					reverse(next_read_layout.begin(), next_read_layout.end());
					for (int cc = 0; cc < next_read_layout.size(); ++cc)
					{
						next_read_layout[cc].contig_no = -next_read_layout[cc].contig_no;
						next_read_layout[cc].coord2 = RLen - next_read_layout[cc].coord2;
					}

				}

				if (reads_info->AdaptiveTh>0.0)
				{
					align_info->CalculateOffset = 1;
					align_info->ref_len = reads_info->LenVec[abs(current_read)];
					align_info->qry_len = reads_info->LenVec[abs(next_read)];
					

				}
				align_info->fix_orientation = 1;
				align_info->force_flip = 0;
				if (align_info->SparseAlign)
				{
					sparse_semi_global_align(contigs_info, current_read_layout, next_read_layout, &align_matrices, align_info);

				}
				else
				{
					approximate_semi_global_align(contigs_info, current_read_layout, next_read_layout, &align_matrices, align_info);

				}
				
				int result = classify_alignment(align_info);
				bool debug_flag = 0;
				if (result != 1)
				{

					debug_flag = 1;
					offset_vt.push_back(1);
					debug_warnings << "Warning. Alignment error: " << result << endl;					
					debug_warnings << "ref: " << endl;
					for (int l = 0; l < current_read_layout.size(); ++l)
					{
						debug_warnings << current_read_layout[l].coord << ", ";
						debug_warnings << current_read_layout[l].contig_no << ", ";
						debug_warnings << current_read_layout[l].cov << ", ";
						debug_warnings << current_read_layout[l].coord2 << ", ";
						debug_warnings << endl;
					}
					debug_warnings << "qry: " << endl;
					for (int l = 0; l < next_read_layout.size(); ++l)
					{
						debug_warnings << next_read_layout[l].coord << ", ";
						debug_warnings << next_read_layout[l].contig_no << ", ";
						debug_warnings << next_read_layout[l].cov << ", ";
						debug_warnings << next_read_layout[l].coord2 << ", ";
						debug_warnings << endl;
					}
					debug_warnings << "Result:" << result << endl;

					
					debug_warnings << endl;
				}
				result = 1;
				if (1)
				{

					offset_vt.push_back(align_info->ahg);

					if (coord_lst.size() >= align_info->right_mismatch.size())
					{
						for (int cc = 0; cc < align_info->right_mismatch.size(); ++cc)
						{
							coord_lst.pop_back();
							contig_lst.pop_back();
							cov_lst.pop_back();


						}
					}
					else
					{
						debug_warnings << "Pop error." << endl;
						debug_warnings << "Pop " << align_info->right_mismatch.size() << " nodes from size " << coord_lst.size()<<" list." <<endl;
						//cout << "Pop error." << endl;
						//cout<< "Pop " << align_info->right_mismatch.size() << " nodes from size " << coord_lst.size() << " list." << endl;
						debug_warnings << "Mismatch list: ";
						for (int cc = 0; cc < align_info->right_mismatch.size(); ++cc)
						{
							debug_warnings << align_info->right_mismatch[cc] << " ";
						}
						list<int>::iterator lit = contig_lst.begin();
						debug_warnings << "contig list: ";
						for (int cc = 0; cc < contig_lst.size(); ++cc)
						{
							debug_warnings << *lit << " ";
							lit++;
						}
						debug_warnings << endl;
					}

					for (int cc = 0; cc < align_info->right_ext.size(); ++cc)
					{
						len_est += contigs_info->contig_sz_vt[abs(align_info->right_ext[cc])];
						coord += align_info->right_dist[cc];
						coord_lst.push_back(coord);
						contig_lst.push_back(align_info->right_ext[cc]);
						cov_lst.push_back(contigs_info->contig_sz_vt[abs(align_info->right_ext[cc])]);
						contigs_used_vt[abs(align_info->right_ext[cc])] = 1;
					}

				}

				o_ContigPaths <<"id: "<< current_read << endl;

				for (int cc = 0; cc < align_info->ref_aligned.size(); ++cc)
				{
					o_ContigPaths << align_info->ref_aligned[cc] << ", ";
				}
				o_ContigPaths << endl;
				o_ContigPaths << "id: " << next_read << endl;

				for (int cc = 0; cc < align_info->qry_aligned.size(); ++cc)
				{
					o_ContigPaths << align_info->qry_aligned[cc] << ", ";
				}
				o_ContigPaths << endl;



			}

			if (reads_info->AdaptiveTh > 0.0)
			{
				for (int c = 0; c < ReadsList.size(); ++c)
				{
					frag_sum += reads_info->LenVec[abs(ReadsList[c])];
				}
				for (int c = 0; c < offset_vt.size(); ++c)
				{

					offset_sum += offset_vt[c];
				}

			}
			// output compressed contig
			list<int>::iterator coord_it = coord_lst.begin();
			list<int>::iterator cov_it = cov_lst.begin();
			for (list<int>::iterator contig_it = contig_lst.begin(); contig_it != contig_lst.end(); ++contig_it)
			{
				o_contigs << (*contig_it) << ", ";

				o_supercontig << (*coord_it) << ", " << (*contig_it) << ", " << *cov_it << endl;
				cov_it++;
				coord_it++;
			}

			o_contigs << endl;
			o_contigs_len << len_est << endl;

			int front_hang_sz = 0, back_hang_sz = 0;

			if (reads_info->AdaptiveTh>0.0)
			{
				//get contig coverage to pick the confident ones
				map<int, int> ctg_cov;
				for (int r = 0; r < ReadsList.size(); ++r)
				{
					vector<Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[abs(ReadsList[r])];
					for (int cc = 0; cc < current_read_layout.size(); ++cc)
					{
						int ctg_id;
						if (ReadsList[r] < 0)
						{
							ctg_id = -current_read_layout[cc].contig_no;
						}
						else
						{
							ctg_id = current_read_layout[cc].contig_no;
						}

						if (ctg_cov[ctg_id] < current_read_layout[cc].cov)
						{
							ctg_cov[ctg_id] = current_read_layout[cc].cov;
						}
					}
				}

				if (reads_info->AdaptiveTh > 0.0)
				{
					o_layout << ">Backbone_" << total_ctgs << endl;
					o_layout2 << ">Backbone_" << total_ctgs << endl;
					o_layout << ReadsList.size() << endl; 
					o_layout2 << ReadsList.size() << endl;


					for (int r = 0; r < ReadsList.size(); ++r)
					{
						n_layout++;
						string raw_read_str = reads_info->selected_long_reads_seq[abs(ReadsList[r])];
						if (ReadsList[r]<0)
						{
							reverse_complement_str(raw_read_str);
						}
						

						vector<Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[abs(ReadsList[r])];
						if (ReadsList[r]<0)
						{
							int RLen = reads_info->LenVec[abs(ReadsList[r])];
							if (RLen == 0)
							{
								RLen = current_read_layout[current_read_layout.size() - 1].coord2;
							}
							reverse(current_read_layout.begin(), current_read_layout.end());
							for (int cc = 0; cc < current_read_layout.size(); ++cc)
							{
								current_read_layout[cc].contig_no = -current_read_layout[cc].contig_no;

								current_read_layout[cc].coord2 = RLen - current_read_layout[cc].coord2;

								current_read_layout[cc].coord = RLen - current_read_layout[cc].coord;
							}


						}
						o_layout << reads_info->tag_vec[abs(ReadsList[r])] << endl;
						o_layout << raw_read_str << endl;
						o_layout2 << reads_info->tag_vec[abs(ReadsList[r])] << endl;
						
						for (int l = 0; l < current_read_layout.size(); ++l)
						{
							o_layout2 << current_read_layout[l].coord << ", ";
							o_layout2 << current_read_layout[l].contig_no << ", ";
							o_layout2 << current_read_layout[l].cov << ", ";
							o_layout2 << current_read_layout[l].coord2 << ", ";
							o_layout2 << endl;
						}


						/*
						if (r > 0)
						{
							//debug_layout2 << "Length before: " << backbone_str.size() << endl;
							extend_long_read_alignment(contigs_info, reads_info, align_info, backbone_str1, raw_read_str, ReadsList[r - 1], ReadsList[r]);
							extend_long_read_alignment_corr(contigs_info, reads_info, align_info, backbone_str, raw_read_str, ReadsList[r - 1], ReadsList[r], ctg_cov);
							//extend_long_read_alignment_coords(contigs_info, reads_info, align_info, backbone_str, raw_read_str, ReadsList[r - 1], ReadsList[r]);
							//debug_layout2 << "Length after: " << backbone_str.size() << endl;
						}
						else
						{
							backbone_str = raw_read_str;
							backbone_str1 = raw_read_str;
							if (0)
							if (reads_info->Debug)
							{
								ctg_cnt++;
								ext_log << endl;
								ext_log << "asmContig_" << ctg_cnt << endl;
								ext_log << reads_info->tag_vec[abs(ReadsList[0])] << endl;

								if (ReadsList[0] > 0)
								{
									ext_log << 0;
									ext_log << " " << reads_info->LenVec[abs(ReadsList[0])] << endl;
								}
								else
								{

									ext_log << reads_info->LenVec[abs(ReadsList[0])];
									ext_log << " " << 0 << endl;
								}
							}

						}
						*/

					}
					//o_backbone << backbone_str << endl;
					//o_backbone1 << backbone_str1 << endl;
				}




				o_consensus_info << ">Backbone_" << total_ctgs << endl;
				map<int, int> ctg_used;
				for (list<int>::iterator contig_it = contig_lst.begin(); contig_it != contig_lst.end(); ++contig_it)
				{
					if ((ctg_cov.count(*contig_it) == 0) || (ctg_cov[*contig_it] < max(0.005,reads_info->AdaptiveTh)*contigs_info->contig_sz_vt[abs(abs(*contig_it))]))
					{

						continue;
					}
					if (ctg_used[*contig_it]>0)
					{
						continue;
					}
					ctg_used[*contig_it] = 1;
					string tag = contigs_info->contig_tag[abs(*contig_it)];
					if (tag[0] == '>' || tag[0] == '@')
					{
						tag = tag.substr(1, tag.size());
					}
					if (contigs_info->contig_sz_vt[abs(*contig_it)] > 30)
					{
						o_consensus_info << tag << endl;
						//contigs
					}

				}


				vector<int> consensus_read_list;
				map<int, int> used_in_backbone;
				for (int c = 0; c < ReadsList.size(); ++c)
				{
					int current_read = abs(ReadsList[c]);
					used_in_backbone[current_read]++;
					//consensus_read_list.push_back(current_read);
					//the reads used in the backbone 
					if (reads_info->Consensus_info.count(current_read))
					{
						for (int cc = 0; cc < reads_info->Consensus_info[current_read].size(); ++cc)
						{
							if (1)//used_in_backbone.count(abs(reads_info->Consensus_info[current_read][cc])) == 0)
							{
								consensus_read_list.push_back(abs(reads_info->Consensus_info[current_read][cc]));
								//only use those not used in backbone
							}
							else
							{
								cout << "skip.";
							}
							
						}
					}

				}
				sort(consensus_read_list.begin(), consensus_read_list.end());
				if (consensus_read_list.size()>0)
				{

					string tag = reads_info->tag_vec[consensus_read_list[0]];
					if (tag[0] == '>' || tag[0] == '@')
					{
						tag = tag.substr(1, tag.size());
					}
					o_consensus_info << tag << endl;
				}
				for (int c = 1; c < consensus_read_list.size(); ++c)
				{

					if (consensus_read_list[c] != consensus_read_list[c - 1])
					{
						if (used_in_backbone.count(consensus_read_list[c]))
						{
							//cout << "err" << endl;
							//continue;
						}
						string tag = reads_info->tag_vec[consensus_read_list[c]];
						if (tag[0] == '>' || tag[0] == '@')
						{
							tag = tag.substr(1, tag.size());
						}
						o_consensus_info << tag << endl;
					}
				}

			}

		}



	}
	for (int i = 1; i <= contigs_info->total_contigs; ++i)
	{
		if (contigs_used_vt[i] == 0 && contigs_info->contig_sz_vt[i]>100)
		{
			total_ctgs++;
			
			if (contigs_info->CallConsensus&&reads_info->AdaptiveTh<0.0)
			{

				o_consensus << ">" << total_ctgs << endl;

				o_consensus << contigs_info->contigs_str[i] << endl;


			}
		}
	}
	cout << "frag sum: " << frag_sum << endl;
	cout << "offset sum: " << offset_sum << endl;


	//output graph

	///
	map<int, map<int, int> > LeftNodes, RightNodes;

	for (int i = 0; i < NodeInReads.size(); ++i)
	{
		int NodeId = i + 1;



		if (NodeInReads[i].size() == 0)
		{
			continue;
		}

		//left
		int read_id = NodeInReads[i][0];

		if (read_id > 0)
		{
			for (map<int, int>::iterator it1 = LeftBestOverlapsTemp[read_id].begin(); it1 != LeftBestOverlapsTemp[read_id].end(); ++it1)
			{
				int next_read = it1->first;
				bool flip1 = 0;
				if (next_read < 0)
				{
					flip1 = 1;
				}
				bool flip2 = 0;
				int next_node = read2node[abs(next_read)];
				if (next_node == 0)
				{
					continue;
				}
				if (next_node<0)
				{
					flip2 = 1;
				}
				bool flip = flip1^flip2;
				next_node = abs(next_node);
				if (flip == 0)
				{
					LeftNodes[NodeId][next_node] = 1;
				}
				else
				{
					LeftNodes[NodeId][-next_node] = 1;
				}

			}

		}
		else
		{
			read_id = -read_id;

			for (map<int, int>::iterator it1 = RightBestOverlapsTemp[read_id].begin(); it1 != RightBestOverlapsTemp[read_id].end(); ++it1)
			{
				int next_read = -it1->first;
				bool flip1 = 0;
				if (next_read < 0)
				{
					flip1 = 1;
				}
				bool flip2 = 0;
				int next_node = read2node[abs(next_read)];
				if (next_node == 0)
				{
					continue;
				}
				if (next_node<0)
				{
					flip2 = 1;
				}
				bool flip = flip1^flip2;
				next_node = abs(next_node);
				if (flip == 0)
				{
					LeftNodes[NodeId][next_node] = 1;
				}
				else
				{
					LeftNodes[NodeId][-next_node] = 1;
				}

			}
		}

		//right
		int len = NodeInReads[i].size();
		read_id = NodeInReads[i][len - 1];
		if (read_id < 0)
		{
			read_id = -read_id;
			for (map<int, int>::iterator it1 = LeftBestOverlapsTemp[read_id].begin(); it1 != LeftBestOverlapsTemp[read_id].end(); ++it1)
			{
				int next_read = -it1->first;
				bool flip1 = 0;
				if (next_read < 0)
				{
					flip1 = 1;
				}
				bool flip2 = 0;
				int next_node = read2node[abs(next_read)];
				if (next_node == 0)
				{
					continue;
				}
				if (next_node<0)
				{
					flip2 = 1;
				}
				bool flip = flip1^flip2;
				next_node = abs(next_node);
				if (flip == 0)
				{
					RightNodes[NodeId][next_node] = 1;
				}
				else
				{
					RightNodes[NodeId][-next_node] = 1;
				}

			}

		}
		else
		{
			for (map<int, int>::iterator it1 = RightBestOverlapsTemp[read_id].begin(); it1 != RightBestOverlapsTemp[read_id].end(); ++it1)
			{
				int next_read = it1->first;
				bool flip1 = 0;
				if (next_read < 0)
				{
					flip1 = 1;
				}
				bool flip2 = 0;
				int next_node = read2node[abs(next_read)];
				if (next_node == 0)
				{
					continue;
				}
				if (next_node<0)
				{
					flip2 = 1;
				}
				bool flip = flip1^flip2;
				next_node = abs(next_node);
				if (flip == 0)
				{
					RightNodes[NodeId][next_node] = 1;
				}
				else
				{
					RightNodes[NodeId][-next_node] = 1;
				}

			}
		}
	}

	o_graph << "strict graph G {" << endl;
	map<int, map<int, int> >::iterator graph_it1;
	for (graph_it1 = LeftNodes.begin(); graph_it1 != LeftNodes.end(); ++graph_it1)
	{
		map<int, int>::iterator graph_it2;
		for (graph_it2 = graph_it1->second.begin(); graph_it2 != graph_it1->second.end(); ++graph_it2)
		{
			o_graph << abs(graph_it1->first) << " -- " << abs(graph_it2->first) << ";" << endl;
			o_graph2 << -graph_it1->first << " " << graph_it2->first << endl;
		}
	}
	for (graph_it1 = RightNodes.begin(); graph_it1 != RightNodes.end(); ++graph_it1)
	{
		map<int, int>::iterator graph_it2;
		for (graph_it2 = graph_it1->second.begin(); graph_it2 != graph_it1->second.end(); ++graph_it2)
		{
			o_graph << abs(graph_it1->first) << " -- " << abs(graph_it2->first) << ";" << endl;
			o_graph2 << graph_it1->first << " " << graph_it2->first << endl;
		}
	}
	o_graph << "}" << endl;


	
}


void RemoveConflictingJoins(contigs_info *contigs_info, reads_info *reads_info, align_info *align_info_default)
{
	//cout << "Removing Conflicting joins." << endl;
	align_matrices align_matrices;
	align_info left_align_info = *align_info_default;
	align_info right_align_info = *align_info_default;
	
	for (int i = 1; i < reads_info->LenVec.size(); ++i)
	{
		bool remove = 0;
		if ((!reads_info->CloseGaps) && (reads_info->chimeric[i]) && (reads_info->LeftOverlaps[i].size()>0 || reads_info->RightOverlaps[i].size()>0))
		{
			remove = 1;
		}
		if (reads_info->LeftOverlaps[i].size() == 1 && reads_info->RightOverlaps[i].size() == 1)
		{
			int left_read = reads_info->LeftOverlaps[i].begin()->first;
			int right_read = reads_info->RightOverlaps[i].begin()->first;

			left_align_info.fix_orientation = 1;
			right_align_info.fix_orientation = 1;
			left_align_info.force_flip = 0;
			if (left_read < 0)
			{
				left_align_info.force_flip = 1;
			}
			right_align_info.force_flip = 0;
			if (right_read < 0)
			{
				right_align_info.force_flip = 1;
			}
			vector<Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[i];
			vector<Coord_CTG_Cov> left_read_layout = reads_info->LongReadIndexVec[abs(left_read)];
			vector<Coord_CTG_Cov> right_read_layout = reads_info->LongReadIndexVec[abs(right_read)];
			if (align_info_default->SparseAlign)
			{
				sparse_semi_global_align(contigs_info, current_read_layout, left_read_layout, &align_matrices, &left_align_info);

				sparse_semi_global_align(contigs_info, current_read_layout, right_read_layout, &align_matrices, &right_align_info);

			}
			else
			{
				approximate_semi_global_align(contigs_info, current_read_layout, left_read_layout, &align_matrices, &left_align_info);

				approximate_semi_global_align(contigs_info, current_read_layout, right_read_layout, &align_matrices, &right_align_info);

			}
			
			int left_min_match = 0,right_max_match=0;
			for (int j = 0; j<=left_align_info.min_match; ++j)
			{
				if (left_align_info.ref_aligned[j] != 0)
				{
					left_min_match++;
				}
			}
			for (int j = 0; j<=right_align_info.max_match; ++j)
			{
				if (right_align_info.ref_aligned[j] != 0)
				{
					right_max_match++;
				}
			}
			if (left_min_match > right_max_match)
			{
				remove = 1;
			}
			
		}
		if (remove)
		{
			reads_info->contained[i] = 1; //bad read.
			for (map<int, int>::iterator lit = reads_info->LeftOverlaps[i].begin(); lit != reads_info->LeftOverlaps[i].end();++lit)
			{
				int left_read = lit->first;
				if (left_read > 0)
				{
					reads_info->RightOverlaps[left_read].erase(i);
				}
				else
				{
					reads_info->LeftOverlaps[-left_read].erase(-i);
				}

			}

			for (map<int, int>::iterator rit = reads_info->RightOverlaps[i].begin(); rit != reads_info->RightOverlaps[i].end(); ++rit)
			{
				int right_read = rit->first;

				if (right_read > 0)
				{
					reads_info->LeftOverlaps[right_read].erase(i);
				}
				else
				{
					reads_info->RightOverlaps[-right_read].erase(-i);
				}
			}
			reads_info->LeftOverlaps[i].clear();
			reads_info->RightOverlaps[i].clear();
		
		}
	}
	//cout << "done." << endl;

}


#endif