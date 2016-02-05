#ifndef __GRAPH_SEARCH_H
#define __GRAPH_SEARCH_H
#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <sstream>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
#include "time.h"
#include "BasicDataStructure.h"
#include "GraphConstruction.h"
#include "BuildContigs.h"

using namespace std;

struct BFS_path_info_LongRead
{
	int cov;
	int depth;
	int futile_depth;
	int len;
	int RepeatVisits;
	int last_read;
	int last_used_source;
	int last_idx;
	vector<int> path;

};

struct BFS_path_info_ShortRead
{
	int cov;
	int depth;
	int len;
	int RepeatVisits;
	int last_read;
	

};
struct BFS_path_info_Contig
{
	int cov;
	int depth;
	int futile_depth;
	int len;
	int RepeatVisits;
	int last_ctg;
	int last_used_source;
	vector<int> path;

};

struct LongRead_corrected
{
	int PathsCnt;
	uint64_t idx;
	int PathsLength;
	bool PathFound;
	uint16_t depth;
	vector< vector<int> > Paths;
};

struct search_info
{
	bool Right;
	int max_depth, max_futile_depth, max_dist,max_stack;

};
struct NodeCorr_info
{
	map<int, vector<int> > SortedShortReads;
	map<int, int> ShortReads_Score, ShortReads_MinPos, ShortReads_MaxPos;
	map<int, vector<int> > ShortReads_CoverPosition, max_aligned_pos;
	map<int, bool> ShortRead_Used,Contig_Used,Position_Covered;
	vector<Coord_CTG_Cov> RawLongRead,Aligned_Read;
	vector<int> optimal_path,path_score;
	vector< vector<int> > ShortReadsBySource,AllPaths,SelectedPaths;
	map<int, vector< vector<int> > > aligned_ref_vec,aligned_qry_vec,aligned_path;
	int NodeIdx;
	vector<int> gap_idx;
	
};

void LoadContigGraph(contigs_info *contigs_info, string filename)
{
	ifstream cg_in;
	cg_in.open(filename.c_str());

	cg_in >> contigs_info->total_contigs;
	int contig_idx;
	contigs_info->contig_adjacency_left.clear();
	contigs_info->contig_adjacency_left.resize(contigs_info->total_contigs + 1);
	contigs_info->contig_adjacency_right.clear();
	contigs_info->contig_adjacency_right.resize(contigs_info->total_contigs + 1);

	for (int i = 1; i <= contigs_info->total_contigs; ++i)
	{
		cg_in >> contig_idx;
		int left_branches, right_branches;
		cg_in >> left_branches;


		for (int b = 0; b<left_branches; ++b)
		{
			int adj_contig, overlap_sz;
			string tmp_str;
			cg_in >> adj_contig >> overlap_sz >> tmp_str;

			contigs_info->contig_adjacency_left[i][adj_contig].dist_sum = overlap_sz;
			contigs_info->contig_adjacency_left[i][adj_contig].cov = 1;
			if (adj_contig>0)
			{
				contigs_info->contig_adjacency_right[adj_contig][i].dist_sum = overlap_sz;
				contigs_info->contig_adjacency_right[adj_contig][i].cov = 1;
			}
			else
			{
				contigs_info->contig_adjacency_left[-adj_contig][-i].dist_sum = overlap_sz;
				contigs_info->contig_adjacency_left[-adj_contig][-i].cov = 1;
			}
		}
		cg_in >> right_branches;
		for (int b = 0; b<right_branches; ++b)
		{
			int adj_contig, overlap_sz;
			string tmp_str;
			cg_in >> adj_contig >> overlap_sz >> tmp_str;

			contigs_info->contig_adjacency_right[i][adj_contig].dist_sum = overlap_sz;
			contigs_info->contig_adjacency_right[i][adj_contig].cov = 1;

			if (adj_contig>0)
			{
				contigs_info->contig_adjacency_left[adj_contig][i].dist_sum = overlap_sz;
				contigs_info->contig_adjacency_left[adj_contig][i].cov = 1;
			}
			else
			{
				contigs_info->contig_adjacency_right[-adj_contig][-i].dist_sum = overlap_sz;
				contigs_info->contig_adjacency_right[-adj_contig][-i].cov = 1;
			}
		}



	}


}


void OutputContigGraph(contigs_info *contigs_info)
{
	ofstream o_contig_graph("ContigGraph.dot");
	o_contig_graph << "strict graph G {" << endl;
	for (int i = 1; i < contigs_info->contig_adjacency_left.size(); ++i)
	{
		map<int, struct adjacent_contig_info>::iterator tmp_it;
		for (tmp_it = contigs_info->contig_adjacency_left[i].begin(); tmp_it != contigs_info->contig_adjacency_left[i].end(); ++tmp_it)
		{
			//if (i<10000 && abs(tmp_it->first)<10000)
			o_contig_graph << abs(i) << " -- " << abs(tmp_it->first) << ";" << endl;
		}
		
	}
	for (int i = 1; i < contigs_info->contig_adjacency_right.size(); ++i)
	{
		map<int, struct adjacent_contig_info>::iterator tmp_it;
		for (tmp_it = contigs_info->contig_adjacency_right[i].begin(); tmp_it != contigs_info->contig_adjacency_right[i].end(); ++tmp_it)
		{

			//if (i<10000 && abs(tmp_it->first)<10000)
			o_contig_graph << abs(i) << " -- " << abs(tmp_it->first) << ";" << endl;
		}

	}

	o_contig_graph << "}" << endl;
}

void LoadCompressedOverlapGraph(struct contigs_info *contigs_info, reads_info *reads_info, string filename)
{
	string ContigPathGraphName = filename + "CompressedOverlapGraph.txt";
	ifstream graph_in(ContigPathGraphName.c_str());
	string NodesInfoName = filename + "Nodes_info.txt";
	ifstream nodes_info_in(NodesInfoName.c_str());
	
	reads_info->LeftNodes.clear();
	reads_info->RightNodes.clear();
	reads_info->LongReadIndexVec.clear();
	reads_info->Contig2LongRead.clear();
	reads_info->Contig2LongRead.resize(contigs_info->total_contigs+1);
	int node1, node2;
	while (graph_in >> node1 >> node2)
	{
		if (node1 < 0)
		{
			reads_info->LeftNodes[abs(node1)][node2] = 1;
		}
		else
		{
			reads_info->RightNodes[abs(node1)][node2] = 1;
		}
	}
	
	cout << "Loading file: " << NodesInfoName << endl;
	string str,tag;
	int n_reads = 0;
	vector<Coord_CTG_Cov> TempIndex;

	while (getline(nodes_info_in, str))
	{
		if (str[0] == '>' || str[0] == '@')
		{
			tag = str;
			n_reads++;

			reads_info->LongReadIndexVec.push_back(TempIndex);
			TempIndex.clear();

		}
		else
		{
			stringstream ss(str);
			int coord, contig, cov;
			string tmp1, tmp2;
			Coord_CTG_Cov Coord_CTG_Cov;
			ss >> coord >> tmp1 >> contig >> tmp2 >> cov;
			
			if (cov<reads_info->KmerCovTh)
			{
				continue;
			}
			Coord_CTG_Cov.coord = coord;
			Coord_CTG_Cov.contig_no = contig;
			Coord_CTG_Cov.cov = cov;
			if ((reads_info->Contig2LongRead[abs(contig)].size() == 0) || (reads_info->Contig2LongRead[abs(contig)].back() != n_reads))
			{
				if (contig>0)
				{
					reads_info->Contig2LongRead[abs(contig)].push_back(n_reads);
				}
				if (contig < 0)
				{
					reads_info->Contig2LongRead[abs(contig)].push_back(-n_reads);
				}
				
			}

			TempIndex.push_back(Coord_CTG_Cov);
		}
		//cout << n_reads << endl;
	}
	
	reads_info->LongReadIndexVec.push_back(TempIndex);
	
}


void OutputBackbone(struct contigs_info *contigs_info, reads_info *reads_info, string backbone_name)
{
	ofstream o_backbone(backbone_name.c_str());
	//ofstream o_backbone_len("backbone_len.txt");
	//string consensus_name2 = "test_" + consensus_name;
	//ofstream o_backbone_2(consensus_name2.c_str());
	//string consensus_name3 = "test2_" + consensus_name;
	//ofstream o_backbone_3(consensus_name3.c_str());
	//consensus_info consensus_info0;
	//consensus_info0.compare_len = 100;
	string ctgs_in_bb = "Contigs_in_"+backbone_name;
	ofstream o_backbone_ctgs(ctgs_in_bb.c_str());

	vector<bool> used_vec;
	used_vec.resize(contigs_info->contigs_str.size() );

	for (int i = 1; i < reads_info->LongReadIndexVec.size(); ++i)
	{
		
		string backbone_str;
		vector<int> ctg_vec;
		o_backbone_ctgs << "#Backbone_" << i << endl;
		for (int j = 0; j < reads_info->LongReadIndexVec[i].size();++j)
		{
			int ctg_idx = reads_info->LongReadIndexVec[i][j].contig_no;
			ctg_vec.push_back(ctg_idx);
			used_vec[abs(ctg_idx)] = 1;
			string contig_str = contigs_info->contigs_str[abs(ctg_idx)];
			if (ctg_idx < 0)
			{
				reverse(contig_str.begin(), contig_str.end());
				complement_str(contig_str);
			}
			if (contig_str.size() < 10)
			{
				for (int n = 0; n<contig_str.size(); ++n)
				{
					contig_str[n] = 'N';
				}
			}
			o_backbone_ctgs << ">Contig_" << abs(ctg_idx) << endl;
			o_backbone_ctgs << contig_str << endl;
			if (j>0)
			{
				int Ns = abs(reads_info->LongReadIndexVec[i][j].coord - reads_info->LongReadIndexVec[i][j - 1].coord);
				int dist = (contigs_info->contig_sz_vt[abs(reads_info->LongReadIndexVec[i][j].contig_no)]+contigs_info->contig_sz_vt[abs(reads_info->LongReadIndexVec[i][j-1].contig_no)])/2;
				Ns -= dist;
				string N_str;
				N_str.push_back('N');
				if (Ns > 1000)
				{
					//cout << Ns << endl;
				}
				for (int n = 0; n < Ns;++n)
				{
					N_str+='N';
				}
				backbone_str += N_str;

			}
			backbone_str+=contig_str;
		}
		
		if (1)//backbone_str.size()>500)
		{
			o_backbone << ">Backbone_" << i << endl;
			o_backbone << backbone_str << endl;
		}


		








	}
	/*
	for (int i = 1; i < used_vec.size(); ++i)
	{
		if (used_vec[i] == 0)
		{
			if (contigs_info->contigs_str[i].size()>0)
			{
				o_backbone << ">Contig_" << i << endl;
				o_backbone << contigs_info->contigs_str[i] << endl;

			}
		}
	}

	*/

	/*
	used_vec.clear();
	used_vec.resize(contigs_info->contigs_str.size() );
	int ctg_id = 0;
	for (int i = 0; i < reads_info->aligned_nodes.size(); ++i)
	{
		if (reads_info->aligned_nodes[i].size() == 0)
		{
			continue;
		}
		int gap_begin = 0;
		int gap_end = reads_info->aligned_nodes[i].size();
		for (int g = 0; g < reads_info->gap_idx_vec[i].size(); ++g)
		{
			gap_end = reads_info->gap_idx_vec[i][g];
			ctg_id++;
			o_consensus_2 << ">Consensus_" << ctg_id << endl;
			string consensus_str;
			for (int j = gap_begin; j < gap_end; ++j)
			{
				string contig_str = contigs_info->contigs_str[abs(reads_info->aligned_nodes[i][j])];

				if (reads_info->aligned_nodes[i][j] < 0)
				{
					reverse(contig_str.begin(), contig_str.end());
					complement_str(contig_str);
				}
				consensus_str += contig_str;
			}
			o_consensus_2 << consensus_str << endl;
			gap_begin = gap_end;
			
		}

		if (gap_begin != reads_info->aligned_nodes[i].size())
		{
			ctg_id++;
			o_consensus_2 << ">Consensus_" << ctg_id << endl;
			string consensus_str;
			gap_end = reads_info->aligned_nodes[i].size();
			for (int j = gap_begin; j < gap_end; ++j)
			{
				string contig_str = contigs_info->contigs_str[abs(reads_info->aligned_nodes[i][j])];

				if (reads_info->aligned_nodes[i][j] < 0)
				{
					reverse(contig_str.begin(), contig_str.end());
					complement_str(contig_str);
				}
				consensus_str += contig_str;

			}

			o_consensus_2 << consensus_str << endl;
		}
		

		for (int j = 0; j < reads_info->aligned_nodes[i].size(); ++j)
		{
			used_vec[abs(reads_info->aligned_nodes[i][j])] = 1;
		}
	}

	for (int i = 1; i < used_vec.size(); ++i)
	{
		if (used_vec[i] == 0)
		{
			if (contigs_info->contigs_str[i].size()>0)
			{
				o_consensus_2 << ">Contig_" << i << endl;
				o_consensus_2 << contigs_info->contigs_str[i] << endl;

			}
		}
	}




	used_vec.clear();
	used_vec.resize(contigs_info->contigs_str.size() );
	ctg_id = 0;
	for (int i = 0; i < reads_info->aligned_nodes.size(); ++i)
	{
		int ref_pos = 0;
		ctg_id++;
		o_consensus_3 << ">Consensus_" << ctg_id << endl;
		if (reads_info->aligned_nodes[i].size() == 0)
		{
			continue;
		}

		map<int, vector<int> > ctg_coords_index;

		for (int j = 0; j < reads_info->LongReadIndexVec[i+1].size(); ++j)
		{
			ctg_coords_index[reads_info->LongReadIndexVec[i+1][j].contig_no].push_back(reads_info->LongReadIndexVec[i+1][j].coord);
		}

		int gap_begin = 0;
		int gap_end = reads_info->aligned_nodes[i].size();
		for (int g = 0; g < reads_info->gap_idx_vec[i].size(); ++g)
		{
			gap_end = reads_info->gap_idx_vec[i][g];
			
			string consensus_str;
			for (int j = gap_begin; j < gap_end; ++j)
			{
				
				string contig_str = contigs_info->contigs_str[abs(reads_info->aligned_nodes[i][j])];

				if (reads_info->aligned_nodes[i][j] < 0)
				{
					reverse(contig_str.begin(), contig_str.end());
					complement_str(contig_str);
				}
				consensus_str += contig_str;
			}
			o_consensus_3 << consensus_str;


			if (gap_end>0&&gap_end < reads_info->aligned_nodes[i].size())
			{
				
				int ctg1 = reads_info->aligned_nodes[i][gap_end - 1];
				int ctg2 = reads_info->aligned_nodes[i][gap_end];
				if (ctg_coords_index[ctg1].size() == 1 && ctg_coords_index[ctg2].size() == 1)
				{
					int coord1 = ctg_coords_index[ctg1][0];
					int coord2 = ctg_coords_index[ctg2][0];
					int Ns = abs(coord1 - coord2);
					int dist = (contigs_info->contig_sz_vt[abs(ctg1)] + contigs_info->contig_sz_vt[abs(ctg2)]) / 2;
					Ns -= dist;
					string N_str;
					if (Ns > 1000)
					{
						//cout << Ns << endl;
					}
					for (int n = 0; n < Ns; ++n)
					{
						N_str += 'N';
					}
					o_consensus_3 << N_str;
				}
				else
				{
					o_consensus_3 << "N";
				}

				
			}



			gap_begin = gap_end;

		}

		if (gap_begin != reads_info->aligned_nodes[i].size())
		{
			
			string consensus_str;
			gap_end = reads_info->aligned_nodes[i].size();
			for (int j = gap_begin; j < gap_end; ++j)
			{
				string contig_str = contigs_info->contigs_str[abs(reads_info->aligned_nodes[i][j])];

				if (reads_info->aligned_nodes[i][j] < 0)
				{
					reverse(contig_str.begin(), contig_str.end());
					complement_str(contig_str);
				}
				consensus_str += contig_str;

			}

			o_consensus_3 << consensus_str << endl;
		}


		for (int j = 0; j < reads_info->aligned_nodes[i].size(); ++j)
		{
			used_vec[abs(reads_info->aligned_nodes[i][j])] = 1;
		}
	}

	for (int i = 1; i < used_vec.size(); ++i)
	{
		if (used_vec[i] == 0)
		{
			if (contigs_info->contigs_str[i].size()>0)
			{
				o_consensus_3 << ">Contig_" << i << endl;
				o_consensus_3 << contigs_info->contigs_str[i] << endl;

			}
		}
	}
	*/
}


#endif