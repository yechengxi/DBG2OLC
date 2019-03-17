#ifndef __GRAPH_SIMPLIFICATION_H
#define __GRAPH_SIMPLIFICATION_H
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


struct BFS_reads_info
{
	int cov;
	int depth;
	int len;
	int last_read;
	//vector<int> edge;
};


bool isSimplePath_read(reads_info *reads_info, int current_read, map<int, struct BFS_reads_info > &Visited_Path, map<int, int > &stacked_nodes)
{
	//return 1;
	int last_read = abs(current_read);
	int dep = Visited_Path[abs(current_read)].depth;
	int node = last_read;
	for (int l = dep; l >= 2; --l)//l>2
	{

		node = abs(Visited_Path[abs(node)].last_read);
		if (node == 0)
		{
			return 1;
		}

		if (Visited_Path[abs(node)].depth == 1)
		{
			return 1;// it is a simple bubble;
		}


		if (stacked_nodes[node]>0 && reads_info->LeftOverlaps[abs(node)].size()>1)
		{
			return 0;
		}//backward branch

		if (stacked_nodes[node]<0 && reads_info->RightOverlaps[abs(node)].size()>1)
		{
			return 0;
		}//backward branch

		if (abs(stacked_nodes[node])>2)// it is a complicated bubble;
		{
			//return 1;
			return 0;
		}//{return 1;}// only backtrack to here so return 1.

	}
	return 1;
}

void BreakLinks_read(reads_info *reads_info, map<int, int > &stacked_nodes, int node1, int node2)
{
	
	node1 = abs(node1);
	node2 = abs(node2);
	if (stacked_nodes[node1]>0)
	{
		int temp = abs(node2);
		if (stacked_nodes[node2]<0)
		{
			temp = -temp;
		}
		if (reads_info->RightOverlaps[abs(node1)].count(temp))
		{
			reads_info->RightOverlaps[abs(node1)].erase(temp);
			if (stacked_nodes[node1]>1)
			{
				stacked_nodes[abs(node1)]--;
			}
		}

	}


	if (stacked_nodes[node1]<0)
	{
		int temp = abs(node2);
		if (stacked_nodes[node2]>0)
		{
			temp = -temp;
		}

		if (reads_info->LeftOverlaps[abs(node1)].count(temp))
		{
			reads_info->LeftOverlaps[abs(node1)].erase(temp);
			if (stacked_nodes[node1]<-1)
			{
				stacked_nodes[abs(node1)]++;
			}
		}


	}


	if (stacked_nodes[node2]<0)
	{
		int temp = abs(node1);
		if (stacked_nodes[node1]>0)
		{
			temp = -temp;
		}
		if (reads_info->RightOverlaps[abs(node2)].count(temp))
		{
			reads_info->RightOverlaps[abs(node2)].erase(temp);

		}
	}

	if (stacked_nodes[node2]>0)
	{
		int temp = abs(node1);
		if (stacked_nodes[node1]<0)
		{
			temp = -temp;
		}
		if (reads_info->LeftOverlaps[abs(node2)].count(temp))
		{
			reads_info->LeftOverlaps[abs(node2)].erase(temp);

		}
	}
	reads_info->n_deleted_edges++;
}

void BacktrackBubbleRemoval_read(reads_info *reads_info, int last_read, int beg_read, map<int, struct BFS_reads_info > & Visited_Path, map<int, int > &stacked_nodes)
{
	beg_read = abs(beg_read);
	last_read = abs(last_read);
	int current_read = last_read;
	int dep = Visited_Path[last_read].depth;
	for (int l = dep; l>1; --l)
	{
		int previous_read = current_read;

		current_read = Visited_Path[current_read].last_read;
		if (current_read == 0)
		{
			return;
		}

		if (stacked_nodes[current_read] >= 1 || stacked_nodes[current_read] <= -1)
		{

			if (abs(stacked_nodes[current_read]) == 1)
			{
				//already removed somehow
				break;
			}

			if (abs(stacked_nodes[current_read])>2)
			{
				BreakLinks_read(reads_info, stacked_nodes, current_read, previous_read);
				break;
			}
			//else

			BreakLinks_read(reads_info, stacked_nodes, current_read, previous_read);

			if (Visited_Path[current_read].last_read == 0)
			{
				break;
			}

			int free_read = current_read;

			if (beg_read == free_read)
			{
				break;
			}

			if (free_read != last_read)
			{

				stacked_nodes[free_read] = stacked_nodes[free_read] / abs(stacked_nodes[free_read]);
				reads_info->contained[abs(free_read)] = 1;
				//freebkt->kmer_info.removed=1;			
			}
		}
	}
}

bool BackCheckLoop_read(int current_read, map<int, struct BFS_reads_info > &Visited_Path)
{


	int last_read = abs(current_read);
	int dep = Visited_Path[abs(current_read)].depth;
	int node = last_read;
	for (int l = dep; l >= 1; --l)//l>2
	{

		node = abs(Visited_Path[abs(node)].last_read);
		if (node == 0)
		{
			return 0;
		}
		if (node == last_read)
		{
			return 1;
		}
	}
	return 0;

}

void BFSearchBubbleRemoval_read(reads_info *reads_overlap_info, int beg_read, int max_depth, int max_dist, bool RemoveTipsOnly)
{


	if (abs(beg_read) == 16140 || abs(beg_read) == 14823 || abs(beg_read) == 13905)
	{
		cout << "";
	}
	map<int, struct BFS_reads_info > Visited_Path;
	map<int, int > stacked_nodes;
	int max_stack = 300;
	int DepthTh = max_depth;
	int LenTh = 30;
	bool RIGHT = 0;
	map<int, list<int> > dist_reads;//neighborset
	dist_reads[0].push_back(beg_read);
	int NBs = 1;
	int dist_searched = 0;

	int new_node = abs(beg_read);
	if (beg_read>0)
	{
		stacked_nodes[abs(beg_read)] = 1;
	}
	else
	{
		stacked_nodes[abs(beg_read)] = -1;
	}

	Visited_Path[new_node].cov = 0;
	Visited_Path[new_node].depth = 1;
	Visited_Path[new_node].last_read = 0;

	map<int, list<int> >::iterator NB_it = dist_reads.begin();

	while (1)
	{
		NB_it = dist_reads.begin();
		if (NB_it == dist_reads.end())
		{
			break;
		}
		if (NB_it->second.size() == 0)
		{
			dist_reads.erase(NB_it->first); continue;
		}

		if (NBs>max_stack)
		{
			break;
		}
		new_node = NB_it->second.front();

		NB_it->second.pop_front();
		NBs--;
		if (NB_it->second.size() == 0)
		{
			dist_reads.erase(NB_it->first);
		}
		if (new_node>0)
		{
			RIGHT = 1;
		}
		else
		{
			RIGHT = 0;
		}

		new_node = abs(new_node);

		if (Visited_Path[new_node].depth>DepthTh || Visited_Path[new_node].len>LenTh)
		{
			continue;
		}


		if (RIGHT)
		{
			int rb = reads_overlap_info->RightOverlaps[new_node].size();
			if (stacked_nodes[new_node] == 1 && rb>0)
			{
				stacked_nodes[new_node] = 1 + rb;
			}
			
			if (rb == 0  )// RemoveTipsOnly
			{
				stacked_nodes[new_node] = 2;

				if (Visited_Path[abs(new_node)].depth>((max_depth+1)/2))// let's define a tip to be non-extendable (2)
				{
					continue;
				}
				if (abs(new_node) == abs(beg_read))
				{
					continue;
				}
				//tip end reached so backtrack to the branching position.
				if (!isSimplePath_read(reads_overlap_info, new_node, Visited_Path, stacked_nodes))
				{
					continue;
				}
				//reads_overlap_info->contained[abs(new_node)] = 1;
				BacktrackBubbleRemoval_read(reads_overlap_info, new_node, beg_read, Visited_Path, stacked_nodes);
				if (beg_read>0)
				{
					reads_overlap_info->used_vt_right[abs(beg_read)] = 1;
				}
				else
				{
					reads_overlap_info->used_vt_left[abs(beg_read)] = 1;
				}
				stacked_nodes[new_node] = 1;
				///////////////latest modification
				if (reads_overlap_info->LeftOverlaps[new_node].size())
				{
					//dist_reads[0].push_back(-new_node);
				}
				break;

			}

			map<int, int >::iterator tmp_it, tmp_it_n;

			for (tmp_it = reads_overlap_info->RightOverlaps[new_node].begin(); tmp_it != reads_overlap_info->RightOverlaps[new_node].end();)
			{
				tmp_it_n = tmp_it;
				tmp_it_n++;
				int next_read = tmp_it->first;
				// not in stack
				if (stacked_nodes[abs(next_read)] == 0)
				{
					//					Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
					//Visited_Path[abs(next_read)].cov=(int)(Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]);
					Visited_Path[abs(next_read)].cov = Visited_Path[abs(new_node)].cov+tmp_it->second; //(int)(Visited_Path[abs(new_node)].cov);// +reads_overlap_info->right_overlap_size[new_node][next_read]);
					Visited_Path[abs(next_read)].depth = (Visited_Path[abs(new_node)].depth + 1);
					//int cum_len=(int)(Visited_Path[abs(next_read)].len+1);
					int cum_len = Visited_Path[abs(next_read)].depth;
					//Visited_Path[abs(next_read)].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);
					Visited_Path[abs(next_read)].last_read = new_node;

					if (next_read<0)
					{
						stacked_nodes[abs(next_read)] = -1;
					}
					else
					{
						stacked_nodes[abs(next_read)] = 1;
					}
					dist_reads[cum_len].push_back(next_read);
					NBs++;

				}
				else
				{
					if (!RemoveTipsOnly)
					if ((stacked_nodes[abs(next_read)]>0 && next_read>0) || (stacked_nodes[abs(next_read)]<0 && next_read<0))
					{

						//if((Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]<=Visited_Path[abs(next_read)].cov))//||(BackCheckLoop(*ptr,new_node,Visited_Path)==1)//loop
						if ((Visited_Path[abs(new_node)].cov + tmp_it->second <= Visited_Path[abs(next_read)].cov) || (BackCheckLoop_read(new_node, Visited_Path) == 1))//||(BackCheckLoop(*ptr,new_node,Visited_Path)==1)//loop
						{
							//backtrack if the same direction is found

							if (!isSimplePath_read(reads_overlap_info, new_node, Visited_Path, stacked_nodes))
							{
								tmp_it = tmp_it_n;
								continue;
							}

							//backtrack the current path, common operation in this search
							if (stacked_nodes[new_node]>2)
							{
								BreakLinks_read(reads_overlap_info, stacked_nodes, new_node, next_read);
								tmp_it = tmp_it_n;
								break;
							}
							else
							{
								if (stacked_nodes[new_node]<-2)
								{
									BreakLinks_read(reads_overlap_info, stacked_nodes, new_node, next_read);

									tmp_it = tmp_it_n;
									break;
								}
								else
								{
									int free_node = (new_node);

									if (free_node == beg_read)
									{
										BreakLinks_read(reads_overlap_info, stacked_nodes, beg_read, new_node);
									}

									//reads_overlap_info->contained[abs(new_node)] = 1;
									//free the node and edge.
									BreakLinks_read(reads_overlap_info, stacked_nodes, new_node, next_read);
									BacktrackBubbleRemoval_read(reads_overlap_info, new_node, beg_read, Visited_Path, stacked_nodes);


									tmp_it = tmp_it_n;
									break;
									//
								}
							}

						}
						else
						{
							//backtrack the original path, rare operation in this search, can lead to errors

							if (!isSimplePath_read(reads_overlap_info, next_read, Visited_Path, stacked_nodes))
							{
								tmp_it = tmp_it_n;
								continue;
							}

							BacktrackBubbleRemoval_read(reads_overlap_info, next_read, beg_read, Visited_Path, stacked_nodes);

							//marginal 
							//Visited_Path[abs(next_read)].cov=(int)(Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]);
							Visited_Path[abs(next_read)].cov = (int)(Visited_Path[abs(new_node)].cov) +tmp_it->second;

							Visited_Path[abs(next_read)].depth = Visited_Path[abs(new_node)].depth + 1;
							Visited_Path[abs(next_read)].len = Visited_Path[abs(new_node)].depth;
							//marginal 

							Visited_Path[abs(next_read)].last_read = abs(new_node);
							tmp_it = tmp_it_n;
							break;
						}

					}
					else
					{
						reads_overlap_info->both_stand_used[abs(next_read)] = 1;
						//don't do anything,since both strands are visited.

					}

				}
				tmp_it = tmp_it_n;
			}

		}
		else
		{

			int lb = reads_overlap_info->LeftOverlaps[new_node].size();
			if (stacked_nodes[new_node] == -1 && lb>0)
			{
				stacked_nodes[new_node] = -1 - lb;
			}
			if (lb == 0  )//&& RemoveTipsOnly
			{
				stacked_nodes[new_node] = -2;
				if (Visited_Path[abs(new_node)].depth>((max_depth+1)/2))// let's define a tip to be non-extendable (2)
				{
					continue;
				}

				if (abs(new_node) == abs(beg_read))
				{
					continue;
				}
				//tip end reached so backtrack to the branching position.
				if (!isSimplePath_read(reads_overlap_info, new_node, Visited_Path, stacked_nodes))
				{
					continue;
				}
				//reads_overlap_info->contained[abs(new_node)] = 1;
				BacktrackBubbleRemoval_read(reads_overlap_info, new_node, beg_read, Visited_Path, stacked_nodes);
				if (beg_read>0)
				{
					reads_overlap_info->used_vt_right[abs(beg_read)] = 1;
				}
				else
				{
					reads_overlap_info->used_vt_left[abs(beg_read)] = 1;
				}
				stacked_nodes[new_node] = -1;


				///////////////latest modification
				if (reads_overlap_info->RightOverlaps[new_node].size())
				{
					//dist_reads[0].push_back(new_node);
				}

				break;

			}

			map<int,int >::iterator tmp_it, tmp_it_n;

			for (tmp_it = reads_overlap_info->LeftOverlaps[new_node].begin(); tmp_it != reads_overlap_info->LeftOverlaps[new_node].end();)
			{
				tmp_it_n = tmp_it;
				tmp_it_n++;
				int next_read = tmp_it->first;
				// not in stack
				if (stacked_nodes[abs(next_read)] == 0)
				{
					//					Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
					//Visited_Path[abs(next_read)].cov=(int)(Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]);
					Visited_Path[abs(next_read)].cov = Visited_Path[abs(new_node)].cov + tmp_it->second;//(int)(Visited_Path[abs(new_node)].cov);// +reads_overlap_info->left_overlap_size[new_node][next_read]);

					Visited_Path[abs(next_read)].depth = (Visited_Path[abs(new_node)].depth + 1);
					//int cum_len=(int)(Visited_Path[abs(next_read)].len+1);
					int cum_len = Visited_Path[abs(next_read)].depth;
					//Visited_Path[abs(next_read)].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);
					Visited_Path[abs(next_read)].last_read = new_node;

					if (next_read<0)
					{
						stacked_nodes[abs(next_read)] = 1;
					}
					else
					{
						stacked_nodes[abs(next_read)] = -1;
					}
					dist_reads[cum_len].push_back(-next_read);
					NBs++;

				}
				else
				{
					if (!RemoveTipsOnly)
					if ((stacked_nodes[abs(next_read)]<0 && next_read>0) || (stacked_nodes[abs(next_read)]>0 && next_read<0))
					{

						//if((Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]<=Visited_Path[abs(next_read)].cov))//||(BackCheckLoop(*ptr,new_node,Visited_Path)==1)//loop
						if ((Visited_Path[abs(new_node)].cov + tmp_it->second <= Visited_Path[abs(next_read)].cov) || (BackCheckLoop_read(new_node, Visited_Path) == 1))//||(BackCheckLoop(*ptr,new_node,Visited_Path)==1)//loop
						{

							//backtrack if the same direction is found

							if (!isSimplePath_read(reads_overlap_info, new_node, Visited_Path, stacked_nodes))
							{
								tmp_it = tmp_it_n;
								continue;
							}

							//backtrack the current path, common operation in this search
							if (stacked_nodes[new_node]>2)
							{
								BreakLinks_read(reads_overlap_info, stacked_nodes, new_node, next_read);
								tmp_it = tmp_it_n;
								break;
							}
							else
							{
								if (stacked_nodes[new_node]<-2)
								{
									BreakLinks_read(reads_overlap_info, stacked_nodes, new_node, next_read);
									tmp_it = tmp_it_n;
									break;
								}
								else
								{
									int free_node = (new_node);

									if (free_node == beg_read)
									{
										BreakLinks_read(reads_overlap_info, stacked_nodes, beg_read, new_node);
									}
									//reads_overlap_info->contained[abs(new_node)]=1;
									//free the node and edge.
									BreakLinks_read(reads_overlap_info, stacked_nodes, new_node, next_read);
									BacktrackBubbleRemoval_read(reads_overlap_info, new_node, beg_read, Visited_Path, stacked_nodes);

									tmp_it = tmp_it_n;
									break;
								}
							}
						}
						else
						{
							
							//backtrack the original path, rare operation in this search, can lead to errors
							if (!isSimplePath_read(reads_overlap_info, next_read, Visited_Path, stacked_nodes))
							{
								tmp_it = tmp_it_n;
								continue;
							}
							BacktrackBubbleRemoval_read(reads_overlap_info, next_read, beg_read, Visited_Path, stacked_nodes);

							//marginal 
							//Visited_Path[abs(next_read)].cov=(int)(Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]);
							Visited_Path[abs(next_read)].cov = (int)(Visited_Path[abs(new_node)].cov) + tmp_it->second;

							Visited_Path[abs(next_read)].depth = Visited_Path[abs(new_node)].depth + 1;
							Visited_Path[abs(next_read)].len = Visited_Path[abs(new_node)].depth;
							//marginal 

							Visited_Path[abs(next_read)].last_read = abs(new_node);
							tmp_it = tmp_it_n;
							break;
						}

					}
					else
					{
						reads_overlap_info->both_stand_used[abs(next_read)] = 1;
						//don't do anything,since both strands are visited.

					}

				}

				tmp_it = tmp_it_n;

			}

		}

	}
	return;
}



void aggressive_cleaning(reads_info *reads_info )
{

	//just a trick to clean the graph
	
	int n_removed_nodes = 0;
	/*
	for (int i = 1; i < reads_info->LeftOverlapsBest.size(); ++i)
	{
		int count = 0;
		if (i == 140270 || i==8290)
		{
			cout<<"";

		}
		if (reads_info->LeftOverlapsBest[i] != 0)
		{
			int next_read = reads_info->LeftOverlapsBest[i];

			if (next_read>0)
			{
				if (reads_info->RightOverlapsBest[next_read] != i)
				{
					count++;
				}
			}
			else
			{

				if (reads_info->LeftOverlapsBest[-next_read] != -i)
				{
					count++;
				}

			}
		}

		if (reads_info->RightOverlapsBest[i] != 0)
		{
			int next_read = reads_info->RightOverlapsBest[i];
			if (next_read>0)
			{

				if (reads_info->LeftOverlapsBest[next_read] != i)
				{
					count++;
				}
			}
			else
			{
				if (reads_info->RightOverlapsBest[-next_read] != -i)
				{
					count++;
				}

			}
		}
		if (count == 2)
		{
			reads_info->contained[i] = 1;
			reads_info->LeftOverlapsBest[i] = 0;
			reads_info->RightOverlapsBest[i] = 0;
			//cout << i << ", ";
			n_removed_nodes++;
			reads_info->degree_vec[abs(i)] = 0;
		}

	}
	*/
	//cout << endl;
	//cout << n_removed_nodes << " nodes removed in aggressive cleaning." << endl;

	n_removed_nodes = 0;
	
	for (int i = 1; i < reads_info->LeftOverlapsBestWithTies.size(); ++i)
	{
		int count = 0;

		
		for (int j = 0; j < reads_info->LeftOverlapsBestWithTies[i].size(); ++j)
		{
			
			int next_read = reads_info->LeftOverlapsBestWithTies[i][j];

			if (next_read>0)
			{
				bool kill = 1;
				for (int k = 0; k != reads_info->RightOverlapsBestWithTies[next_read].size(); ++k)
				{
					if (reads_info->RightOverlapsBestWithTies[next_read][k] == i)
					{
						kill = 0;
						break;
					}
				}
				if (kill)
				{
					count++;
				}
			}
			else
			{
				bool kill = 1;
				for (int k = 0; k != reads_info->LeftOverlapsBestWithTies[-next_read].size(); ++k)
				{
					if (reads_info->LeftOverlapsBestWithTies[-next_read][k] == -i)
					{
						kill = 0;
						break;
					}
				}
				if (kill)
				{
					count++;
				}

			
			}
		}
		
		for (int j = 0; j < reads_info->RightOverlapsBestWithTies[i].size(); ++j)
		{
			
			
			int next_read = reads_info->RightOverlapsBestWithTies[i][j];
			if (next_read>0)
			{
				bool kill = 1;
				for (int k = 0; k != reads_info->LeftOverlapsBestWithTies[next_read].size(); ++k)
				{
					if (reads_info->LeftOverlapsBestWithTies[next_read][k] == i)
					{
						kill = 0;
						break;
					}
				}
				if (kill)
				{
					count++;
				}

			}
			else
			{

				bool kill = 1;
				for (int k = 0; k != reads_info->RightOverlapsBestWithTies[-next_read].size(); ++k)
				{
					if (reads_info->RightOverlapsBestWithTies[-next_read][k] == -i)
					{
						kill = 0;
						break;
					}
				}
				if (kill)
				{
					count++;
				}


			}
		}
		if (count >= 2)
		{
			reads_info->contained[i] = 1;
			reads_info->LeftOverlapsBest[i] = 0;
			reads_info->RightOverlapsBest[i] = 0;
			reads_info->LeftOverlapsBestWithTies[i].clear();
			reads_info->RightOverlapsBestWithTies[i].clear();
			
			n_removed_nodes++;
			reads_info->degree_vec[abs(i)] = 0;
			//cout << i << endl;
		}

	}

	int n_removed_nodes2 = 0;
	for (int i = 1; i < reads_info->both_stand_used.size(); ++i)
	{
		if (reads_info->both_stand_used[i])
		{
			n_removed_nodes2++;
			//reads_info->contained[i] = 1;//latest modification
			//cout << i << endl;
		}
	}
	cout << n_removed_nodes << " bad nodes removed in aggressive cleaning." << endl;
	//cout << n_removed_nodes2 << " nodes removed because both strands are used." << endl;
}



#endif