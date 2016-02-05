#ifndef __ALIGH_H
#define __ALIGH_H


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
using namespace std;

struct align_matrices
{
	int score[152][152];
	int16_t path[152][152];
};



struct align_maps
{
	map<int, map<int, int> > score, path;
};
struct align_info
{

	bool SparseAlign;
	int max_len;
	bool Debug;
	bool flip,success,force_flip;
	vector<int> ref_aligned, qry_aligned;
	vector<int> right_ext,right_dist,right_mismatch;
	int ref_len,qry_len, match_len,mismatch_len[2],frontal_mismatch_len[4],rear_mismatch_len[4],max_score,ref_len_compressed,qry_len_compressed;
	//FR_mismatch_len :0 ref mismatch score, 1: qry mismatch score, 2: ref mismatch region len, 3: qry mismatch region len
	int ref_aligned_center,qry_aligned_center;
	int aligned_nodes;
	int min_overlap, max_mismatch_FR_score, min_extension_FR, min_match_qry, max_match_qry, min_match_ref, max_match_ref;
	int offset;
	int AbsOffsetSum, OffsetSum, MaxOffset, MaxOffsetRQ, MaxMismatchFR;
	double mm_ratio,mm_penalty;
	int band_width;
	char adj;int ahg;int bhg;//AMOS,Celera format
	int ref_idx;
	int top_reads;
	bool fix_orientation;
	int min_match, max_match;
	int aligned_range_ref, aligned_range_qry;
	bool FAST;
	bool CalculateOffset;
	bool Inward, Outward;
	int matching_method;
	int scoring_method;

	
};

int classify_alignment(align_info *align_info)
{
	/*
	classify into 5 cases:
	case 1: ref->qry overlap
	f_mm1>>f_mm2, r_mm2>>r_mm1, m>>mm
	case 2: qry->ref overlap
	case 3: qry in ref

	case 4: ref in qry
	case 5: no relation
	*/

	if (align_info->matching_method == 3)
	{
		if (align_info->scoring_method == 1)
		{
			align_info->max_score = max(align_info->aligned_range_ref, align_info->aligned_range_qry);
		}
		if (align_info->scoring_method == 2)
		{
			align_info->max_score = min(align_info->aligned_range_ref, align_info->aligned_range_qry);
		}

		if (align_info->scoring_method == 3)
		{
			align_info->max_score = (align_info->aligned_range_ref + align_info->aligned_range_qry) / 2;
		}
	}
	

	int result = 0;
	if ((align_info->match_len < align_info->min_overlap) || (((double)align_info->match_len*align_info->mm_ratio) < max(align_info->mismatch_len[0], align_info->mismatch_len[1])))
	{
		result = 5;
		return result;
	}
	
	if ((align_info->MaxMismatchFR)<min(align_info->frontal_mismatch_len[2], align_info->frontal_mismatch_len[3]) || (align_info->MaxMismatchFR)<min(align_info->rear_mismatch_len[2], align_info->rear_mismatch_len[3]))
	{
		//result = 5;
		//return result;
	}

	int AbsOffsetSum = align_info->AbsOffsetSum; 
	int OffsetSum = align_info->OffsetSum;
	int MaxOffsetRQ = align_info->MaxOffsetRQ;
	if (AbsOffsetSum>1)
	{
		cout << "";
	}
	if (OffsetSum >align_info->MaxOffset)// align_info->mm_ratio*min(align_info->aligned_range_ref, align_info->aligned_range_qry)
	{
		//cout << AbsOffsetSum << endl;
		//cout << "R: " << align_info->aligned_range_ref<<" ";// << endl;
		//cout <<"Q: "<< align_info->aligned_range_qry << endl;
		//cout << align_info->MaxOffset <<  endl;
		result = 5;
		//cout << "MM:" << align_info->aligned_range_ref << " " << align_info->aligned_range_qry << endl;
		return result;
	}
	if (((double)align_info->match_len*align_info->mm_ratio)<min(align_info->frontal_mismatch_len[0], align_info->frontal_mismatch_len[1]) || ((double)align_info->match_len*align_info->mm_ratio)<min(align_info->rear_mismatch_len[0], align_info->rear_mismatch_len[1]))
	{
		result = 5;
		return result;
	}

	if ((align_info->max_mismatch_FR_score)<min(align_info->frontal_mismatch_len[0], align_info->frontal_mismatch_len[1]) || (align_info->max_mismatch_FR_score)<min(align_info->rear_mismatch_len[0], align_info->rear_mismatch_len[1]))
	{
		result = 5;
		return result;
	}
	if (align_info->success == 0)
	{
		result = 5;
		return result;
	}

	if (max(align_info->frontal_mismatch_len[0], align_info->frontal_mismatch_len[1])<align_info->min_extension_FR)
	{
		align_info->frontal_mismatch_len[0] = 0;
		align_info->frontal_mismatch_len[1] = 0;
	}

	if (max(align_info->rear_mismatch_len[0], align_info->rear_mismatch_len[1])<align_info->min_extension_FR)
	{
		align_info->rear_mismatch_len[0] = 0;
		align_info->rear_mismatch_len[1] = 0;
	}

	//case 1 check:
	if ((align_info->frontal_mismatch_len[0] > align_info->frontal_mismatch_len[1]) && (align_info->frontal_mismatch_len[1]<align_info->max_mismatch_FR_score) && (align_info->frontal_mismatch_len[0]>align_info->min_extension_FR))
	{
		if ((align_info->rear_mismatch_len[0] < align_info->rear_mismatch_len[1]) && (align_info->rear_mismatch_len[0]<align_info->max_mismatch_FR_score) && (align_info->rear_mismatch_len[1]>align_info->min_extension_FR))
		{
			result = 1;
			return result;
		}
	}
	//case 2 check:
	if ((align_info->frontal_mismatch_len[0] < align_info->frontal_mismatch_len[1]) && (align_info->frontal_mismatch_len[0]<align_info->max_mismatch_FR_score) && (align_info->frontal_mismatch_len[1]>align_info->min_extension_FR))
	{
		if ((align_info->rear_mismatch_len[1] < align_info->rear_mismatch_len[0]) && (align_info->rear_mismatch_len[1]<align_info->max_mismatch_FR_score) && (align_info->rear_mismatch_len[0]>align_info->min_extension_FR))
		{
			result = 2;
			return result;
		}
	}
	//case 3 check:
	//may modify condition to enforce strict contain
	if ((align_info->frontal_mismatch_len[0] >= align_info->frontal_mismatch_len[1]) && (align_info->rear_mismatch_len[0] >= align_info->rear_mismatch_len[1]))// this is almost strict contain
	{
		if ((max(align_info->rear_mismatch_len[1], align_info->frontal_mismatch_len[1])<align_info->max_mismatch_FR_score))
		{

			result = 3;
			return result;
		}
	}
	//case 4 check:
	if ((align_info->frontal_mismatch_len[1] >= align_info->frontal_mismatch_len[0]) && (align_info->rear_mismatch_len[1] >= align_info->rear_mismatch_len[0]))
	{
		if ((max(align_info->rear_mismatch_len[0], align_info->frontal_mismatch_len[0])<align_info->max_mismatch_FR_score))
		{
			//may add condition to enforce strict contain
			result = 4;
			return result;
		}
	}

	return 0;
}



bool sparse_semi_global_align(contigs_info *contigs_info, vector<Coord_CTG_Cov>  ref, vector<Coord_CTG_Cov>  qry, align_matrices *align_matrices, align_info *align_info)
{
	memset(align_matrices->score, 0, sizeof(align_matrices->score));
	memset(align_matrices->path, 0, sizeof(align_matrices->path));
	align_info->flip = 0;
	align_info->success = 0;
	align_info->ref_aligned.clear();
	align_info->qry_aligned.clear();
	align_info->right_ext.clear();
	align_info->right_mismatch.clear();
	int max_len = align_info->max_len;
	align_info->CalculateOffset = 1;

	if (align_info->force_flip)
	{
		align_info->flip = 1;
	}
	
	/*
	semi global alignment,
	if two contigs match, the score is the contig length,
	mismatch: not defined.
	gap: penalty score is the contig size
	*/
	map<int, int> shared_index, ref_index, qry_index;
	int ref_sz = ref.size(), qry_sz = qry.size();

	for (int i = 0; i < ref_sz; ++i)
	{
		shared_index[ref[i].contig_no] = 1;
	}

	if (align_info->fix_orientation == 0)
	{
		
		int flip = 0, non_flip = 0;

		for (int i = 0; i < qry_sz; ++i)
		{
			if (shared_index.count(qry[i].contig_no))
			{
				non_flip++;
			}
			if (shared_index.count(-qry[i].contig_no))
			{
				flip++;
			}

		}


		if (flip == non_flip)
		{
			return 0;
		}

		if (flip > non_flip)
		{
			align_info->flip = 1;

		}
	}

	if (align_info->flip)
	{
		reverse(qry.begin(), qry.end());
		for (int i = 0; i < qry_sz; ++i)
		{
			qry[i].contig_no = -qry[i].contig_no;
			//also change coords
		}
		for (int i = 0; i < qry_sz; ++i)
		{
			qry[i].coord2 = align_info->qry_len - qry[i].coord2;
		}
	}

	for (int i = 0; i < qry_sz; ++i)
	{
		if (shared_index.count(qry[i].contig_no))
		{
			shared_index[qry[i].contig_no] = 2;
		}
	}



	//put matches idx into new structure, calculate match/gap score;
	vector<Coord_CTG_Cov>  ref_compressed, qry_compressed,ref_gap_compressed,qry_gap_compressed;
	ref_compressed.reserve(max_len);
	qry_compressed.reserve(max_len);
	ref_gap_compressed.reserve(max_len);
	qry_gap_compressed.reserve(max_len);
	
	Coord_CTG_Cov Coord_CTG_Cov_null;
	Coord_CTG_Cov_null.contig_no = 0;
	Coord_CTG_Cov_null.cov = 0;
	Coord_CTG_Cov_null.coord = 0;
	Coord_CTG_Cov_null.coord2 = 0;
	ref_gap_compressed.push_back(Coord_CTG_Cov_null);
	qry_gap_compressed.push_back(Coord_CTG_Cov_null);

	int n_ref_match = 0;
	ref_index.clear();

	for (int i = 0; i < ref_sz; ++i)
	{
		if (shared_index[ref[i].contig_no]==2)
		{
			ref_compressed.push_back(ref[i]);
			ref_gap_compressed.push_back(Coord_CTG_Cov_null);
			ref_index[n_ref_match] = i;
			n_ref_match++;
			
		}
		else
		{
			if (align_info->matching_method == 1)
			{
				ref_gap_compressed[n_ref_match].cov += contigs_info->contig_sz_vt[abs(ref[i].contig_no)];
			}
			if (align_info->matching_method == 2 || align_info->matching_method == 3)
			{
				ref_gap_compressed[n_ref_match].cov += ref[i].cov;
			}			
		}		
	}
	ref_compressed.push_back(Coord_CTG_Cov_null);// last one, stands for gaps;
	//ref_compressed[ref_compressed.size() - 1].cov = -ref_compressed[ref_compressed.size() - 1].cov;
	ref_index[n_ref_match] = ref_sz-1;
	ref_index[-1] = - 1;

	int n_qry_match = 0;
	qry_index.clear();
	for (int i = 0; i < qry_sz; ++i)
	{
		if (shared_index[qry[i].contig_no]==2)
		{
			qry_compressed.push_back(qry[i]);
			qry_gap_compressed.push_back(Coord_CTG_Cov_null);
			qry_index[n_qry_match] = i;
			n_qry_match++;
		}
		else
		{
			if (align_info->matching_method == 1)
			{
				qry_gap_compressed[n_qry_match].cov += contigs_info->contig_sz_vt[abs(qry[i].contig_no)];
			}
			if (align_info->matching_method == 2 || align_info->matching_method == 3)
			{
				qry_gap_compressed[n_qry_match].cov += qry[i].cov;
			}
		}
	}

	qry_compressed.push_back(Coord_CTG_Cov_null);// last one, stands for gaps;
	//qry_compressed[qry_compressed.size() - 1].cov = -qry_compressed[qry_compressed.size() - 1].cov;
	qry_index[n_qry_match] = qry_sz - 1;
	qry_index[-1] = - 1;


	if (ref_compressed.size() > max_len || qry_compressed.size() > max_len)
	{
		return 0;// no result
	}

	//calculate the sparse alignment using compressed reads;
	ref_sz = ref_compressed.size();
	qry_sz = qry_compressed.size();

	align_info->ref_len_compressed = ref_sz;
	align_info->qry_len_compressed = qry_sz;

	if (ref_sz > 2 && qry_sz > 2)
	{
		cout << "";
	}


	for (int r = 0; r < ref_sz; ++r)
	{
		align_matrices->path[r][0] = 2;
		for (int c = 0; c < qry_sz; ++c)
		{
			//
			int max_score=0;
			if (r>0)
			{
				max_score = align_matrices->score[r - 1][c] - align_info->mm_penalty*(double)(ref_compressed[r].cov + ref_gap_compressed[r].cov);
			}
			align_matrices->path[r][c] = 1;
			int tmp_score;
			if (c > 0)
			{
				tmp_score = align_matrices->score[r][c - 1] - align_info->mm_penalty*(double)(qry_compressed[c].cov + qry_gap_compressed[c].cov);
			}
			else
			{
				
				tmp_score = 0;
			}
			
			if (tmp_score > max_score)
			{
				align_matrices->path[r][c] = 2;
				max_score = tmp_score;
			}

			if (qry_sz-1!=c&&ref_sz-1!=r)
			if (ref_compressed[r ].contig_no == qry_compressed[c].contig_no)
			{				
				// here use min, max_len use max
				//use average resolves more ties
				if (r > 0 && c > 0)
				{
					tmp_score = align_matrices->score[r - 1][c - 1];
				}
				else
				{
					tmp_score = 0;
				}
				if (align_info->scoring_method == 1)
				{
					tmp_score =  tmp_score+ max(ref_compressed[r].cov, qry_compressed[c].cov);
				}
				if (align_info->scoring_method == 2)
				{

					tmp_score = tmp_score + min(ref_compressed[r].cov, qry_compressed[c].cov);
				}
				if (align_info->scoring_method == 3)
				{

					tmp_score = tmp_score + (ref_compressed[r].cov + qry_compressed[c].cov) / 2;
				}
				if (r > 0 && c > 0)
				{
					tmp_score -= align_info->mm_penalty*(double)(ref_gap_compressed[r].cov + qry_gap_compressed[c].cov);
				}
				
			}

			if (tmp_score > max_score)
			{
				align_matrices->path[r][c] = 3;
				max_score = tmp_score;
			}
			align_matrices->score[r][c] = max_score;
		}
	}


	//return alignement result
	int max_score = -1000000;
	int opt_r = 0, opt_c = 0;
	for (int r = 0; r < ref_sz; ++r)
	{
		if (align_matrices->score[r][qry_sz-1] > max_score)
		{
			max_score = align_matrices->score[r][qry_sz-1];
			opt_r = r;
			opt_c = qry_sz-1;
		}
	}
	for (int c = 0; c < qry_sz; ++c)
	{
		if (align_matrices->score[ref_sz-1][c] > max_score)
		{
			max_score = align_matrices->score[ref_sz-1][c];
			opt_r = ref_sz-1;
			opt_c = c;
		}
	}

	if(ref_sz-1 > opt_r)
	{
		for (int pp = ref_index[ref_sz-1]; pp > ref_index[opt_r]; --pp)
		{
			align_info->qry_aligned.push_back(0);
			int contig_no = ref[pp].contig_no;
			align_info->ref_aligned.push_back(contig_no);
		}		
	}

	if(qry_sz-1 > opt_c)
	{
		for (int pp = qry_index[qry_sz-1]; pp > qry_index[opt_c]; --pp)
		{
			align_info->ref_aligned.push_back(0);
			int contig_no = qry[pp].contig_no;
			align_info->qry_aligned.push_back(contig_no);
		}
	}

	while (1)
	{
		if (opt_r == -1 || opt_c == -1)
		{
			if (opt_r>=0)
			for (int pp = ref_index[opt_r] ; pp >= 0; --pp)
			{

				align_info->qry_aligned.push_back(0);
				int contig_no = ref[pp].contig_no;
				align_info->ref_aligned.push_back(contig_no);
			}

			if (opt_c >= 0)
			for (int pp = qry_index[opt_c] ; pp >= 0; --pp)
			{
				align_info->ref_aligned.push_back(0);
				int contig_no = qry[pp].contig_no;
				align_info->qry_aligned.push_back(contig_no);
			}

			break;
		}

		if (align_matrices->path[opt_r][opt_c] == 1)
		{
			
			for (int pp = ref_index[opt_r]; pp > ref_index[opt_r - 1]; --pp)
			{
				align_info->qry_aligned.push_back(0);
				int contig_no = ref[pp].contig_no;
				align_info->ref_aligned.push_back(contig_no);
			}
			opt_r--;

		}
		if (align_matrices->path[opt_r][opt_c] == 2)
		{
			
			
			for (int pp = qry_index[opt_c]; pp > qry_index[opt_c-1]; --pp)
			{
				align_info->ref_aligned.push_back(0);
				int contig_no = qry[pp].contig_no;
				align_info->qry_aligned.push_back(contig_no);
			}
			opt_c--;

		}
		if (align_matrices->path[opt_r][opt_c] == 3)
		{
			
			int contig_no = ref[ref_index[opt_r]].contig_no;
			align_info->ref_aligned.push_back(contig_no);
			align_info->qry_aligned.push_back(contig_no);
			
			for (int pp = ref_index[opt_r]-1; pp > ref_index[opt_r - 1]; --pp)
			{
				align_info->qry_aligned.push_back(0);
				int contig_no = ref[pp].contig_no;
				align_info->ref_aligned.push_back(contig_no);
			}
			for (int pp = qry_index[opt_c]-1; pp > qry_index[opt_c - 1]; --pp)
			{
				align_info->ref_aligned.push_back(0);
				int contig_no = qry[pp].contig_no;
				align_info->qry_aligned.push_back(contig_no);
			}

			opt_r--; opt_c--;
		}
	}
	reverse(align_info->ref_aligned.begin(), align_info->ref_aligned.end());

	reverse(align_info->qry_aligned.begin(), align_info->qry_aligned.end());

	
	/*

	report statistics
	ref length
	qry length
	matching length
	matching vs non matching
	mismatch in the beginning and end

	*/
	ref_index.clear();
	qry_index.clear();
	int n_match = 0;
	for (int i = 0; i < align_info->ref_aligned.size(); ++i)
	{
		if (align_info->ref_aligned[i] != 0)
		{
			ref_index[i] = n_match;
			n_match++;
		}
	}

	n_match = 0;
	for (int i = 0; i < align_info->qry_aligned.size(); ++i)
	{
		if (align_info->qry_aligned[i] != 0)
		{
			qry_index[i] = n_match;
			n_match++;
		}
	}


	align_info->max_score = max_score;
	int min_match = -1, max_match = -1;
	for (int i = 0; i < align_info->ref_aligned.size(); ++i)
	{
		if (align_info->ref_aligned[i] == align_info->qry_aligned[i])
		{
			if (min_match < 0)
			{
				min_match = i;
			}
			if (max_match < i)
			{
				max_match = i;
			}
		}
	}
	if (min_match < 0 || max_match < 0)
	{
		align_info->success = 0;
		return 0;
	}
	align_info->success = 1;

	align_info->match_len = 0;
	align_info->mismatch_len[0] = 0;
	align_info->mismatch_len[1] = 0;
	align_info->rear_mismatch_len[0] = 0;
	align_info->rear_mismatch_len[1] = 0;
	align_info->rear_mismatch_len[2] = 0;
	align_info->rear_mismatch_len[3] = 0;
	align_info->frontal_mismatch_len[0] = 0;
	align_info->frontal_mismatch_len[1] = 0;
	align_info->frontal_mismatch_len[2] = 0;
	align_info->frontal_mismatch_len[3] = 0;
	align_info->min_match = min_match;
	align_info->max_match = max_match;
	int temp_range = contigs_info->contig_sz_vt[abs(ref[ref_index[min_match]].contig_no)] + contigs_info->contig_sz_vt[abs(ref[ref_index[max_match]].contig_no)];
	temp_range /= 2;
	align_info->aligned_range_ref = abs(ref[ref_index[max_match]].coord2 - ref[ref_index[min_match]].coord2) + temp_range;
	align_info->aligned_range_qry = abs(qry[qry_index[max_match]].coord2 - qry[qry_index[min_match]].coord2) + temp_range;


	if (abs(align_info->aligned_range_ref - align_info->aligned_range_qry) > align_info->mm_ratio*min(align_info->aligned_range_ref, align_info->aligned_range_qry))
	{
		//cout << "MM:" << align_info->aligned_range_ref << " " << align_info->aligned_range_qry << endl;
	}
	vector<int> offset_vec0, offset_vec;
	vector<int> ref_matched_positions0, qry_matched_positions0, ref_matched_positions, qry_matched_positions;
	align_info->ref_aligned_center = 0;
	align_info->qry_aligned_center = 0;
	align_info->aligned_nodes = 0;
	int max_cov = 0;
	int last_match = -1;
	align_info->MaxOffsetRQ = -1;
	align_info->AbsOffsetSum = 0;
	align_info->OffsetSum = abs(align_info->aligned_range_ref - align_info->aligned_range_qry);

	for (int i = min_match; i <= max_match; ++i)
	{
		if (align_info->ref_aligned[i] == align_info->qry_aligned[i])
		{
			if (last_match >= 0)
			{
				int abs_offset = abs((int)abs(ref[ref_index[i]].coord2 - ref[ref_index[last_match]].coord2) - (int)abs(qry[qry_index[i]].coord2 - qry[qry_index[last_match]].coord2));
				if (abs_offset > align_info->MaxOffsetRQ)
				{
					align_info->MaxOffsetRQ = abs_offset;
					if (abs_offset > 500)
					{
						//	cout << abs(ref[ref_index[i]].coord2 - ref[ref_index[last_match]].coord2) << endl;

						//	cout << abs(qry[qry_index[i]].coord2 - qry[qry_index[last_match]].coord2) << endl;
					}
				}
				align_info->AbsOffsetSum += abs_offset;

			}
			last_match = i;
			if (align_info->matching_method == 1)
			{
				align_info->match_len += contigs_info->contig_sz_vt[abs(align_info->ref_aligned[i])];

			}
			if (align_info->matching_method == 2 || align_info->matching_method == 3)
			{

				if (ref_index.count(i) && qry_index.count(i))
				{
					//use average resolves more ties
					if (align_info->scoring_method == 1)
					{
						//cout << "1";
						align_info->match_len += max(ref[ref_index[i]].cov, qry[qry_index[i]].cov);
					}
					if (align_info->scoring_method == 2)
					{
						//cout << "2";
						align_info->match_len += min(ref[ref_index[i]].cov, qry[qry_index[i]].cov);
					}
					if (align_info->scoring_method == 3)
					{
						//cout << "3";
						align_info->match_len += (ref[ref_index[i]].cov + qry[qry_index[i]].cov) / 2;
					}

				}
				else
				{
					cout << "Warning. Unexpected mismatch." << endl;
				}

			}
			if (align_info->CalculateOffset)
			{

				//int ctg_matched = align_info->qry_aligned[i];

				if (ref[ref_index[i]].cov > 1 && qry[qry_index[i]].cov > 1)
				{

					if (max_cov == 0 || max_cov< min(ref[ref_index[i]].cov, qry[qry_index[i]].cov))//very important. use > to minimize risk, use < to have best robustness 
					{
						if (max(ref[ref_index[i]].cov, qry[qry_index[i]].cov) <= contigs_info->contig_sz_vt[abs(ref[ref_index[i]].contig_no)])
						{
							//cout << "o: " << max_cov << "n: " << min(ref[ref_index2[ctg_matched]].cov, qry[qry_index2[ctg_matched]].cov) << endl;
							offset_vec0.push_back(ref[ref_index[i]].coord2 - qry[qry_index[i]].coord2);
							max_cov = min(ref[ref_index[i]].cov, qry[qry_index[i]].cov);
							ref_matched_positions0.push_back(ref[ref_index[i]].coord2);
							qry_matched_positions0.push_back(qry[qry_index[i]].coord2);
						}

					}
					align_info->aligned_nodes++;
				}
				offset_vec.push_back(ref[ref_index[i]].coord2 - qry[qry_index[i]].coord2);
				ref_matched_positions.push_back(ref[ref_index[i]].coord2);
				qry_matched_positions.push_back(qry[qry_index[i]].coord2);

			}

		}
		else
		{
			if (align_info->matching_method == 1)
			{
				align_info->mismatch_len[0] += contigs_info->contig_sz_vt[abs(align_info->ref_aligned[i])];
				align_info->mismatch_len[1] += contigs_info->contig_sz_vt[abs(align_info->qry_aligned[i])];
			}
			if (align_info->matching_method == 2 || align_info->matching_method == 3)
			{
				if (ref_index.count(i))
				{
					align_info->mismatch_len[0] += ref[ref_index[i]].cov;
				}
				if (qry_index.count(i))
				{
					align_info->mismatch_len[1] += qry[qry_index[i]].cov;
				}
			}
		}
	}

	int ref_beg = ref[0].coord2 + contigs_info->contig_sz_vt[abs(ref[0].contig_no)] / 2;
	int qry_beg = qry[0].coord2 + contigs_info->contig_sz_vt[abs(qry[0].contig_no)] / 2;
	align_info->frontal_mismatch_len[2] = ref[ref_index[min_match]].coord2 - contigs_info->contig_sz_vt[abs(align_info->ref_aligned[min_match])] / 2 - ref_beg;
	align_info->frontal_mismatch_len[3] = qry[qry_index[min_match]].coord2 - contigs_info->contig_sz_vt[abs(align_info->qry_aligned[min_match])] / 2 - qry_beg;
	int ref_len = ref[ref.size() - 1].coord2 - contigs_info->contig_sz_vt[abs(ref[ref.size() - 1].contig_no)] / 2;
	int qry_len = qry[qry.size() - 1].coord2 - contigs_info->contig_sz_vt[abs(qry[qry.size() - 1].contig_no)] / 2;
	align_info->rear_mismatch_len[2] = ref_len - ref[ref_index[max_match]].coord2 - contigs_info->contig_sz_vt[abs(align_info->ref_aligned[max_match])] / 2;
	align_info->rear_mismatch_len[3] = qry_len - qry[qry_index[max_match]].coord2 - contigs_info->contig_sz_vt[abs(align_info->qry_aligned[max_match])] / 2;
	/*
	if (ref_index[min_match] == 0)
	{
	align_info->frontal_mismatch_len[2] = 0;
	}
	if (qry_index[min_match] == 0)
	{
	align_info->frontal_mismatch_len[3] = 0;
	}
	if (ref_index[max_match]+1 == ref.size())
	{
	align_info->rear_mismatch_len[2] = 0;
	}
	if (qry_index[max_match] + 1 == qry.size())
	{
	align_info->rear_mismatch_len[3] = 0;
	}
	*/
	for (int i = 0; i < min_match; ++i)
	{
		if (align_info->matching_method == 1)
		{
			align_info->frontal_mismatch_len[0] += contigs_info->contig_sz_vt[abs(align_info->ref_aligned[i])];
			align_info->frontal_mismatch_len[1] += contigs_info->contig_sz_vt[abs(align_info->qry_aligned[i])];
		}
		if (align_info->matching_method == 2 || align_info->matching_method == 3)
		{
			if (ref_index.count(i))
			{
				align_info->frontal_mismatch_len[0] += ref[ref_index[i]].cov;
			}
			if (qry_index.count(i))
			{
				align_info->frontal_mismatch_len[1] += qry[qry_index[i]].cov;
			}
		}

	}
	for (int i = max_match + 1; i < align_info->ref_aligned.size(); ++i)
	{
		if (align_info->matching_method == 1)
		{
			align_info->rear_mismatch_len[0] += contigs_info->contig_sz_vt[abs(align_info->ref_aligned[i])];
			align_info->rear_mismatch_len[1] += contigs_info->contig_sz_vt[abs(align_info->qry_aligned[i])];
		}
		if (align_info->matching_method == 2 || align_info->matching_method == 3)
		{
			if (ref_index.count(i))
			{
				align_info->rear_mismatch_len[0] += ref[ref_index[i]].cov;
			}
			if (qry_index.count(i))
			{
				align_info->rear_mismatch_len[1] += qry[qry_index[i]].cov;
			}
		}
	}
	for (int i = max_match + 1; i < align_info->ref_aligned.size(); ++i)
	{
		if (align_info->qry_aligned[i] != 0)
		{
			align_info->right_ext.push_back(align_info->qry_aligned[i]);
		}
		if (align_info->ref_aligned[i] != 0)
		{
			align_info->right_mismatch.push_back(align_info->ref_aligned[i]);
		}

	}
	int ext_len = align_info->right_ext.size();
	align_info->right_dist.resize(ext_len);
	if (ext_len<qry_sz)
	{
		for (int i = 0; i < ext_len; ++i)
		{
			int dist = (qry[qry_sz - ext_len + i].coord2 - qry[qry_sz - ext_len - 1 + i].coord2);
			align_info->right_dist[i] = abs(dist);
		}
	}


	if (align_info->CalculateOffset)
	{
		if (offset_vec0.size()>0)
		{
			//cout << "catch." << endl;
			//cout << offset_vec0.size() << endl;
			//cout << ref_matched_positions0.size() << endl;
			offset_vec = offset_vec0;
			ref_matched_positions = ref_matched_positions0;
			qry_matched_positions = qry_matched_positions0;


		}

		align_info->aligned_nodes = ref_matched_positions.size();

		sort(offset_vec.begin(), offset_vec.end());



		if (align_info->aligned_nodes>0)
		{

			align_info->ref_aligned_center = ref_matched_positions[ref_matched_positions.size() - 1];
			align_info->qry_aligned_center = qry_matched_positions[qry_matched_positions.size() - 1];
			if (align_info->aligned_nodes > 1)
			{
				for (int jj = 0; jj < align_info->aligned_nodes; ++jj)
				{
					//cout << jj<<" ref: " << ref_matched_positions[jj] << " qry: " << qry_matched_positions [jj ]<< endl;

				}
			}
			align_info->aligned_nodes = 1;
		}


		align_info->offset = offset_vec[offset_vec.size() / 2];

		align_info->ahg = align_info->offset;
		align_info->bhg = align_info->qry_len - align_info->ref_len + align_info->offset;
		if (align_info->flip)
		{
			align_info->adj = 'I';
		}
		else
		{
			align_info->adj = 'N';
		}
		//if (align_info->ref_len == 10322 || align_info->qry_len == 10322 || align_info->ref_len == 10625 || align_info->qry_len == 10625)
		//cout << "offset: " << align_info->offset<<endl;
	}


	//align_info->max_score = align_info->max_score - align_info->mm_penalty*((double)min(align_info->frontal_mismatch_len[0], align_info->frontal_mismatch_len[1]));
	//align_info->max_score = align_info->max_score - align_info->mm_penalty*((double)min(align_info->rear_mismatch_len[0], align_info->rear_mismatch_len[1]));


	return 1;

}


bool approximate_semi_global_align(contigs_info *contigs_info, vector<Coord_CTG_Cov>  ref, vector<Coord_CTG_Cov>  qry, align_matrices *align_matrices, align_info *align_info)
{
	memset(align_matrices->score, 0, sizeof(align_matrices->score));
	memset(align_matrices->path, 0, sizeof(align_matrices->path));
	align_info->flip = 0;
	align_info->success = 0;
	align_info->ref_aligned.clear();
	align_info->qry_aligned.clear();
	align_info->right_ext.clear();
	align_info->right_mismatch.clear();
	int max_len = align_info->max_len;
	align_info->CalculateOffset = 1;

	if (align_info->force_flip)
	{
		align_info->flip = 1;
	}
	/*
	TBD: preprocess if input fragments are too long, keep the best ones.
	*/

	if (ref.size() > max_len || qry.size() > max_len)
	{
		return 0;// no result
	}
	/*
	semi global alignment,
	if two contigs match, the score is the contig length,
	mismatch: not defined.
	gap: penalty score is the contig size
	*/
	map<int, int> local_index, ref_index, qry_index;
	int ref_sz = ref.size(), qry_sz = qry.size();
	for (int i = 0; i < ref_sz; ++i)
	{
		local_index[ref[i].contig_no] = ref[i].coord;
	}
	
	

	if (align_info->fix_orientation == 0)
	{
		int flip = 0, non_flip = 0;

		for (int i = 0; i < qry_sz; ++i)
		{
			if (local_index.count(qry[i].contig_no))
			{
				non_flip++;
			}
			if (local_index.count(-qry[i].contig_no))
			{
				flip++;
			}

		}

		
		if (flip == non_flip)
		{
			return 0;
		}

		if (flip > non_flip)
		{
			align_info->flip = 1;
			
		}
	}

	if (align_info->flip)
	{
		reverse(qry.begin(), qry.end());
		for (int i = 0; i < qry_sz; ++i)
		{
			qry[i].contig_no = -qry[i].contig_no;
			//also change coords
		}
		for (int i = 0; i < qry_sz; ++i)
		{
			qry[i].coord2 = align_info->qry_len - qry[i].coord2;
		}
	}


	for (int r = 1; r <= ref_sz; ++r)
	{

		for (int c = 1; c <= qry_sz; ++c)
		{
			//
			int max_score;
			if (align_info->matching_method == 1)
			{
				max_score = align_matrices->score[r - 1][c] - align_info->mm_penalty*(double)contigs_info->contig_sz_vt[abs(ref[r - 1].contig_no)];
			}
			if (align_info->matching_method == 2 || align_info->matching_method == 3)
			{
				max_score = align_matrices->score[r - 1][c] - align_info->mm_penalty*(double)abs(ref[r - 1].cov);
			}
			align_matrices->path[r][c] = 1;

			int tmp_score;
			if (align_info->matching_method == 1)
			{
				tmp_score = align_matrices->score[r][c - 1] - align_info->mm_penalty*(double)contigs_info->contig_sz_vt[abs(qry[c - 1].contig_no)];
			}
			if (align_info->matching_method == 2 || align_info->matching_method == 3)
			{
				tmp_score = align_matrices->score[r][c - 1] - align_info->mm_penalty*(double)abs(qry[c - 1].cov);
			}

			if (tmp_score > max_score)
			{
				align_matrices->path[r][c] = 2;
				max_score = tmp_score;
			}

			if (ref[r - 1].contig_no == qry[c - 1].contig_no)
			{
				if (align_info->matching_method == 1)
				{
					tmp_score = align_matrices->score[r - 1][c - 1] + contigs_info->contig_sz_vt[abs(ref[r - 1].contig_no)];
				}
				if (align_info->matching_method == 2 || align_info->matching_method == 3)
				{
					// here use min, max_len use max
					//use average resolves more ties
					if (align_info->scoring_method == 1)
					{
						
						tmp_score = align_matrices->score[r - 1][c - 1] + max(ref[r - 1].cov, qry[c - 1].cov);
					}
					if (align_info->scoring_method == 2)
					{
						
						tmp_score = align_matrices->score[r - 1][c - 1] + min(ref[r - 1].cov, qry[c - 1].cov);
					}
					if (align_info->scoring_method == 3)
					{
						
						tmp_score = align_matrices->score[r - 1][c - 1] + (ref[r - 1].cov + qry[c - 1].cov) / 2;
					}
				}
			}

			if (tmp_score > max_score)
			{
				align_matrices->path[r][c] = 3;
				max_score = tmp_score;
			}

			align_matrices->score[r][c] = max_score;

		}

	}


	//return alignement result
	int max_score = -1000000;
	int opt_r = 0, opt_c = 0;
	for (int r = 1; r <= ref_sz; ++r)
	{
		if (align_matrices->score[r][qry_sz] > max_score)
		{
			max_score = align_matrices->score[r][qry_sz];
			opt_r = r;
			opt_c = qry_sz;
		}
	}
	for (int c = 1; c <= qry_sz; ++c)
	{
		if (align_matrices->score[ref_sz][c] > max_score)
		{
			max_score = align_matrices->score[ref_sz][c];
			opt_r = ref_sz;
			opt_c = c;
		}
	}

	for (int pad = ref_sz - opt_r; pad > 0; --pad)
	{
		align_info->qry_aligned.push_back(0);
		int contig_no = ref[opt_r + pad - 1].contig_no;
		align_info->ref_aligned.push_back(contig_no);
	}
	for (int pad = qry_sz - opt_c; pad > 0; --pad)
	{
		align_info->ref_aligned.push_back(0);
		int contig_no = qry[opt_c + pad - 1].contig_no;
		align_info->qry_aligned.push_back(contig_no);
	}
	while (1)
	{
		if (opt_r == 0 || opt_c == 0)
		{
			for (int pad = opt_r; pad > 0; --pad)
			{
				align_info->qry_aligned.push_back(0);
				int contig_no = ref[pad - 1].contig_no;
				align_info->ref_aligned.push_back(contig_no);
			}
			for (int pad = opt_c; pad > 0; --pad)
			{
				align_info->ref_aligned.push_back(0);
				int contig_no = qry[pad - 1].contig_no;
				align_info->qry_aligned.push_back(contig_no);
			}

			break;
		}

		if (align_matrices->path[opt_r][opt_c] == 1)
		{
			align_info->qry_aligned.push_back(0);
			int contig_no = ref[opt_r - 1].contig_no;
			align_info->ref_aligned.push_back(contig_no);

			opt_r--;
		}
		if (align_matrices->path[opt_r][opt_c] == 2)
		{
			align_info->ref_aligned.push_back(0);
			int contig_no = qry[opt_c - 1].contig_no;
			align_info->qry_aligned.push_back(contig_no);
			opt_c--;
		}
		if (align_matrices->path[opt_r][opt_c] == 3)
		{

			int contig_no = ref[opt_r - 1].contig_no;
			align_info->ref_aligned.push_back(contig_no);
			align_info->qry_aligned.push_back(contig_no);
			opt_r--; opt_c--;
		}
	}
	reverse(align_info->ref_aligned.begin(), align_info->ref_aligned.end());

	reverse(align_info->qry_aligned.begin(), align_info->qry_aligned.end());



	/*

	report statistics
	ref length
	qry length
	matching length
	matching vs non matching
	mismatch in the beginning and end

	*/

	int n_match = 0;
	for (int i = 0; i < align_info->ref_aligned.size(); ++i)
	{
		if (align_info->ref_aligned[i] != 0)
		{
			ref_index[i] = n_match;
			n_match++;
		}
	}

	n_match = 0;
	for (int i = 0; i < align_info->qry_aligned.size(); ++i)
	{
		if (align_info->qry_aligned[i] != 0)
		{
			qry_index[i] = n_match;
			n_match++;
		}
	}


	align_info->max_score = max_score;
	int min_match = -1, max_match = -1;
	for (int i = 0; i < align_info->ref_aligned.size(); ++i)
	{
		if (align_info->ref_aligned[i] == align_info->qry_aligned[i])
		{
			if (min_match < 0)
			{
				min_match = i;
			}
			if (max_match < i)
			{
				max_match = i;
			}
		}
	}
	if (min_match < 0 || max_match < 0)
	{
		align_info->success = 0;
		return 0;
	}
	align_info->success = 1;

	align_info->match_len = 0;
	align_info->mismatch_len[0] = 0;
	align_info->mismatch_len[1] = 0;
	align_info->rear_mismatch_len[0] = 0;
	align_info->rear_mismatch_len[1] = 0;
	align_info->rear_mismatch_len[2] = 0;
	align_info->rear_mismatch_len[3] = 0;
	align_info->frontal_mismatch_len[0] = 0;
	align_info->frontal_mismatch_len[1] = 0;
	align_info->frontal_mismatch_len[2] = 0;
	align_info->frontal_mismatch_len[3] = 0;
	align_info->min_match = min_match;
	align_info->max_match = max_match;
	int temp_range = contigs_info->contig_sz_vt[abs(ref[ref_index[min_match]].contig_no)] + contigs_info->contig_sz_vt[abs(ref[ref_index[max_match]].contig_no)];
	temp_range /= 2;
	align_info->aligned_range_ref = abs(ref[ref_index[max_match]].coord2 - ref[ref_index[min_match]].coord2) + temp_range;
	align_info->aligned_range_qry = abs(qry[qry_index[max_match]].coord2 - qry[qry_index[min_match]].coord2) + temp_range;

	
	if (abs(align_info->aligned_range_ref - align_info->aligned_range_qry) > align_info->mm_ratio*min(align_info->aligned_range_ref, align_info->aligned_range_qry))
	{
		//cout << "MM:" << align_info->aligned_range_ref << " " << align_info->aligned_range_qry << endl;
	}
	vector<int> offset_vec0, offset_vec;
	vector<int> ref_matched_positions0, qry_matched_positions0, ref_matched_positions, qry_matched_positions;
	align_info->ref_aligned_center = 0;
	align_info->qry_aligned_center = 0;
	align_info->aligned_nodes = 0;
	int max_cov = 0;
	int last_match = -1;
	align_info->MaxOffsetRQ = -1;
	align_info->AbsOffsetSum = 0;
	align_info->OffsetSum = abs(align_info->aligned_range_ref - align_info->aligned_range_qry);

	for (int i = min_match; i <= max_match; ++i)
	{
		if (align_info->ref_aligned[i] == align_info->qry_aligned[i])
		{
			if (last_match >= 0)
			{
				int abs_offset = abs((int)abs(ref[ref_index[i]].coord2 - ref[ref_index[last_match]].coord2) - (int)abs(qry[qry_index[i]].coord2 - qry[qry_index[last_match]].coord2));
				if (abs_offset > align_info->MaxOffsetRQ)
				{
					align_info->MaxOffsetRQ = abs_offset;
					if (abs_offset > 500)
					{
					//	cout << abs(ref[ref_index[i]].coord2 - ref[ref_index[last_match]].coord2) << endl;

					//	cout << abs(qry[qry_index[i]].coord2 - qry[qry_index[last_match]].coord2) << endl;
					}
				}
				align_info->AbsOffsetSum += abs_offset;
				
			}
			last_match = i;
			if (align_info->matching_method == 1)
			{
				align_info->match_len += contigs_info->contig_sz_vt[abs(align_info->ref_aligned[i])];

			}
			if (align_info->matching_method == 2 || align_info->matching_method == 3)
			{

				if (ref_index.count(i) && qry_index.count(i))
				{
					//use average resolves more ties
					if (align_info->scoring_method == 1)
					{
						//cout << "1";
						align_info->match_len += max(ref[ref_index[i]].cov, qry[qry_index[i]].cov);
					}
					if (align_info->scoring_method == 2)
					{
						//cout << "2";
						align_info->match_len += min(ref[ref_index[i]].cov, qry[qry_index[i]].cov);
					}
					if (align_info->scoring_method == 3)
					{
						//cout << "3";
						align_info->match_len += (ref[ref_index[i]].cov + qry[qry_index[i]].cov) / 2;
					}

				}
				else
				{
					cout << "Warning. Unexpected mismatch." << endl;
				}

			}
			if (align_info->CalculateOffset)
			{

				//int ctg_matched = align_info->qry_aligned[i];

				if (ref[ref_index[i]].cov > 1 && qry[qry_index[i]].cov > 1)
				{

					if (max_cov == 0 || max_cov< min(ref[ref_index[i]].cov, qry[qry_index[i]].cov))//very important. use > to minimize risk, use < to have best robustness 
					{
						if (max(ref[ref_index[i]].cov, qry[qry_index[i]].cov)<=contigs_info->contig_sz_vt[abs(ref[ref_index[i]].contig_no)])
						{
							//cout << "o: " << max_cov << "n: " << min(ref[ref_index2[ctg_matched]].cov, qry[qry_index2[ctg_matched]].cov) << endl;
							offset_vec0.push_back(ref[ref_index[i]].coord2 - qry[qry_index[i]].coord2);
							max_cov = min(ref[ref_index[i]].cov, qry[qry_index[i]].cov);
							ref_matched_positions0.push_back(ref[ref_index[i]].coord2);
							qry_matched_positions0.push_back(qry[qry_index[i]].coord2);
						}

					}
					align_info->aligned_nodes++;
				}
				offset_vec.push_back(ref[ref_index[i]].coord2 - qry[qry_index[i]].coord2);
				ref_matched_positions.push_back(ref[ref_index[i]].coord2);
				qry_matched_positions.push_back(qry[qry_index[i]].coord2);

			}

		}
		else
		{
			if (align_info->matching_method == 1)
			{
				align_info->mismatch_len[0] += contigs_info->contig_sz_vt[abs(align_info->ref_aligned[i])];
				align_info->mismatch_len[1] += contigs_info->contig_sz_vt[abs(align_info->qry_aligned[i])];
			}
			if (align_info->matching_method == 2 || align_info->matching_method == 3)
			{
				if (ref_index.count(i))
				{
					align_info->mismatch_len[0] += ref[ref_index[i]].cov;
				}
				if (qry_index.count(i))
				{
					align_info->mismatch_len[1] += qry[qry_index[i]].cov;
				}
			}
		}
	}
	
	int ref_beg = ref[0].coord2 + contigs_info->contig_sz_vt[abs(ref[0].contig_no)] / 2;
	int qry_beg = qry[0].coord2 + contigs_info->contig_sz_vt[abs(qry[0].contig_no)] / 2;
	align_info->frontal_mismatch_len[2] = ref[ref_index[min_match]].coord2 - contigs_info->contig_sz_vt[abs(align_info->ref_aligned[min_match])] / 2 - ref_beg;
	align_info->frontal_mismatch_len[3] = qry[qry_index[min_match]].coord2 - contigs_info->contig_sz_vt[abs(align_info->qry_aligned[min_match])] / 2 - qry_beg;
	int ref_len = ref[ref.size() - 1].coord2 - contigs_info->contig_sz_vt[abs(ref[ref.size() - 1].contig_no)] / 2;
	int qry_len = qry[qry.size() - 1].coord2 - contigs_info->contig_sz_vt[abs(qry[qry.size() - 1].contig_no)] / 2;
	align_info->rear_mismatch_len[2] = ref_len - ref[ref_index[max_match]].coord2 - contigs_info->contig_sz_vt[abs(align_info->ref_aligned[max_match])] / 2;
	align_info->rear_mismatch_len[3] = qry_len - qry[qry_index[max_match]].coord2 - contigs_info->contig_sz_vt[abs(align_info->qry_aligned[max_match])] / 2;
	/*
	if (ref_index[min_match] == 0)
	{
		align_info->frontal_mismatch_len[2] = 0;
	}
	if (qry_index[min_match] == 0)
	{
		align_info->frontal_mismatch_len[3] = 0;
	}
	if (ref_index[max_match]+1 == ref.size())
	{
		align_info->rear_mismatch_len[2] = 0;
	}
	if (qry_index[max_match] + 1 == qry.size())
	{
		align_info->rear_mismatch_len[3] = 0;
	}
	*/
	for (int i = 0; i < min_match; ++i)
	{
		if (align_info->matching_method == 1)
		{
			align_info->frontal_mismatch_len[0] += contigs_info->contig_sz_vt[abs(align_info->ref_aligned[i])];
			align_info->frontal_mismatch_len[1] += contigs_info->contig_sz_vt[abs(align_info->qry_aligned[i])];
		}
		if (align_info->matching_method == 2 || align_info->matching_method == 3)
		{
			if (ref_index.count(i))
			{
				align_info->frontal_mismatch_len[0] += ref[ref_index[i]].cov;
			}
			if (qry_index.count(i))
			{
				align_info->frontal_mismatch_len[1] += qry[qry_index[i]].cov;
			}
		}

	}
	for (int i = max_match + 1; i < align_info->ref_aligned.size(); ++i)
	{
		if (align_info->matching_method == 1)
		{
			align_info->rear_mismatch_len[0] += contigs_info->contig_sz_vt[abs(align_info->ref_aligned[i])];
			align_info->rear_mismatch_len[1] += contigs_info->contig_sz_vt[abs(align_info->qry_aligned[i])];
		}
		if (align_info->matching_method == 2 || align_info->matching_method == 3)
		{
			if (ref_index.count(i))
			{
				align_info->rear_mismatch_len[0] += ref[ref_index[i]].cov;
			}
			if (qry_index.count(i))
			{
				align_info->rear_mismatch_len[1] += qry[qry_index[i]].cov;
			}
		}
	}
	for (int i = max_match + 1; i < align_info->ref_aligned.size(); ++i)
	{
		if (align_info->qry_aligned[i] != 0)
		{
			align_info->right_ext.push_back(align_info->qry_aligned[i]);
		}
		if (align_info->ref_aligned[i] != 0)
		{
			align_info->right_mismatch.push_back(align_info->ref_aligned[i]);
		}

	}
	int ext_len = align_info->right_ext.size();
	align_info->right_dist.resize(ext_len);
	if (ext_len<qry_sz)
	{
		for (int i = 0; i < ext_len; ++i)
		{
			int dist = (qry[qry_sz - ext_len + i].coord2 - qry[qry_sz - ext_len - 1 + i].coord2);
			align_info->right_dist[i] = abs(dist);
		}
	}


	if (align_info->CalculateOffset)
	{
		if (offset_vec0.size()>0)
		{
			//cout << "catch." << endl;
			//cout << offset_vec0.size() << endl;
			//cout << ref_matched_positions0.size() << endl;
			offset_vec = offset_vec0;
			ref_matched_positions = ref_matched_positions0;
			qry_matched_positions = qry_matched_positions0;


		}

		align_info->aligned_nodes = ref_matched_positions.size();

		sort(offset_vec.begin(), offset_vec.end());



		if (align_info->aligned_nodes>0)
		{

			align_info->ref_aligned_center = ref_matched_positions[ref_matched_positions.size() - 1];
			align_info->qry_aligned_center = qry_matched_positions[qry_matched_positions.size() - 1];
			if (align_info->aligned_nodes > 1)
			{
				for (int jj = 0; jj < align_info->aligned_nodes; ++jj)
				{
					//cout << jj<<" ref: " << ref_matched_positions[jj] << " qry: " << qry_matched_positions [jj ]<< endl;

				}
			}
			align_info->aligned_nodes = 1;
		}


		align_info->offset = offset_vec[offset_vec.size() / 2];

		align_info->ahg = align_info->offset;
		align_info->bhg = align_info->qry_len - align_info->ref_len + align_info->offset;
		if (align_info->flip)
		{
			align_info->adj = 'I';
		}
		else
		{
			align_info->adj = 'N';
		}
		//if (align_info->ref_len == 10322 || align_info->qry_len == 10322 || align_info->ref_len == 10625 || align_info->qry_len == 10625)
		//cout << "offset: " << align_info->offset<<endl;
	}


	//align_info->max_score = align_info->max_score - align_info->mm_penalty*((double)min(align_info->frontal_mismatch_len[0], align_info->frontal_mismatch_len[1]));
	//align_info->max_score = align_info->max_score - align_info->mm_penalty*((double)min(align_info->rear_mismatch_len[0], align_info->rear_mismatch_len[1]));


	return 1;

}


bool approximate_local_align_banded(contigs_info *contigs_info, vector<Coord_CTG_Cov>  ref, vector<Coord_CTG_Cov>  qry, align_maps *align_matrices, align_info *align_info)
{
	align_matrices->score.clear();
	align_matrices->path.clear();

	align_info->flip = 0;
	align_info->success = 0;
	align_info->ref_aligned.clear();
	align_info->qry_aligned.clear();
	align_info->right_ext.clear();
	int max_len = align_info->max_len;
	int band_width = align_info->band_width;
	/*
	TBD: preprocess if input fragments are too long, keep the best ones.
	*/

	if (ref.size() > max_len || qry.size() > max_len)
	{
		return 0;// no result
	}
	/*
	semi global alignment,
	if two contigs match, the score is the contig length,
	mismatch: not defined.
	gap: penalty score is the contig size
	*/
	map<int, int> local_index;
	int ref_sz = ref.size(), qry_sz = qry.size();
	for (int i = 0; i < ref_sz; ++i)
	{
		local_index[ref[i].contig_no] = ref[i].coord;
	}
	int flip = 0, non_flip = 0;

	for (int i = 0; i < qry_sz; ++i)
	{
		if (local_index.count(qry[i].contig_no))
		{
			non_flip++;
		}
		if (local_index.count(-qry[i].contig_no))
		{
			flip++;
		}

	}

	if (align_info->fix_orientation == 0)
	{
		if (flip == non_flip)
		{
			return 0;
		}

		if (flip > non_flip)
		{
			align_info->flip = 1;
			reverse(qry.begin(), qry.end());
			for (int i = 0; i < qry_sz; ++i)
			{
				qry[i].contig_no = -qry[i].contig_no;
				//also change coords
			}
		}
	}



	int idx1 = -1, idx2 = -1;

	for (int i = 0; i <= band_width; ++i)
	{
		align_matrices->score[i][0] = 0;

		align_matrices->score[0][i] = 0;
	}

	int opt_r = -1, opt_c = -1, opt_score = -1000000;
	for (int r = 1; r <= ref_sz; ++r)
	{

		for (int c = max(1, r - band_width); c <= min(qry_sz, r + band_width); ++c)
		{
			if (idx2<0 && r == ref_sz)
			{
				idx2 = c;
			}
			if (idx1 < 0 && c == qry_sz)
			{
				idx1 = r;
			}
			//
			int max_score = 0;

			int tmp_score = -10000000;
			align_matrices->path[r][c] = 0;

			if (align_matrices->score[r - 1].count(c))
			{
				tmp_score = align_matrices->score[r - 1][c] - align_info->mm_penalty*(double)contigs_info->contig_sz_vt[abs(ref[r - 1].contig_no)];
			}
			if (tmp_score > max_score)
			{
				align_matrices->path[r][c] = 1;
				max_score = tmp_score;
			}

			if (align_matrices->score[r].count(c - 1))
			{

				tmp_score = align_matrices->score[r][c - 1] - align_info->mm_penalty*(double)contigs_info->contig_sz_vt[abs(qry[c - 1].contig_no)];

			}
			if (tmp_score > max_score)
			{
				align_matrices->path[r][c] = 2;
				max_score = tmp_score;
			}

			if (ref[r - 1].contig_no == qry[c - 1].contig_no)
			{
				tmp_score = align_matrices->score[r - 1][c - 1] + contigs_info->contig_sz_vt[abs(ref[r - 1].contig_no)];
			}

			if (tmp_score > max_score)
			{
				align_matrices->path[r][c] = 3;
				max_score = tmp_score;
			}

			align_matrices->score[r][c] = max_score;

			if (max_score > opt_score)
			{
				opt_score = max_score;
				opt_r = r;
				opt_c = c;
			}

		}

	}

	//return alignement result

	for (int pad = ref_sz - opt_r; pad > 0; --pad)
	{
		align_info->qry_aligned.push_back(0);
		int contig_no = ref[opt_r + pad - 1].contig_no;
		align_info->ref_aligned.push_back(contig_no);
	}
	//delete
	/*
	for (int pad = qry_sz - opt_c; pad > 0; --pad)
	{
	align_info->ref_aligned.push_back(0);
	int contig_no = qry[opt_c + pad - 1].contig_no;
	align_info->qry_aligned.push_back(contig_no);
	}
	*/
	while (1)
	{
		if (opt_r == 0 || opt_c == 0 || align_matrices->score[opt_r][opt_c] == 0)
		{
			for (int pad = opt_r; pad > 0; --pad)
			{
				align_info->qry_aligned.push_back(0);
				int contig_no = ref[pad - 1].contig_no;
				align_info->ref_aligned.push_back(contig_no);
			}
			//delete
			/*
			for (int pad = opt_c; pad > 0; --pad)
			{
			align_info->ref_aligned.push_back(0);
			int contig_no = qry[pad - 1].contig_no;
			align_info->qry_aligned.push_back(contig_no);
			}
			*/
			break;
		}
		if (align_matrices->path[opt_r][opt_c] == 0)
		{
			break;
		}
		if (align_matrices->path[opt_r][opt_c] == 1)
		{
			align_info->qry_aligned.push_back(0);
			int contig_no = ref[opt_r - 1].contig_no;
			align_info->ref_aligned.push_back(contig_no);

			opt_r--;
		}
		if (align_matrices->path[opt_r][opt_c] == 2)
		{
			align_info->ref_aligned.push_back(0);
			int contig_no = qry[opt_c - 1].contig_no;
			align_info->qry_aligned.push_back(contig_no);
			opt_c--;
		}
		if (align_matrices->path[opt_r][opt_c] == 3)
		{

			int contig_no = ref[opt_r - 1].contig_no;
			align_info->ref_aligned.push_back(contig_no);
			align_info->qry_aligned.push_back(contig_no);
			opt_r--; opt_c--;
		}

	}
	reverse(align_info->ref_aligned.begin(), align_info->ref_aligned.end());

	reverse(align_info->qry_aligned.begin(), align_info->qry_aligned.end());


	/*

	report statistics
	ref length
	qry length
	matching length
	matching vs non matching
	mismatch in the beginning and end

	*/
	align_info->max_score = opt_score;
	int min_match = -1, max_match = -1;
	for (int i = 0; i < align_info->ref_aligned.size(); ++i)
	{
		if (align_info->ref_aligned[i] == align_info->qry_aligned[i])
		{
			if (min_match < 0)
			{
				min_match = i;
			}
			if (max_match < i)
			{
				max_match = i;
			}
		}
	}
	if (min_match < 0 || max_match < 0)
	{
		align_info->success = 0;
		return 0;
	}
	align_info->success = 1;

	align_info->match_len = 0;
	align_info->mismatch_len[0] = 0;
	align_info->mismatch_len[1] = 0;
	align_info->rear_mismatch_len[0] = 0;
	align_info->rear_mismatch_len[1] = 0;
	align_info->frontal_mismatch_len[0] = 0;
	align_info->frontal_mismatch_len[1] = 0;
	align_info->min_match = min_match;
	align_info->max_match = max_match;

	for (int i = min_match; i <= max_match; ++i)
	{
		if (align_info->ref_aligned[i] == align_info->qry_aligned[i])
		{
			align_info->match_len += contigs_info->contig_sz_vt[abs(align_info->ref_aligned[i])];
		}
		else
		{
			align_info->mismatch_len[0] += contigs_info->contig_sz_vt[abs(align_info->ref_aligned[i])];
			align_info->mismatch_len[1] += contigs_info->contig_sz_vt[abs(align_info->qry_aligned[i])];
		}
	}
	for (int i = 0; i < min_match; ++i)
	{

		align_info->frontal_mismatch_len[0] += contigs_info->contig_sz_vt[abs(align_info->ref_aligned[i])];
		align_info->frontal_mismatch_len[1] += contigs_info->contig_sz_vt[abs(align_info->qry_aligned[i])];

	}
	for (int i = max_match + 1; i < align_info->ref_aligned.size(); ++i)
	{

		align_info->rear_mismatch_len[0] += contigs_info->contig_sz_vt[abs(align_info->ref_aligned[i])];
		align_info->rear_mismatch_len[1] += contigs_info->contig_sz_vt[abs(align_info->qry_aligned[i])];

	}


	return 1;

}



bool align2ref(contigs_info *contigs_info, reads_info *reads_info, vector<Coord_CTG_Cov>  qry, align_maps *align_matrices, vector<struct align_info> &align_info_vec)
{
	align_info align_info0 = align_info_vec.front();
	align_info_vec.clear();
	map<int, vector<int> > read_oi;// begin_position, end_position;
	map<int, int > read_score;

	map<int, int> local_index;
	vector< vector<int> > matching_contigs_vec;
	int qry_sz = qry.size();
	//cout <<"q"<< qry_sz << endl;
	for (int i = 0; i < qry_sz; ++i)
	{
		if (local_index.count(qry[i].contig_no) == 0)
		{
			local_index[qry[i].contig_no] = i;
		}
		else
		{
			local_index[qry[i].contig_no] = -i;
		}
	}


	for (int i = 0; i < qry.size(); ++i)
	{
		bool flip1 = 0;
		int contig = qry[i].contig_no;
		if (contig < 0)
		{
			flip1 = 1;
			contig = -contig;
		}
		
		for (int j = 0; j < reads_info->Contig2Ref[contig].size(); j += 2)
		{
			//bool flip2 = 0;
			int read = reads_info->Contig2Ref[contig][j];
			int position = reads_info->Contig2Ref[contig][j + 1];
			if (flip1)
			{
				read = -read;
			}
			
			read_oi[read].push_back(position);
			read_score[read] += contigs_info->contig_sz_vt[abs(contig)];

		}

	}
	
	map<int, vector<int> > sorted_reads;
	for (map<int, int >::iterator it = read_score.begin(); it != read_score.end(); ++it)
	{
		sorted_reads[it->second].push_back(it->first);
	}
	int total_ref = 0;
	//cout << sorted_reads.size() << endl;
	for (map<int, vector<int> >::reverse_iterator rit = sorted_reads.rbegin(); rit != sorted_reads.rend(); ++rit)
	{
		vector<int> ref_vec = rit->second;
		for (int j = 0; j != ref_vec.size(); ++j)
		{
			int ref_no = ref_vec[j];
			total_ref++;
			if (total_ref > 2)
			{
				break;
			}
			vector<Coord_CTG_Cov> ref2,ref=reads_info->RefIndexVec[abs(ref_no)];
			
			align_maps align_maps;
			int min_match = -1, max_match = -1;
			for (int k = 0; k != read_oi[ref_no].size(); ++k)
			{
				if (min_match < 0 || read_oi[ref_no][k] < min_match)
				{
					min_match = read_oi[ref_no][k];
				}
				if (read_oi[ref_no][k]>max_match)
				{
					max_match = read_oi[ref_no][k];
				}
			}
			ref2.clear();
			for (int k = min_match; k <= max_match; ++k)
			{
				ref2.push_back(ref[k]);
			}
			ref = ref2;
			if (ref_no < 0)
			{
				reverse(ref.begin(), ref.end());
				for (int k = 0; k < ref.size(); ++k)
				{
					ref[k].contig_no = -ref[k].contig_no;
				}
			}
			
			int ref_sz = ref.size(), qry_sz = qry.size();
			vector<int> tmp_vec;
			for (int k = 0; k < ref_sz; ++k)
			{
				tmp_vec.push_back(ref[k].contig_no);
			}
			matching_contigs_vec.push_back(tmp_vec);

			struct align_info align_info2 = align_info0;
			align_info2.ref_idx = ref_no;
			align_info2.band_width = abs(ref_sz - qry_sz) + 200;
			align_info2.max_len = max(ref_sz, qry_sz)+10;
			approximate_local_align_banded(contigs_info, qry, ref, &align_maps, &align_info2);

			//approximate_semi_global_align_banded(contigs_info, ref, qry, &align_maps, &align_info2);
			if (align_info2.max_score > align_info2.min_overlap)
			{
				
				int min_match2 = -1, max_match2 = -1,position=0;
				for (int q = 0; q < align_info2.qry_aligned.size(); ++q)
				{

					if (align_info2.qry_aligned[q] == align_info2.ref_aligned[q])
					{
						if (min_match2<0)
						{
							min_match2 = position;
						}
						max_match2 = position;
						
					}
					if (align_info2.qry_aligned[q]!=0)
					{
						position++;
					}
				}
				align_info2.min_match_qry = min_match2;
				align_info2.max_match_qry = max_match2;

				min_match2 = -1, max_match2 = -1, position = 0;
				for (int r = 0; r < align_info2.ref_aligned.size(); ++r)
				{
					if (align_info2.qry_aligned[r] == align_info2.ref_aligned[r])
					{
						if (min_match2<0)
						{
							min_match2 = position;
						}
						max_match2 = position;

					}

					if (align_info2.ref_aligned[r] != 0)
					{
						position++;
					}
				}
				align_info2.min_match_ref = min_match2;
				align_info2.max_match_ref = max_match2;
				align_info_vec.push_back(align_info2);
			}
			
		}
	}

	return 1;
}





#endif