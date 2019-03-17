#ifndef __GRAPH_CONSTRUCTION_H
#define __GRAPH_CONSTRUCTION_H


#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "time.h"
#include "BasicDataStructure.h"

#include "Align.h"
#include "GraphSearch.h"
using namespace std;


// initialize a hashtable
void Init_HT(struct hashtable* ht, size_t ht_sz)
{
	ht->ht_sz = ht_sz;
	ht->store_pos = (struct bucket**)myCalloc(ht_sz, sizeof(struct bucket*));

	for (size_t i = 0; i<ht_sz; ++i)
	{
		ht->store_pos[i] = NULL;
	}
}

void Init_HT2(struct hashtable2* ht, size_t ht_sz)
{
	ht->ht_sz = ht_sz;
	ht->store_pos = (struct bucket2**)myCalloc(ht_sz, sizeof(struct bucket2*));
	for (size_t i = 0; i<ht_sz; ++i)
	{
		ht->store_pos[i] = NULL;
	}
}

void Init_HT3(struct hashtable3* ht, size_t ht_sz)
{
	ht->ht_sz = ht_sz;
	ht->store_pos = (struct bucket3**)myCalloc(ht_sz, sizeof(struct bucket3*));
	for (size_t i = 0; i<ht_sz; ++i)
	{
		ht->store_pos[i] = NULL;
	}
}

void Init_HT4(struct hashtable4* ht, size_t ht_sz)
{
	ht->ht_sz = ht_sz;
	ht->store_pos = (struct bucket4**)myCalloc(ht_sz, sizeof(struct bucket4*));
	for (size_t i = 0; i<ht_sz; ++i)
	{
		ht->store_pos[i] = NULL;
	}
}






//look up for a k-mer in a hashtable, if exists: 1, otherwise: 0. for round 2
bool look_up_in_a_list(uint64_t seq, struct bucket *** ptr)
{
	bool found = 0;
	//struct bucket ** bktptr;


	while ((**ptr) != NULL)
	{
		if ((**ptr)->kmer_t.kmer == seq)
		{
			//	bktptr=ptr;
			break;
		}

		(*ptr) = &((**ptr)->nxt_bucket);
	}
	if ((**ptr) == NULL)
	{
		found = 0;
		//bktptr=ptr;
	}
	else
	{
		found = 1;
	}
	return found;
}

bool look_up_in_a_list2(struct kmer_t2 *seq, struct bucket2 *** ptr)
{
	bool found = 0;
	//struct bucket ** bktptr;

	while ((**ptr) != NULL)
	{
		if (memcmp(&((**ptr)->kmer_t2.kmer), &(seq->kmer), sizeof(uint64_t)* 2) == 0)
		{
			break;
		}

		(*ptr) = &((**ptr)->nxt_bucket);
	}
	if ((**ptr) == NULL)
	{
		found = 0;
	}
	else
	{
		found = 1;
	}
	return found;
}

bool look_up_in_a_list3(struct kmer_t3 *seq, struct bucket3 *** ptr)
{
	bool found = 0;
	//struct bucket ** bktptr;

	while ((**ptr) != NULL)
	{
		if (memcmp(&((**ptr)->kmer_t3.kmer), &(seq->kmer), sizeof(uint64_t)* 3) == 0)
		{
			break;
		}

		(*ptr) = &((**ptr)->nxt_bucket);
	}
	if ((**ptr) == NULL)
	{
		found = 0;
	}
	else
	{
		found = 1;
	}
	return found;
}

bool look_up_in_a_list4(struct kmer_t4 *seq, struct bucket4 *** ptr)
{
	bool found = 0;
	//struct bucket ** bktptr;

	while ((**ptr) != NULL)
	{
		if (memcmp(&((**ptr)->kmer_t4.kmer), &(seq->kmer), sizeof(uint64_t)* 4) == 0)
		{
			break;
		}

		(*ptr) = &((**ptr)->nxt_bucket);
	}
	if ((**ptr) == NULL)
	{
		found = 0;
	}
	else
	{
		found = 1;
	}
	return found;
}


bool look_up_in_a_list_rm(uint64_t seq, struct bucket_rm *** ptr)
{
	bool found = 0;
	//struct bucket ** bktptr;


	while ((**ptr) != NULL)
	{
		if ((**ptr)->kmer_t.kmer == seq)
		{
			//	bktptr=ptr;
			break;
		}

		(*ptr) = &((**ptr)->nxt_bucket);
	}
	if ((**ptr) == NULL)
	{
		found = 0;
		//bktptr=ptr;
	}
	else
	{
		found = 1;
	}
	return found;
}



bool look_up_in_a_list_rm2(struct kmer_t2 *seq, struct bucket_rm2 *** ptr)
{
	bool found = 0;
	//struct bucket ** bktptr;

	while ((**ptr) != NULL)
	{
		if (memcmp(&((**ptr)->kmer_t2.kmer), &(seq->kmer), sizeof(uint64_t)* 2) == 0)
		{
			break;
		}

		(*ptr) = &((**ptr)->nxt_bucket);
	}
	if ((**ptr) == NULL)
	{
		found = 0;
	}
	else
	{
		found = 1;
	}
	return found;
}

bool look_up_in_a_list_rm3(struct kmer_t3 *seq, struct bucket_rm3 *** ptr)
{
	bool found = 0;
	//struct bucket ** bktptr;

	while ((**ptr) != NULL)
	{
		if (memcmp(&((**ptr)->kmer_t3.kmer), &(seq->kmer), sizeof(uint64_t)* 3) == 0)
		{
			break;
		}

		(*ptr) = &((**ptr)->nxt_bucket);
	}
	if ((**ptr) == NULL)
	{
		found = 0;
	}
	else
	{
		found = 1;
	}
	return found;
}

bool look_up_in_a_list_rm4(struct kmer_t4 *seq, struct bucket_rm4 *** ptr)
{
	bool found = 0;
	//struct bucket ** bktptr;

	while ((**ptr) != NULL)
	{
		if (memcmp(&((**ptr)->kmer_t4.kmer), &(seq->kmer), sizeof(uint64_t)* 4) == 0)
		{
			break;
		}

		(*ptr) = &((**ptr)->nxt_bucket);
	}
	if ((**ptr) == NULL)
	{
		found = 0;
	}
	else
	{
		found = 1;
	}
	return found;
}



// compare 2 arrays.
int uint64_t_cmp(uint64_t* A, uint64_t* B, int Kmer_arr_sz)
{
	int flag = 0;
	for (int jj = 0; jj<Kmer_arr_sz; ++jj)
	{
		if (A[jj]>B[jj])
		{
			flag = 1;
			break;
		}
		if (A[jj]<B[jj])
		{
			flag = -1;
			break;
		}
		if (A[jj] == B[jj])
		{
			continue;
		}
	}
	return flag;


}


//look up for a k-mer in a hashtable, if exists: 1, otherwise: 0. for round 1
bool look_up_in_a_list_r1(uint64_t seq, struct bucket_r1 *** ptr)
{
	bool found = 0;

	while ((**ptr) != NULL)
	{
		if ((**ptr)->kmer_t.kmer == seq)
		{
			//	bktptr=ptr;	
			break;
		}

		(*ptr) = &((**ptr)->nxt_bucket);
	}
	if ((**ptr) == NULL)
	{
		found = 0;
		//bktptr=ptr;
	}
	else
	{
		found = 1;
	}
	return found;
}

bool look_up_in_a_list2_r1(struct kmer_t2 *seq, struct bucket2_r1 *** ptr)
{
	bool found = 0;

	while ((**ptr) != NULL)
	{
		if (memcmp(&((**ptr)->kmer_t2.kmer), &(seq->kmer), sizeof(uint64_t)* 2) == 0)
		{
			break;
		}

		(*ptr) = &((**ptr)->nxt_bucket);
	}
	if ((**ptr) == NULL)
	{
		found = 0;
	}
	else
	{
		found = 1;
	}
	return found;
}

bool look_up_in_a_list3_r1(struct kmer_t3 *seq, struct bucket3_r1 *** ptr)
{
	bool found = 0;

	while ((**ptr) != NULL)
	{
		if (memcmp(&((**ptr)->kmer_t3.kmer), &(seq->kmer), sizeof(uint64_t)* 3) == 0)
		{
			break;
		}

		(*ptr) = &((**ptr)->nxt_bucket);
	}
	if ((**ptr) == NULL)
	{
		found = 0;
	}
	else
	{
		found = 1;
	}
	return found;
}

bool look_up_in_a_list4_r1(struct kmer_t4 *seq, struct bucket4_r1 *** ptr)
{
	bool found = 0;

	while ((**ptr) != NULL)
	{
		if (memcmp(&((**ptr)->kmer_t4.kmer), &(seq->kmer), sizeof(uint64_t)* 4) == 0)
		{
			break;
		}

		(*ptr) = &((**ptr)->nxt_bucket);
	}
	if ((**ptr) == NULL)
	{
		found = 0;
	}
	else
	{
		found = 1;
	}
	return found;
}


void SwitchBuckets(hashtable *ht)
{
	size_t ht_sz;
	
	ht_sz = ht->ht_sz;
	bucket_r1 *store_pos_o, *store_pos_t;
	bucket *store_pos_n;
	bucket **bktp2p;
	for (size_t i = 0; i<ht_sz; ++i)
	{
		
		bktp2p = &(ht->store_pos[i]);
		store_pos_o = (bucket_r1*)ht->store_pos[i];
		while (store_pos_o != NULL)
		{
			if (store_pos_o->kmer_t.cov==1)
			{
				store_pos_n = (bucket*)myMalloc(sizeof(struct bucket));
				memset(store_pos_n, 0, sizeof(struct bucket));
				store_pos_n->kmer_t = store_pos_o->kmer_t;
				store_pos_n->kmer_info.cov1 = store_pos_o->kmer_t.cov;
				*bktp2p = store_pos_n;
				bktp2p = &(store_pos_n->nxt_bucket);
			}

			store_pos_t = store_pos_o;
			store_pos_o = store_pos_o->nxt_bucket;
			free(store_pos_t);
		}
		*bktp2p = NULL;
	}

}


//reference graph




void Contig_Kmer_Index(struct ref_read_t *read, struct hashtable *ht, int K_size, int * bucket_count)
{
	int readLen = read->readLen;
	int OverlappingKmers = readLen - K_size + 1;
	if (OverlappingKmers <= 0)
	{
		return;
	}
	int Read_arr_sz = readLen / 32 + 1;
	int rem = readLen % 32;
	if (rem == 0)
	{
		Read_arr_sz--;
	}
	int tot_bits = Read_arr_sz * 64;

	size_t ht_sz = ht->ht_sz;
	bool flip, found;
	size_t hash_idx;
	uint64_t seq, f_seq, hv;
	bucket **bktptr;


	for (int j = 0; j<OverlappingKmers; j++)
	{
		get_sub_arr(read->read_bits, read->readLen, j, K_size, &seq);
		f_seq = get_rev_comp_seq(seq, K_size);
		flip = 0;
		if (seq>f_seq)
		{
			uint64_t t = seq;
			seq = f_seq;
			f_seq = t;
			flip = 1;
		}

		hv = MurmurHash64A(&seq, sizeof(seq), 0);

		hash_idx = (size_t)(hv%ht_sz);

		bktptr = &(ht->store_pos[hash_idx]);
		if (ht->round == 1)
		{
			found = look_up_in_a_list_r1(seq, (bucket_r1***)&bktptr);
		}
		else
		{
			found = look_up_in_a_list(seq, &bktptr);
		}


		if (ht->round == 1)
		{
			if (found == 0)
			{
				*(bktptr) = (struct bucket*)myMalloc(sizeof(struct bucket_r1));

				memset(*bktptr, 0, sizeof(struct bucket_r1));

				((struct bucket_r1*) *(bktptr))->nxt_bucket = NULL;

				((struct bucket_r1*) *(bktptr))->kmer_t.kmer = seq;

				((struct bucket_r1*) *(bktptr))->kmer_t.cov++;


				(*bucket_count)++;
			}
			else
			{
				if (((struct bucket_r1*) *(bktptr))->kmer_t.cov <= 1)
				{
					((struct bucket_r1*) *(bktptr))->kmer_t.cov++;
				}
			}


		}

		if (ht->round == 2 && found == 1)
		{

			(*(bktptr))->kmer_info.contig_no = read->contig_no;

			(*(bktptr))->kmer_info.cod = j;

			(*(bktptr))->kmer_info.flip = flip;
			(*bucket_count)++;
		}

	}

}


void ReadCompression(struct ref_read_t *read, struct hashtable *ht, contigs_info* contigs_info, int K_size, struct LongReadContigIndex *LongReadContigIndex)
{
	//ofstream out_positions_debug("CenterEstimates.txt",std::ios_base::app | std::ios_base::out);
	LongReadContigIndex->LR2CTG.clear();
	LongReadContigIndex->CTG2LR.clear();
	LongReadContigIndex->CTG2Coords.clear();
	int readLen = read->readLen;
	int OverlappingKmers = readLen - K_size + 1;
	if (OverlappingKmers <= 1)
	{
		return;
	}
	int Read_arr_sz = readLen / 32 + 1;
	int rem = readLen % 32;
	if (rem == 0)
	{
		Read_arr_sz--;
	}
	int tot_bits = Read_arr_sz * 64;

	size_t ht_sz = ht->ht_sz;
	bool flip, found;
	size_t hash_idx;
	uint64_t seq, f_seq, hv;
	bucket **bktptr;


	int first_found = -1;

	for (int j = 0; j<OverlappingKmers; j++)
	{
		get_sub_arr(read->read_bits, read->readLen, j, K_size, &seq);
		f_seq = get_rev_comp_seq(seq, K_size);
		flip = 0;
		if (seq>f_seq)
		{
			uint64_t t = seq;
			seq = f_seq;
			f_seq = t;
			flip = 1;
		}

		hv = MurmurHash64A(&seq, sizeof(seq), 0);

		hash_idx = (size_t)(hv%ht_sz);

		bktptr = &(ht->store_pos[hash_idx]);

		found = look_up_in_a_list(seq, &bktptr);

		if (found == 0)
		{

		}
		else
		{
			if ((*(bktptr))->kmer_info.cov1 > 1)
			{
				continue;
			}
			else
			{
				LongReadContigIndex->nMatches++;
				uint32_t contig_no = (*(bktptr))->kmer_info.contig_no;
				int cod = (*(bktptr))->kmer_info.cod;

				if (cod > readLen && cod < (contigs_info->contig_sz_vt[contig_no] - readLen))
				{
					//continue;//contained hit
				}
				LongReadContigIndex->LR2CTG[j].contig_no = contig_no;
				LongReadContigIndex->LR2CTG[j].pos = cod;
				LongReadContigIndex->LR2CTG[j].flip = flip ^ (*(bktptr))->kmer_info.flip;
				if (LongReadContigIndex->LR2CTG[j].flip == 0)
				{
					int offset = j - (*bktptr)->kmer_info.cod + contigs_info->contig_sz_vt[contig_no] / 2;
					LongReadContigIndex->CTG2LR[(*(bktptr))->kmer_info.contig_no].push_back(j);
					LongReadContigIndex->CTG2LR_2[(*(bktptr))->kmer_info.contig_no].push_back(offset);
					if (LongReadContigIndex->FastMap == 0)
					{
						LongReadContigIndex->CTG2Coords[(*(bktptr))->kmer_info.contig_no][j] = cod;
						LongReadContigIndex->CTG2Offsets[(*(bktptr))->kmer_info.contig_no][j] = offset;

					}
					
				}
				else
				{
					int offset = j + (*bktptr)->kmer_info.cod - contigs_info->contig_sz_vt[contig_no] / 2 + (K_size + 1) / 2;
					LongReadContigIndex->CTG2LR[-(*(bktptr))->kmer_info.contig_no].push_back(j);
					LongReadContigIndex->CTG2LR_2[-(*(bktptr))->kmer_info.contig_no].push_back(offset);
					if (LongReadContigIndex->FastMap == 0)
					{
						LongReadContigIndex->CTG2Coords[-(*(bktptr))->kmer_info.contig_no][j] =cod+K_size;
						LongReadContigIndex->CTG2Offsets[-(*(bktptr))->kmer_info.contig_no][j] = offset;
					}
					
				}
				

			}

		}

	}
	//some ugly post processing.
	bool debug = 0;
	map<int, vector<int> >::iterator tmp_it, tmp_it2;
	tmp_it2 = LongReadContigIndex->CTG2LR_2.begin();
	int best_ctg = 0, max_cov = 0, ctg_center = 0, coord1 = 0, coord2 = 0;
	for (tmp_it = LongReadContigIndex->CTG2LR.begin(); tmp_it != LongReadContigIndex->CTG2LR.end(); ++tmp_it)
	{

		sort(tmp_it->second.begin(), tmp_it->second.end());
		sort(tmp_it2->second.begin(), tmp_it2->second.end());

		//	LongReadContigIndex->layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
		//	LongReadContigIndex->layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
		if (contigs_info->contig_sz_vt[abs(tmp_it->first)] < LongReadContigIndex->BlockSize || tmp_it2->second.size()<60)
		{
			int center_est = tmp_it->second[tmp_it->second.size() / 2];
			if (LongReadContigIndex->layout.count(center_est) == 0)
			{
				LongReadContigIndex->layout[center_est].ctg_no = tmp_it->first;//ctg_no
				LongReadContigIndex->layout[center_est].cov = tmp_it->second.size();//cov
				LongReadContigIndex->layout[center_est].coord2 = tmp_it2->second[tmp_it2->second.size() / 2];//coord2
			}
		}
		else
		{
			//resolving repeats
			int64_t center_sum = 0, cnt = 0;
			vector<int> center_vec;
			center_vec.reserve(2000);
			bool resolved = 0;
			for (int p = 10; p<=tmp_it2->second.size()-10; ++p)
			{
				if ((abs(tmp_it2->second[p + 1] - tmp_it2->second[p])<contigs_info->contig_sz_vt[abs(tmp_it->first)] / 2) && (p != tmp_it2->second.size() - 10))
				{
					center_sum += tmp_it2->second[p + 1];
					center_vec.push_back(tmp_it2->second[p + 1]);
					cnt++;
				}
				else
				{
					if (cnt>0)
					{
						int center_est = center_vec[center_vec.size()/2];// (int)(center_sum / cnt);
						
						if (cnt>(0.002*contigs_info->contig_sz_vt[abs(tmp_it->first)]) || cnt>20)
						{
							cnt += 20;
							resolved = 1;
							if (LongReadContigIndex->layout.count(center_est) == 0)
							{
								
								LongReadContigIndex->layout[center_est].ctg_no=tmp_it->first;//ctg_no
								LongReadContigIndex->layout[center_est].cov=cnt;//cov
								LongReadContigIndex->layout[center_est].coord2=center_est;//coord2
							}
							else
							{
								
								//resolve conflicts (same center coord)
								double ratio1 = (double)cnt / (double)(contigs_info->contig_sz_vt[abs(tmp_it->first)] + 1);
								double ratio2 = (double)LongReadContigIndex->layout[center_est].cov / (double)(contigs_info->contig_sz_vt[abs(LongReadContigIndex->layout[center_est].ctg_no)] + 1);
								if (ratio1 > ratio2)
								{
									
									LongReadContigIndex->layout[center_est].ctg_no=(tmp_it->first);//ctg_no
									LongReadContigIndex->layout[center_est].cov=(cnt);//cov
									LongReadContigIndex->layout[center_est].coord2=(center_est);//coord2
								}
							}
							//out_positions_debug << "ctg: " << tmp_it->first << " cov: " << cnt << " coord: " << center_est << endl;
							


						}
						else
						{
							//cout <<  "weak hit. " << "ctg: " << tmp_it->first << " cov: " << cnt << " coord: " << center_est << endl;;
							//out_positions_debug << "weak hit. " << "ctg: " << tmp_it->first << " cov: " << cnt << " coord: " << center_est << endl;
						}
					}
					else
					{
						/*
						out_positions_debug << "zero hit. " << "ctg: " << tmp_it->first << endl;
						for (vector<int>::iterator vit = tmp_it2->second.begin(); vit != tmp_it2->second.end(); ++vit)
						{
							out_positions_debug << *vit << ", ";
						}
						out_positions_debug << endl;
						*/
					}
					center_vec.clear();
					center_sum = 0;
					cnt = 0;
				}

			}
			if (resolved == 0)
			{
				int center_est = tmp_it->second[tmp_it->second.size() / 2];
				if (LongReadContigIndex->layout.count(center_est) == 0)
				{
					LongReadContigIndex->layout[center_est].ctg_no = tmp_it->first;//ctg_no
					LongReadContigIndex->layout[center_est].cov = tmp_it->second.size();//cov
					LongReadContigIndex->layout[center_est].coord2 = tmp_it2->second[tmp_it2->second.size() / 2];//coord2
				}
			}

		}
		
		if (tmp_it->second.size() > max_cov)
		{
			max_cov = tmp_it->second.size();
			best_ctg = tmp_it->first;
			ctg_center = tmp_it2->second[tmp_it2->second.size() / 2];
			coord1 = tmp_it->second[tmp_it->second.size() / 2];
			coord2 = tmp_it2->second[tmp_it2->second.size() / 2];
		}
		tmp_it2++;


		
	}
	
	//the greedy cleaning part: assume best ctg is correct
	if (LongReadContigIndex->FastMap == 1)
	{
		
		if ((ctg_center + contigs_info->contig_sz_vt[abs(best_ctg)] / 2 > readLen) && (ctg_center - contigs_info->contig_sz_vt[abs(best_ctg)] / 2 < 0))
		{
			LongReadContigIndex->layout.clear();
			LongReadContigIndex->layout[coord1].ctg_no = (best_ctg);
			LongReadContigIndex->layout[coord1].cov = (max_cov);
			LongReadContigIndex->layout[coord1].coord2 = (coord2);


		}
	}

}



void SaveContigKmerIndex(hashtable *ht, string &fname)
{
	ofstream  o_ht_idx, o_ht_content;
	//File *ht_content;
	map<int, int> cov_hist, edge_cov_hist;

	string ht_idx_name = fname + "HT_idx.txt", ht_content_name = fname + "HT_content";
	o_ht_idx.open(ht_idx_name.c_str(), ios_base::out | ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(), ios_base::out | ios_base::binary);
	o_ht_idx << "Hashtable size: " << endl << ht->ht_sz << endl;

	bucket * bktptr = NULL;

	for (size_t i = 0; i<ht->ht_sz; ++i)
	{
		size_t list_sz = 0;
		bktptr = ht->store_pos[i];
		while (bktptr != NULL)
		{
			o_ht_content.write((char*)bktptr, sizeof(struct bucket));

			int cov = bktptr->kmer_info.cov1;

			cov_hist[cov]++;
			bktptr = bktptr->nxt_bucket;
			list_sz++;
		}
		o_ht_idx << list_sz << endl;
	}

}

//load the saved information
void LoadContigKmerIndex(hashtable *ht, string &fname)
{
	string ht_idx_name = fname + "HT_idx.txt", ht_content_name = fname + "HT_content";
	ifstream in_ht_idx(ht_idx_name.c_str(), ios_base::in | ios_base::binary), in_ht_content(ht_content_name.c_str(), ios_base::in | ios_base::binary);

	size_t ht_sz;
	string s;
	getline(in_ht_idx, s);
	getline(in_ht_idx, s);
	//ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT(ht, ht_sz);
	for (size_t i = 0; i<ht_sz; ++i)
	{
		int list_sz;
		getline(in_ht_idx, s);
		if (s[s.size() - 1] == '\r' || s[s.size() - 1] == '\n')
		{
			s.resize(s.size() - 1);
		}
		list_sz = atoi(s.c_str());
		struct bucket **bktp2p = &(ht->store_pos[i]);
		*bktp2p = NULL;
		for (int j = 0; j<list_sz; ++j)
		{
			*bktp2p = (struct bucket*)myMalloc(sizeof(struct bucket));
			in_ht_content.read((char*)(*bktp2p), sizeof(struct bucket));

			//(*bktp2p)->kmer_info.used = 0;
			(*bktp2p)->nxt_bucket = NULL;

			bktp2p = &((*bktp2p)->nxt_bucket);
		}
	}

}


void FreeSparseKmerGraph(struct hashtable *ht)
{

	bucket * bktptr, *n_bktptr;//=NULL;
	edge_node *edgeptr, *n_edgeptr;
	for (size_t i = 0; i<ht->ht_sz; ++i)
	{
		size_t list_sz = 0;
		bktptr = ht->store_pos[i];
		while (bktptr != NULL)
		{
			n_bktptr = bktptr->nxt_bucket;
			

			free(bktptr);
			bktptr = n_bktptr;
		}

	}
	free(ht->store_pos);


}







void CorrectPacBioRead(contigs_info *contigs_info,reads_info *reads_info, vector<Coord_CTG_Cov> &TempIndex)
{
	map<int, int> local_index;
	int vec_sz = TempIndex.size();
	for (int i = 0; i < vec_sz; ++i)
	{
		local_index[TempIndex[i].contig_no]=i+1;
	}
	int trim_head = -1,trim_rear=vec_sz;
	for (int i = 0; i < vec_sz; ++i)
	{
		int ctg = -(TempIndex[i].contig_no);
		if (local_index.count(ctg) &&local_index[ctg]>i && contigs_info->contig_sz_vt[abs(ctg)]>reads_info->MinOverlap)
		{
			trim_head = i;
		}
	}
	for (int i = vec_sz - 1; i >= 0; --i)
	{
		int ctg = -(TempIndex[i].contig_no);
		if (local_index.count(ctg) && local_index[ctg]<= i && contigs_info->contig_sz_vt[abs(ctg)]>reads_info->MinOverlap)
		{
			trim_rear = i;
		}
	}
	vector<Coord_CTG_Cov> TempIndex2;

	if (trim_head != -1 )
	{
		for (int i = 0; i < vec_sz; ++i)
		{
			//cout << TempIndex[i].contig_no << ", ";
		}
		//cout << endl;
		if (trim_head < (vec_sz - 1 - trim_rear))
		{
			for (int i = trim_head + 1; i < vec_sz; ++i)
			{
				TempIndex2.push_back(TempIndex[i]);
			}
			TempIndex = TempIndex2;
		}
		if (trim_head>(vec_sz - 1 - trim_rear))
		{
			for (int i = 0; i < trim_rear; ++i)
			{
				TempIndex2.push_back(TempIndex[i]);
			}
			TempIndex = TempIndex2;
			
		}
		if (trim_head == (vec_sz - 1 - trim_rear))
		{
			for (int i = trim_head + 1; i < trim_rear; ++i)
			{
				TempIndex2.push_back(TempIndex[i]);
			}
			TempIndex = TempIndex2;

		}
		//cout << trim_head <<" "<< trim_rear << ": ";
		for (int i = 0; i < TempIndex.size(); ++i)
		{
			//cout << TempIndex[i].contig_no << ", ";
		}
		//cout << endl;
	}

}

void LoadLongReadIndex(vector<string> filename_vec, reads_info *reads_info, contigs_info *contigs_info)//vector<vector<Coord_CTG_Cov> > &LongReadIndexVec, vector< vector<int> > &Contig2LongRead)
{
	reads_info->LongReadIndexVec.clear();
	reads_info->Contig2LongRead.clear();
	reads_info->Contig2LongRead.resize(contigs_info->total_contigs + 1);
	
	string tag, n_tag, tmp1, tmp2,tmp3,str;
	reads_info->tag_vec.clear();
	reads_info->LenVec.clear();
	int coord, contig, cov,coord2;
	int n_reads = 0;
	int readlen = 0;
	reads_info->MaxCompressedReadLen = 0;
	reads_info->tag_vec.push_back(tag);
	reads_info->LenVec.push_back(readlen);
	vector<Coord_CTG_Cov> TempIndex, FinalIndex;

	reads_info->LongReadIndexVec.push_back(FinalIndex);
	FinalIndex.reserve(1000);
	TempIndex.reserve(1000);
	Coord_CTG_Cov Coord_CTG_Cov;
	uint64_t len_sum=0, n_compressed_reads = 0;
	for (int f = 0; f < filename_vec.size(); ++f)
	{
		string filename = filename_vec[f];
		cout << "Loading file: " << filename << endl;
		ifstream LongReadIndex_in(filename.c_str());
		bool read_success = 1;
		int Len = 0;
		while (read_success = get_a_contig_path(LongReadIndex_in, tag, TempIndex, Len, reads_info->KmerCovTh, n_tag))
		{
			//if (tag==">m130605_000141_42207_c100515142550000001823076608221372_s1_p0/118764/0_11490")
			//cout << "";
			if (reads_info->Clean)
			{
				if (TempIndex.size() == 0)
				{
					continue;
				}
				if ((Len > 0 && Len < reads_info->MinLen))//||(abs( TempIndex[0].coord2 - TempIndex[TempIndex.size()-1].coord2)<reads_info->MinLen)
				{
					continue;
				}
			}
			
			reads_info->tag_vec.push_back(tag);

			reads_info->LenVec.push_back(Len);
			n_reads++;
			if (reads_info->TotalReads>0 && n_reads > reads_info->TotalReads)
			{
				break;
			}
			FinalIndex.clear();


			//basic cleaning
			
			for (int t = 0; t < TempIndex.size(); ++t)
			{
				int contig = TempIndex[t].contig_no;
				int ctg_sz = contigs_info->contig_sz_vt[abs(contig)];
				int beg = TempIndex[t].coord2 - ctg_sz / 2;
				int end = TempIndex[t].coord2 + ctg_sz / 2;
				int contained_sz = ctg_sz;
				if (beg < 0)
				{
					contained_sz += beg;
				}
				if (Len>0 && end > Len)
				{
					contained_sz -= (end - Len);
				}

				if ((!reads_info->Clean)||TempIndex[t].cov >= ((double)ctg_sz)*reads_info->AdaptiveTh)
				{
					FinalIndex.push_back(TempIndex[t]);

					if (reads_info->Redundancy<0.0 || (!reads_info->Clean))//otherwisedo it in the advanced version
					if ((reads_info->Contig2LongRead[abs(contig)].size() == 0) || (reads_info->Contig2LongRead[abs(contig)].back() != n_reads))
					{
						reads_info->Contig2LongRead[abs(contig)].push_back(n_reads);
					}
				}

			}

			if (reads_info->Clean&&reads_info->Redundancy>0.0)
			{
				TempIndex = FinalIndex;
				FinalIndex.clear();

				//advanced cleaning
				map<double, vector<int> > ranked_contigs;
				vector <bool> selected_positions;
				selected_positions.resize(TempIndex.size());
				for (int t = 0; t < TempIndex.size(); ++t)
				{
					int contig = TempIndex[t].contig_no;
					int ctg_sz = contigs_info->contig_sz_vt[abs(contig)];
					int beg = TempIndex[t].coord2 - ctg_sz / 2;
					int end = TempIndex[t].coord2 + ctg_sz / 2;
					int contained_sz = ctg_sz;
					if (beg < 0)
					{
						contained_sz += beg;
					}
					if (Len>0 && end>Len)
					{
						contained_sz -= (end - Len);
					}

					ranked_contigs[(double)TempIndex[t].cov / (double)contained_sz].push_back(t);
				

					
				}
				map<double, vector<int> >::reverse_iterator rit;
				int covered_len = 0;
				for (rit = ranked_contigs.rbegin(); rit != ranked_contigs.rend(); ++rit)
				{
					for (int cc = 0; cc < rit->second.size(); ++cc)
					{
						int t = rit->second[cc];
						int contig = TempIndex[t].contig_no;
						int ctg_sz = contigs_info->contig_sz_vt[abs(contig)];
						int beg = TempIndex[t].coord2 - ctg_sz / 2;
						int end = TempIndex[t].coord2 + ctg_sz / 2;
						int contained_sz = ctg_sz;
						if (beg < 0)
						{
							contained_sz += beg;
						}
						if (Len>0 && end>Len)
						{
							contained_sz -= (end - Len);
						}
						covered_len += contained_sz;
						if (covered_len > reads_info->Redundancy*((double)Len))
						{
							break;
						}
						selected_positions[t] = 1;
					
					}
					
				}

				for (int t = 0; t < TempIndex.size(); ++t)
				{
					if (selected_positions[t])
					{
						FinalIndex.push_back(TempIndex[t]);
						int contig = TempIndex[t].contig_no;
						if ((reads_info->Contig2LongRead[abs(contig)].size() == 0) || (reads_info->Contig2LongRead[abs(contig)].back() != n_reads))
						{
							reads_info->Contig2LongRead[abs(contig)].push_back(n_reads);
						}
					}
					else
					{
						cout << "";
					}
				}
			}
			
			if (FinalIndex.size() > reads_info->MaxCompressedReadLen)
			{
				reads_info->MaxCompressedReadLen = FinalIndex.size();
			}
			reads_info->LongReadIndexVec.push_back(FinalIndex);
			if (FinalIndex.size() > 1)
			{
				len_sum += FinalIndex.size();
				n_compressed_reads++;
			}
			
		}
		

		
	}
	
	cout << n_reads << " reads loaded."<<endl;
	cout << "Average size: " << (len_sum/n_compressed_reads) << endl;
	ofstream o_read_contig_idx("selected_reads.txt");
	int total_reads=reads_info->LongReadIndexVec.size();
	for (int i = 1; i < total_reads; ++i)
	{
		o_read_contig_idx << ">" << i << endl;

		for (int v = 0; v != reads_info->LongReadIndexVec[i].size(); ++v)
		{
			struct Coord_CTG_Cov Coord_CTG_Cov_tmp = reads_info->LongReadIndexVec[i][v];
			
			o_read_contig_idx << Coord_CTG_Cov_tmp.coord << ", ";
			o_read_contig_idx << Coord_CTG_Cov_tmp.contig_no << ", ";
			o_read_contig_idx << Coord_CTG_Cov_tmp.cov<<", ";
			o_read_contig_idx << "[" << Coord_CTG_Cov_tmp.coord2 - contigs_info->contig_sz_vt[abs(Coord_CTG_Cov_tmp.contig_no)]/2 << ", ";
			o_read_contig_idx << Coord_CTG_Cov_tmp.coord2 + contigs_info->contig_sz_vt[abs(Coord_CTG_Cov_tmp.contig_no)] / 2 << "], ";
			o_read_contig_idx << (double)Coord_CTG_Cov_tmp.cov / (double)(contigs_info->contig_sz_vt[abs(Coord_CTG_Cov_tmp.contig_no)] + 1);

			o_read_contig_idx << endl;

		}
	}
	
	

	reads_info->LenVec.resize(reads_info->LongReadIndexVec.size());
	reads_info->LengthRank.clear();
	reads_info->LengthRank.reserve(reads_info->LongReadIndexVec.size() + 10);
	if (0)
	for (int i = 0; i < reads_info->LongReadIndexVec.size(); ++i)
	{
		int Len = 0;

		int vec_sz = reads_info->LongReadIndexVec[i].size();
		for (int j = 0; j < vec_sz; ++j)
		{
			Len = reads_info->LongReadIndexVec[i][vec_sz - 1].coord - reads_info->LongReadIndexVec[i][0].coord;
		}

		reads_info->LenVec[i] = Len;
	}

	vector<int>::iterator vit = reads_info->LenVec.begin();
	int numReads = reads_info->LenVec.size() - 1;
	vit++;
	for (size_t i = 0; i < numReads; ++i)
	{
		reads_info->LengthRank.push_back(vit);
		vit++;
	}

	sort(reads_info->LengthRank.begin(), reads_info->LengthRank.end(), it_cmp);



	contigs_info->LengthRank.clear();
	contigs_info->LengthRank.reserve(reads_info->LongReadIndexVec.size() + 10);


	vit = contigs_info->contig_sz_vt.begin();
	int numContigs = contigs_info->contig_sz_vt.size()-1;//double check
	vit++;
	for (size_t i = 0; i < numContigs; ++i)
	{
		contigs_info->LengthRank.push_back(vit);
		vit++;
	}

	sort(contigs_info->LengthRank.begin(), contigs_info->LengthRank.end(), it_cmp);



	reads_info->contained.clear();
	reads_info->contained.resize(numReads + 10);

	reads_info->LeftOverlaps.resize(numReads + 10);
	reads_info->RightOverlaps.resize(numReads + 10);
	reads_info->LeftOverlapsBestWithTies.resize(numReads + 10);
	reads_info->RightOverlapsBestWithTies.resize(numReads + 10);
	reads_info->used_vt_left.resize(numReads + 10);
	reads_info->used_vt_right.resize(numReads + 10);
	reads_info->LeftOverlapsBest.resize(numReads + 10);
	reads_info->RightOverlapsBest.resize(numReads + 10);
	reads_info->LeftBestOverlapsTemp.resize(numReads + 10);
	reads_info->RightBestOverlapsTemp.resize(numReads + 10);
	reads_info->both_stand_used.resize(numReads + 10);
	reads_info->chimeric.clear();
	reads_info->chimeric.resize(numReads + 10);

}

void LoadLongReadIndexWithMerging(vector<string> filename_vec, reads_info *reads_info, contigs_info *contigs_info)//vector<vector<Coord_CTG_Cov> > &LongReadIndexVec, vector< vector<int> > &Contig2LongRead)
{
	reads_info->LongReadIndexVec.clear();
	reads_info->Contig2LongRead.clear();
	reads_info->Contig2LongRead.resize(contigs_info->total_contigs + 1);
	reads_info->MaxCompressedReadLen = 0;
	string tag, n_tag, tmp1, tmp2, tmp3, str;
	int coord, contig, cov, coord2;
	int n_reads = 0,n_bad_reads=0;
	vector<struct Coord_CTG_Cov> TempIndex;
	vector<int> temp_vec, temp_vec_rc;
	TempIndex.reserve(1000);
	temp_vec.reserve(1000); 
	temp_vec_rc.reserve(1000);
	Coord_CTG_Cov Coord_CTG_Cov;
	map<int, list<vector<int> > > sorted_index;
	map<vector<int>, int> read_map;

	uint64_t len_sum = 0, n_compressed_reads = 0;
	for (int f = 0; f < filename_vec.size(); ++f)
	{
		string filename = filename_vec[f];
		cout << "Loading file: " << filename << endl;
		ifstream LongReadIndex_in(filename.c_str());

		bool read_success = 1;
		int Len = 0;
		while (read_success = get_a_contig_path(LongReadIndex_in, tag, TempIndex, Len, reads_info->KmerCovTh, n_tag))
		{

			n_reads++;

			if (reads_info->TotalReads>0 && n_reads>reads_info->TotalReads)
			{
				break;
			}

			temp_vec.clear();
			temp_vec_rc.clear();

			bool bad_read = 0;
			if (TempIndex.size()>0)
			for (int l = 0; l <TempIndex.size() - 1; ++l)
			{

				if (TempIndex[l].coord2 - TempIndex[l + 1].coord2 > 100)
				{
					bad_read = 1;
					n_bad_reads++;
					//cout << TempIndex[l].coord2 << endl;
					//cout << TempIndex[l+1].coord2 << endl;

					break;
				}
			}

				
			if (bad_read)
			{
				//cout << "bad read detected." << endl;
				TempIndex.clear();
			}

			for (int v = 0; v < TempIndex.size(); ++v)
			{
				temp_vec.push_back(TempIndex[v].contig_no);
				temp_vec_rc.push_back(-TempIndex[TempIndex.size() - 1 - v].contig_no);
			}
			bool take_rc = 0;

			for (int v = 0; v < TempIndex.size(); ++v)
			{
				if (temp_vec[v] < temp_vec_rc[v])
				{
					take_rc = 0;
					break;
				}
				if (temp_vec[v] > temp_vec_rc[v])
				{
					take_rc = 1;
					break;
				}
			}
			if (take_rc)
			{
				temp_vec = temp_vec_rc;
			}
			if (temp_vec.size() > 0)
			{
				read_map[temp_vec] ++;

			}
			TempIndex.clear();

			
		}
		cout << n_bad_reads << " bad reads skipped." << endl;
		cout << n_reads << " reads loaded." << endl;
	
	}
	//reads_info->LongReadIndexVec.push_back(TempIndex);
	reads_info->LongReadIndexVec.clear();
	reads_info->LongReadIndexVec.resize(read_map.size()+1);
	ofstream o_read_contig_idx("selected_reads.txt");
	ofstream o_sorted_read_contig_idx("sorted_selected_reads.txt");


	map<vector<int>, int>::iterator read_map_it;
	int n_selected_read = 0;
	for (read_map_it = read_map.begin(); read_map_it != read_map.end(); ++read_map_it)
	{
		sorted_index[read_map_it->second].push_back(read_map_it->first);
		int temp_len1 = 0,temp_len2=0;
		int path_sz=read_map_it->first.size();
		temp_len1 += contigs_info->contig_sz_vt[abs(read_map_it->first[0])];

		temp_len1 += contigs_info->contig_sz_vt[abs(read_map_it->first[path_sz-1])];

	
		if (read_map_it->second <= reads_info->PathCovTh)  //&&temp_len1<500
		{
			continue;
		}

		n_selected_read++;
		o_read_contig_idx << ">" << n_selected_read << endl;
		temp_vec = read_map_it->first;
		vector<struct Coord_CTG_Cov> temp_index;
		temp_index.reserve(1000);
		//cout << temp_vec.size() << endl;
		if (temp_vec.size() > 0)
		{
			reads_info->Selected_Overlaps[temp_vec].clear();
			//cout << "cleared." << endl;
		}
		for (int v = 0; v != temp_vec.size(); ++v)
		{
			struct Coord_CTG_Cov Coord_CTG_Cov_tmp;
			Coord_CTG_Cov_tmp.coord = v;
			Coord_CTG_Cov_tmp.contig_no = temp_vec[v];
			Coord_CTG_Cov_tmp.cov = contigs_info->contig_sz_vt[abs(temp_vec[v])];
			Coord_CTG_Cov_tmp.coord2 = v;
			o_read_contig_idx << Coord_CTG_Cov_tmp.coord << ", ";
			o_read_contig_idx << Coord_CTG_Cov_tmp.contig_no << ", ";
			o_read_contig_idx << Coord_CTG_Cov_tmp.cov<<", ";
			o_read_contig_idx << Coord_CTG_Cov_tmp.coord2;
			o_read_contig_idx << endl;
			temp_index.push_back(Coord_CTG_Cov_tmp);
		}

		if (temp_index.size() > reads_info->MaxCompressedReadLen)
		{
			reads_info->MaxCompressedReadLen = temp_index.size();
		}

		reads_info->LongReadIndexVec[n_selected_read] = temp_index;
		for (int v = 0; v != temp_index.size(); ++v)
		{
			contig = temp_index[v].contig_no;
			if ((reads_info->Contig2LongRead[abs(contig)].size() == 0) || (reads_info->Contig2LongRead[abs(contig)].back() != n_selected_read))
			{
				reads_info->Contig2LongRead[abs(contig)].push_back(n_selected_read);
			}
		}



		if (temp_index.size() > 1)
		{
			len_sum += temp_index.size();
			n_compressed_reads++;
		}


	}

	cout << n_selected_read << " selected reads." << endl;



	cout << n_reads << " reads loaded." << endl;
	cout << "Average size: " << (len_sum / n_compressed_reads) << endl;

	if (reads_info->Debug)
	{
		map<int, list< vector<int> > >::reverse_iterator sorted_read_map_it;
		int count = 0;
		for (sorted_read_map_it = sorted_index.rbegin(); sorted_read_map_it != sorted_index.rend(); ++sorted_read_map_it)
		{
			list< vector<int> > temp_list = sorted_read_map_it->second;
			list< vector<int> >::iterator lit;

			for (lit = temp_list.begin(); lit != temp_list.end(); ++lit)
			{
				vector<int> temp_vec = *lit;
				count++;
				o_sorted_read_contig_idx << ">" << count << "_cov_" << sorted_read_map_it->first << endl;
				for (int v = 0; v != temp_vec.size(); ++v)
				{

					struct Coord_CTG_Cov Coord_CTG_Cov_tmp;
					Coord_CTG_Cov_tmp.coord = v;
					Coord_CTG_Cov_tmp.contig_no = temp_vec[v];
					Coord_CTG_Cov_tmp.cov = contigs_info->contig_sz_vt[abs(temp_vec[v])];
					Coord_CTG_Cov_tmp.coord2 = v;
					o_sorted_read_contig_idx << Coord_CTG_Cov_tmp.coord << ", ";
					o_sorted_read_contig_idx << Coord_CTG_Cov_tmp.contig_no << ", ";
					o_sorted_read_contig_idx << Coord_CTG_Cov_tmp.cov << ", ";
					o_sorted_read_contig_idx << Coord_CTG_Cov_tmp.coord2 << ", ";
					o_sorted_read_contig_idx << endl;

				}

			}

		}
	}
	

	/**/
	reads_info->LongReadIndexVec.resize(n_selected_read+1);
	reads_info->LenVec.resize(reads_info->LongReadIndexVec.size());
	reads_info->LengthRank.clear();
	reads_info->LengthRank.reserve(reads_info->LongReadIndexVec.size() + 10);

	for (int i = 0; i < reads_info->LongReadIndexVec.size(); ++i)
	{
		int Len = 0;

		int vec_sz = reads_info->LongReadIndexVec[i].size();
		for (int j = 0; j < vec_sz; ++j)
		{
			Len = reads_info->LongReadIndexVec[i][vec_sz - 1].coord - reads_info->LongReadIndexVec[i][0].coord;
		}

		reads_info->LenVec[i] = Len;
	}

	vector<int>::iterator vit = reads_info->LenVec.begin();
	int numReads = reads_info->LenVec.size() - 1;
	vit++;
	for (size_t i = 0; i < numReads; ++i)
	{
		reads_info->LengthRank.push_back(vit);
		vit++;
	}

	sort(reads_info->LengthRank.begin(), reads_info->LengthRank.end(), it_cmp);

//	cout << "";

	reads_info->contained.clear();
	reads_info->contained.resize(numReads + 10);


	reads_info->LeftOverlaps.resize(numReads + 10);
	reads_info->RightOverlaps.resize(numReads + 10);
	reads_info->used_vt_left.resize(numReads + 10);
	reads_info->used_vt_right.resize(numReads + 10);
	reads_info->LeftOverlapsBest.resize(numReads + 10);
	reads_info->RightOverlapsBest.resize(numReads + 10);

	reads_info->LeftOverlapsBestWithTies.resize(numReads + 10);
	reads_info->RightOverlapsBestWithTies.resize(numReads + 10);

	reads_info->both_stand_used.resize(numReads + 10);
	reads_info->chimeric.resize(numReads + 10);

}


void LoadingRawReads(vector<string> filename_vec, vector<string> filename_vec2, reads_info *reads_info, contigs_info *contigs_info)//vector<vector<Coord_CTG_Cov> > &LongReadIndexVec, vector< vector<int> > &Contig2LongRead)
{

	cout << "Loading non-contained sequences." << endl;
	string tag, n_tag, tmp1, tmp2, str, raw_read;
	int coord, contig, cov;
	int n_reads = 0, n_selected = 0,n_updated=0;
	int Len = 0;
	vector<struct Coord_CTG_Cov> TempIndex;
	vector<int> temp_vec, temp_vec_rc;
	TempIndex.reserve(1000);
	temp_vec.reserve(1000);
	temp_vec_rc.reserve(1000);
	Coord_CTG_Cov Coord_CTG_Cov;
	map<int, list<vector<int> > > sorted_index;
	map<vector<int>, int> read_map;
	
	cout<<"Loading sequences for "<<reads_info->Selected_Overlaps.size()<<" non-contained reads."<<endl;
	reads_info->MaxReadLen = 0;
	for (int f = 0; f < filename_vec.size(); ++f)
	{
		string filename = filename_vec[f];
		string filename2 = filename_vec2[f];
		cout << "Loading reads info file: " << filename << endl;
		cout << "Loading reads from: " << filename2 << endl;

		ifstream LongReadIndex_in(filename.c_str());
		ifstream in_NC_Reads(filename2.c_str());
		bool read_success = 1;
		while (read_success)
		{
			read_success = get_a_contig_path(LongReadIndex_in, tag, TempIndex, Len, 0, n_tag);
			n_reads++;
		
			string nc_tag;
			getline(in_NC_Reads, nc_tag);
			getline(in_NC_Reads,raw_read);
			if (raw_read.size() > reads_info->MaxReadLen)
			{
				reads_info->MaxReadLen = raw_read.size();
			}
			temp_vec.clear();
			temp_vec_rc.clear();
			int path_score = 0;
			for (int v = 0; v < TempIndex.size(); ++v)
			{
				temp_vec.push_back(TempIndex[v].contig_no);
				temp_vec_rc.push_back(-TempIndex[TempIndex.size() - 1 - v].contig_no);
				path_score += TempIndex[v].cov;
			}
			bool take_rc = 0;

			for (int v = 0; v < TempIndex.size(); ++v)
			{
				if (temp_vec[v] < temp_vec_rc[v])
				{
					take_rc = 0;
					break;
				}
				if (temp_vec[v] > temp_vec_rc[v])
				{
					take_rc = 1;
					break;
				}
			}
			if (take_rc)
			{
				temp_vec = temp_vec_rc;
				
				reverse_complement_str(raw_read);

			}
			if (reads_info->Selected_Overlaps.count(temp_vec))
			{
				if (reads_info->Selected_Overlaps[temp_vec].size() < 1)//change in the future
				{
					raw_overlap_str raw_overlap_str;
					raw_overlap_str.score = path_score;
					raw_overlap_str.raw_seq = raw_read;
					//cout << raw_read << endl;
					reads_info->Selected_Overlaps[temp_vec].push_back(raw_overlap_str);
					n_selected++;
				}
				if (1)
				if (path_score>reads_info->Selected_Overlaps[temp_vec][0].score)
				{
					reads_info->Selected_Overlaps[temp_vec][0].score = path_score;
					reads_info->Selected_Overlaps[temp_vec][0].raw_seq = raw_read;
					n_updated++;
				}
					
			}
			TempIndex.clear();

			
		}

		cout << n_reads << " reads loaded." << endl;

	}
	//reads_info->LongReadIndexVec.push_back(TempIndex);
	

	cout << n_selected << " reads selected." << endl;
	cout << n_updated << " updated." << endl;
	if (reads_info->Debug)
	{
		n_selected = 0;
		for (int i = 1; i < reads_info->LongReadIndexVec.size(); ++i)
		{
			vector<struct Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[i];
			vector<int> read_vec;
			for (int cc = 0; cc < current_read_layout.size(); ++cc)
			{
				read_vec.push_back(current_read_layout[cc].contig_no);
			}
			if (reads_info->Selected_Overlaps[read_vec].size() > 0)
			{
				n_selected++;
			}
		}
		cout << n_selected << " non-empty." << endl;
	}
	
	
}

void LoadingRefIndex(string filename, reads_info *reads_info, contigs_info *contigs_info)//vector<vector<Coord_CTG_Cov> > &LongReadIndexVec, vector< vector<int> > &Contig2LongRead)
{
	reads_info->RefIndexVec.clear();
	reads_info->Contig2Ref.clear();
	reads_info->Contig2Ref.resize(contigs_info->total_contigs + 1);
	ifstream LongReadIndex_in(filename.c_str());
	string tag, n_tag, tmp1, tmp2, str;
	int coord, contig, cov;
	int n_reads = 0,position_in_read=0;
	vector<Coord_CTG_Cov> TempIndex;
	Coord_CTG_Cov Coord_CTG_Cov;

	while (getline(LongReadIndex_in, str))
	{
		if (str[0] == '>')
		{
			tag = str;
			n_reads++;
			position_in_read = 0;
			reads_info->RefIndexVec.push_back(TempIndex);
			TempIndex.clear();
		}
		else
		{
			stringstream ss(str);
			ss >> coord >> tmp1 >> contig >> tmp2 >> cov;
			if (cov < reads_info->KmerCovTh || cov<(double)contigs_info->contig_sz_vt[abs(contig)]*reads_info->AdaptiveTh)
			{
				continue;
			}
			
			Coord_CTG_Cov.coord = coord;
			Coord_CTG_Cov.contig_no = contig;
			Coord_CTG_Cov.cov = cov;
			if (contig>0)
			{
				reads_info->Contig2Ref[abs(contig)].push_back(n_reads);
			}
			else
			{
				reads_info->Contig2Ref[abs(contig)].push_back(-n_reads);
			}
			reads_info->Contig2Ref[abs(contig)].push_back(position_in_read);
			position_in_read++;
			TempIndex.push_back(Coord_CTG_Cov);
		}
		//cout << n_reads << endl;
	}

	reads_info->RefIndexVec.push_back(TempIndex);
	reads_info->RefLenVec.resize(reads_info->RefIndexVec.size());
	
	for (int i = 0; i < reads_info->RefIndexVec.size(); ++i)
	{
		int Len = 0;
		int vec_sz = reads_info->RefIndexVec[i].size();
		
		for (int j = 0; j < vec_sz; ++j)
		{
			Len = reads_info->RefIndexVec[i][vec_sz - 1].coord - reads_info->RefIndexVec[i][0].coord;
		}

		reads_info->RefLenVec[i] = Len;

	}
}

void OutputOverlapGraph(reads_info *reads_info, string filename)
{
	string graph_name = filename + ".dot";
	ofstream o_graph(graph_name.c_str());
	
	string graph_links = filename + ".txt";
	ofstream o_graph3(graph_links.c_str());
	o_graph3 << reads_info->LeftOverlaps.size() << endl;
	o_graph << "digraph G {" << endl;
	for (int i = 1; i < reads_info->LeftOverlaps.size(); ++i)
	{
		if (reads_info->contained[i]==0)
		for (map<int, int>::iterator tmp_it = reads_info->LeftOverlaps[i].begin(); tmp_it != reads_info->LeftOverlaps[i].end(); ++tmp_it)
		{
			if (reads_info->contained[abs(tmp_it->first)] == 0)
			o_graph << abs(i) << " -> " << abs(tmp_it->first) << ";" << endl;
		}

		for (map<int, int>::iterator tmp_it = reads_info->LeftOverlaps[i].begin(); tmp_it != reads_info->LeftOverlaps[i].end(); ++tmp_it)
		{
			
			o_graph3 << -i << " " << tmp_it->first << " "<< tmp_it->second<< endl;
		}

	}
	for (int i = 1; i < reads_info->RightOverlaps.size(); ++i)
	{
		if (reads_info->contained[i] == 0)
		for (map<int, int>::iterator tmp_it = reads_info->RightOverlaps[i].begin(); tmp_it != reads_info->RightOverlaps[i].end(); ++tmp_it)
		{
			if (reads_info->contained[abs(tmp_it->first)] == 0)
			o_graph << abs(i) << " -> " << abs(tmp_it->first) << ";" << endl;
		}

		for (map<int, int>::iterator tmp_it = reads_info->RightOverlaps[i].begin(); tmp_it != reads_info->RightOverlaps[i].end(); ++tmp_it)
		{

			o_graph3 << i << " " << tmp_it->first << " " << tmp_it->second << endl;
		}

	}
	o_graph << "}" << endl;
}

void OutputBestOverlapGraph(reads_info *reads_info, string filename)
{
	string graph_name = filename + ".dot";
	ofstream o_graph(graph_name.c_str());
	string filename2 = "Strong" + filename + ".dot";
	ofstream o_graph2(filename2.c_str());
	string graph_links = filename + ".txt";
	ofstream o_graph3(graph_links.c_str());

	string graph_name4 = filename + "WithTies.dot";
	ofstream o_graph4(graph_name4.c_str());
	string graph_name5 = filename + "WithTies.txt";
	ofstream o_graph5(graph_name5.c_str());

	o_graph << "digraph G {" << endl;
	o_graph2 << "strict graph G {" << endl;
	o_graph3 << reads_info->LeftOverlapsBest.size() << endl;
	o_graph4 << "digraph G {" << endl;
	for (int i = 1; i < reads_info->LeftOverlapsBest.size(); ++i)
	{
		if (reads_info->LeftOverlapsBest[i] != 0)
		{
			int next_read = reads_info->LeftOverlapsBest[i];
			o_graph << abs(i) << " -> " << abs(next_read) << ";" << endl;

			if (next_read>0 && reads_info->RightOverlapsBest[next_read] == i)
			{
				o_graph2 << abs(i) << " -- " << abs(next_read) << ";" << endl;
			}
			if (next_read<0 && reads_info->LeftOverlapsBest[-next_read] == -i)
			{
				o_graph2 << abs(i) << " -- " << abs(next_read) << ";" << endl;
			}

		}
		if (reads_info->LeftOverlapsBestWithTies[i].size() != 0)
		{
			map<int, vector<int> >::iterator it;

			for (int j = 0; j < reads_info->LeftOverlapsBestWithTies[i].size(); ++j)
			{
				int next_read = reads_info->LeftOverlapsBestWithTies[i][j];
				o_graph4 << abs(i) << " -> " << abs(next_read) << ";" << endl;
				o_graph5 << -i << " " <<next_read << endl;
			}
		}
	}
	for (int i = 1; i < reads_info->RightOverlapsBest.size(); ++i)
	{
		if (reads_info->RightOverlapsBest[i] != 0)
		{
			int next_read = reads_info->RightOverlapsBest[i];
			o_graph << abs(i) << " -> " << abs(next_read) << ";" << endl;

			if (next_read>0 && reads_info->LeftOverlapsBest[next_read] == i)
			{
				o_graph2 << abs(i) << " -- " << abs(next_read) << ";" << endl;
			}
			if (next_read<0 && reads_info->RightOverlapsBest[-next_read] == -i)
			{
				o_graph2 << abs(i) << " -- " << abs(next_read) << ";" << endl;
			}
		}

		if (reads_info->RightOverlapsBestWithTies[i].size() != 0)
		{
			for (int j = 0; j < reads_info->RightOverlapsBestWithTies[i].size(); ++j)
			{
				int next_read = reads_info->RightOverlapsBestWithTies[i][j];
				o_graph4 << abs(i) << " -> " << abs(next_read) << ";" << endl;
				o_graph5 << i << " " << next_read << endl;
			}
		}
	}

	
	for (int i = 1; i < reads_info->RightOverlapsBest.size(); ++i)
	{
		
		if (reads_info->RightOverlapsBest[i] != 0 || reads_info->LeftOverlapsBest[i] != 0)
		{
			if (reads_info->degree_vec.size() == 0||reads_info->degree_vec[i] == 4)
			{
				if (reads_info->chimeric[i])//reads_info->both_stand_used[i]
				{
					o_graph << abs(i) << "[color=blue ];" << endl;
					o_graph2 << abs(i) << "[color=blue ];" << endl;

				}
				else
				{
					o_graph << abs(i)  << endl;
					o_graph2 << abs(i) << endl;
				}
			}
			else
			{

				if( reads_info->chimeric[i])//(reads_info->both_stand_used[i])
				{
					o_graph << abs(i) << "[color=red];" << endl;
					o_graph2 << abs(i) << "[color=red];" << endl;

				}
				else
				{
					o_graph << abs(i) << "[color=red];" << endl;
					o_graph2 << abs(i) << "[color=red];" << endl;

				}
			}
			

		}
	}
	o_graph << "}" << endl;
	o_graph2 << "}" << endl;
	o_graph4 << "}" << endl;


	for (int i = 1; i < reads_info->LeftOverlapsBest.size(); ++i)
	{
		if (reads_info->LeftOverlapsBest[i] != 0 || reads_info->RightOverlapsBest[i] != 0)
		{
			o_graph3 << i << endl;
			o_graph3 << reads_info->LeftOverlapsBest[i] << endl;
			o_graph3 << reads_info->RightOverlapsBest[i] << endl;
		}
	}
}

void LoadBestOverlapGraph(reads_info *reads_info, string filename)
{

	string graph_links = filename + ".txt";
	ifstream in_graph(graph_links.c_str());
	int num_reads = 0;
	in_graph >> num_reads;
	string graph_links2 = filename + "WithTies.txt";
	ifstream in_graph2(graph_links2.c_str());

	reads_info->LeftOverlapsBest.clear();
	reads_info->RightOverlapsBest.clear();
	reads_info->LeftOverlapsBest.resize(num_reads);
	reads_info->RightOverlapsBest.resize(num_reads);
	reads_info->LeftOverlapsBestWithTies.clear();
	reads_info->RightOverlapsBestWithTies.clear();
	reads_info->LeftOverlapsBestWithTies.resize(num_reads);
	reads_info->RightOverlapsBestWithTies.resize(num_reads);
	reads_info->contained.clear();
	reads_info->contained.resize(num_reads);
	reads_info->used.resize(num_reads);
	reads_info->used_vt_left.resize(num_reads);
	reads_info->used_vt_right.resize(num_reads);
	reads_info->both_stand_used.resize(num_reads);
	reads_info->degree_vec.clear();
	reads_info->degree_vec.resize(num_reads + 10);
	int read_idx = 0;
	while (in_graph >> read_idx)
	{
		int next_read = 0;
		in_graph >> next_read;
		reads_info->LeftOverlapsBest[read_idx] = next_read;
		in_graph >> next_read;
		reads_info->RightOverlapsBest[read_idx] = next_read;
	}
	int read1, read2;
	while (in_graph2 >> read1 >> read2)
	{
		if (read1 < 0)
		{
			reads_info->LeftOverlapsBestWithTies[-read1].push_back(read2);
		}
		else
		{
			reads_info->RightOverlapsBestWithTies[read1].push_back(read2);
		}
	}
}



void LoadOverlapGraph(reads_info *reads_info, string filename)
{

	string graph_links = filename + ".txt";
	ifstream in_graph(graph_links.c_str());
	int num_reads = 0;
	in_graph >> num_reads;
	reads_info->LeftOverlaps.clear();
	reads_info->RightOverlaps.clear();
	reads_info->LeftOverlaps.resize(num_reads);
	reads_info->RightOverlaps.resize(num_reads);
	int read_idx = 0;
	int next_read = 0;
	int overlap_score = 0;
	while (in_graph >> read_idx)
	{
		in_graph >> next_read;
		in_graph >> overlap_score;
		if (read_idx < 0)
		{
			reads_info->LeftOverlaps[-read_idx][next_read] = overlap_score;

		}
		if (read_idx>0)
		{
			reads_info->RightOverlaps[read_idx][next_read] = overlap_score;
		}
	}
	
}


void OutputCleanedOverlapGraph(reads_info *reads_info, string filename)
{

	ofstream o_graph(filename.c_str());
	o_graph << "strict graph G {" << endl;
	for (int i = 1; i < reads_info->LeftOverlaps.size(); ++i)
	{
		if (reads_info->LeftOverlaps[i].size())
		{
			for (map<int, int>::iterator tmp_it = reads_info->LeftOverlaps[i].begin(); tmp_it != reads_info->LeftOverlaps[i].end(); ++tmp_it)
			{
				int next_read = tmp_it->first;
				o_graph << abs(i) << " -- " << abs(next_read) << ";" << endl;

			}
		}

	}
	for (int i = 1; i < reads_info->RightOverlaps.size(); ++i)
	{
		if (reads_info->RightOverlaps[i].size())
		{
			for (map<int, int>::iterator tmp_it = reads_info->RightOverlaps[i].begin(); tmp_it != reads_info->RightOverlaps[i].end(); ++tmp_it)
			{
				int next_read = tmp_it->first;
				o_graph << abs(i) << " -- " << abs(next_read) << ";" << endl;

			}
		}
	}
	o_graph << "}" << endl;
}



void ConstructUndirectedBestOverlapGraph(reads_info *reads_info, string filename)
{

	/*
	mode1: add edges to make mutual overlaps.
	mode2: just consider a subgraph with mutual best overlaps(non-mutual best overlaps are skipped).
	results are recorded in Left/RightBestOverlapsTemp file
	*/
	int n_edges_added = 0;
	ofstream o_graph(filename.c_str());
	o_graph << "strict graph G {" << endl;
	reads_info->LeftBestOverlapsTemp.clear();
	reads_info->RightBestOverlapsTemp.clear();
	reads_info->LeftBestOverlapsTemp.resize(reads_info->LeftOverlapsBestWithTies.size());
	reads_info->RightBestOverlapsTemp.resize(reads_info->RightOverlapsBestWithTies.size());

	for (int i = 1; i < reads_info->LeftOverlapsBestWithTies.size(); ++i)
	{


		if (reads_info->contained[abs(i)])
		{
			continue;
		}
		for (int j = 0; j< reads_info->LeftOverlapsBestWithTies[i].size(); ++j)
		{
			int next_read = reads_info->LeftOverlapsBestWithTies[i][j];
			int overlap_score = reads_info->LeftOverlaps[i][next_read];

			//o_graph << abs(next_read) << " -> " << abs(i) << ";" << endl;
			if (reads_info->contained[abs(next_read)])
			{
				continue;
			}
			if (next_read>0)
			{

				if (reads_info->mode == 2)
				{
					bool skip = 1;
					for (int k = 0; k < reads_info->RightOverlapsBestWithTies[next_read].size(); ++k)
					{
						if (reads_info->RightOverlapsBestWithTies[next_read][k] == i)
						{
							skip = 0;
							break;
						}
					}
					if (skip)
					{
						continue;
					}
				}

				reads_info->RightBestOverlapsTemp[next_read][i] = overlap_score;

			}
			else
			{

				if (reads_info->mode == 2)
				{
					bool skip = 1;
					for (int k = 0; k < reads_info->LeftOverlapsBestWithTies[-next_read].size(); ++k)
					{
						if (reads_info->LeftOverlapsBestWithTies[-next_read][k] == -i)
						{
							skip = 0;
							break;
						}
					}
					if (skip)
					{
						continue;
					}
				}


				reads_info->LeftBestOverlapsTemp[-next_read][-i] = overlap_score;

			}
			reads_info->LeftBestOverlapsTemp[i][next_read] = overlap_score;
			o_graph << abs(i) << " -- " << abs(next_read) << ";" << endl;
		}



	}
	for (int i = 1; i < reads_info->RightOverlapsBestWithTies.size(); ++i)
	{

		if (reads_info->contained[abs(i)])
		{
			continue;
		}
		vector<int> cleaned_overlap_vec;
		for (int j = 0; j< reads_info->RightOverlapsBestWithTies[i].size(); ++j)
		{

			int next_read = reads_info->RightOverlapsBestWithTies[i][j];
			int overlap_score = reads_info->RightOverlaps[i][next_read];
			if (reads_info->contained[abs(next_read)])
			{
				continue;
			}

			if (next_read>0)
			{
				if (reads_info->mode == 2)
				{
					bool skip = 1;
					for (int k = 0; k < reads_info->LeftOverlapsBestWithTies[next_read].size(); ++k)
					{
						if (reads_info->LeftOverlapsBestWithTies[next_read][k] == i)
						{
							skip = 0;
							break;
						}
					}
					if (skip)
					{
						continue;
					}
				}
				reads_info->LeftBestOverlapsTemp[next_read][i] = overlap_score;

			}
			else
			{
				if (reads_info->mode == 2)
				{
					bool skip = 1;
					for (int k = 0; k < reads_info->RightOverlapsBestWithTies[-next_read].size(); ++k)
					{
						if (reads_info->RightOverlapsBestWithTies[-next_read][k] == -i)
						{
							skip = 0;
							break;
						}
					}
					if (skip)
					{
						continue;
					}
				}


				reads_info->RightBestOverlapsTemp[-next_read][-i] = overlap_score;

			}
			reads_info->RightBestOverlapsTemp[i][next_read] = overlap_score;
			o_graph << abs(i) << " -- " << abs(next_read) << ";" << endl;
		}

	}
	o_graph << "}" << endl;
	int n_tips = 0;
	for (int i = 1; i < reads_info->LeftBestOverlapsTemp.size(); ++i)
	{
		if ((reads_info->LeftBestOverlapsTemp[i].size()>0 || reads_info->RightBestOverlapsTemp[i].size() > 0) && (reads_info->LeftBestOverlapsTemp[i].size() == 0 || reads_info->RightBestOverlapsTemp[i].size() == 0))
		{
			n_tips++;
		}
	}

	cout << n_tips << " tips in the graph." << endl;
}


/*this function analyzes the Left/RightOverlaps and calculate the best overlaps and stores into Left/RightBestOverlaps(WithTies)*/
void ConstructCleanedBestOverlapGraph(reads_info *reads_info, contigs_info *contigs_info)
{


	for (int i = 1; i < reads_info->LeftOverlaps.size(); ++i)
	{
		if (reads_info->LeftOverlaps[i].size() == 0 && reads_info->RightOverlaps[i].size() == 0)
		{
			reads_info->contained[i] = 1;//these nodes will not be used
		}
	}
	int n_chimera = 0;
	if (reads_info->RemoveChimera)
	for (int i = 1; i < reads_info->chimeric.size(); ++i)
	{
		if (reads_info->chimeric[i] == 1 && reads_info->mode == 1)
		{
			
			reads_info->contained[i] = 1;//these nodes will not be used
			
			n_chimera++;

		}
	}
	cout << n_chimera << " chimeric reads deleted." << endl;

	int numReads = reads_info->LenVec.size() - 1;
	reads_info->degree_vec.clear();
	reads_info->degree_vec.resize(numReads + 10);
	reads_info->LeftOverlapsBest.clear();
	reads_info->LeftOverlapsBest.resize(numReads + 10);
	reads_info->RightOverlapsBest.clear();
	reads_info->RightOverlapsBest.resize(numReads + 10);
	reads_info->LeftOverlapsBestWithTies.clear();
	reads_info->LeftOverlapsBestWithTies.resize(numReads + 10);
	reads_info->RightOverlapsBestWithTies.clear();
	reads_info->RightOverlapsBestWithTies.resize(numReads + 10);



	for (int i = 0; i < numReads; ++i)
	{
		int current_read = (reads_info->LengthRank[i]) - reads_info->LenVec.begin();


		if (reads_info->contained[current_read])// && reads_info->linked[current_read] == 0)
		{
			continue;
		}
		map<int, int> local_index;

		for (int c = 0; c < reads_info->LongReadIndexVec[abs(current_read)].size(); ++c)
		{
			int ctg_no = reads_info->LongReadIndexVec[abs(current_read)][c].contig_no;
			int coord = reads_info->LongReadIndexVec[abs(current_read)][c].coord;
			if (local_index.count(ctg_no) == 0)//only one copy is saved.
			{
				local_index[ctg_no] = (coord);
			}

		}
		vector<Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[abs(current_read)];

		map<int, int> candidate_reads;
		map<int, vector<int> > overlap_index;
		map<int, vector<Coord_CTG_Cov> > candidate_reads_layouts;

		map<int, vector<int> > sorted_left_overlaps, sorted_right_overlaps;

		/////edit
		map<int, int>::iterator tmp_it1;
		for (tmp_it1 = reads_info->LeftOverlaps[abs(current_read)].begin(); tmp_it1 != reads_info->LeftOverlaps[abs(current_read)].end(); ++tmp_it1)
		{
			if (!reads_info->contained[abs(tmp_it1->first)])
			{
				sorted_left_overlaps[tmp_it1->second].push_back(tmp_it1->first);
			}
		}
		for (tmp_it1 = reads_info->RightOverlaps[abs(current_read)].begin(); tmp_it1 != reads_info->RightOverlaps[abs(current_read)].end(); ++tmp_it1)
		{
			if (!reads_info->contained[abs(tmp_it1->first)])
			{
				sorted_right_overlaps[tmp_it1->second].push_back(tmp_it1->first);
			}
		}

		


		if (sorted_left_overlaps.size() > 0)
		{
			map<int, vector<int> >::reverse_iterator rit = sorted_left_overlaps.rbegin();
			vector<int > best_overlap_candidate = rit->second;
			int next_read = best_overlap_candidate[0];
			for (int c = 0; c < best_overlap_candidate.size(); ++c)
			{
				reads_info->LeftOverlapsBestWithTies[current_read].push_back(best_overlap_candidate[c]);
			}

			reads_info->LeftOverlapsBest[current_read] = next_read;
			reads_info->degree_vec[abs(current_read)]++;
			reads_info->degree_vec[abs(next_read)]++;

		}
		if (sorted_right_overlaps.size() > 0)
		{

			map<int, vector<int> >::reverse_iterator rit = sorted_right_overlaps.rbegin();
			vector<int > best_overlap_candidate = rit->second;
			int next_read = best_overlap_candidate[0];
			for (int c = 0; c < best_overlap_candidate.size(); ++c)
			{
				reads_info->RightOverlapsBestWithTies[current_read].push_back(best_overlap_candidate[c]);
			}

			reads_info->RightOverlapsBest[current_read] = next_read;
			reads_info->degree_vec[abs(current_read)]++;
			reads_info->degree_vec[abs(next_read)]++;



		}

	}



}


void MultipleAlignmentErrorCorrection(reads_info *reads_info, contigs_info *contigs_info, align_info *align_info_default)
{

	/*
	This algorithm isn't acclerated because we want to clean each read (also the contained ones).
	
	This function is critical for pacbio reads.
	1. Remove chimeric reads and reads with adaptors
	2. Remove false positive contigs in each converted read
	*/

	bool Debug = 0;
	Debug = reads_info->Debug;
	align_matrices align_matrices;

	//Debug = 0;
	align_info align_info;
	align_info = *align_info_default;
	align_info.flip = 0;
	align_info.force_flip = 0;
	//align_info.fix_orientation = 1;
	align_info.MaxOffset = 1* align_info.MaxOffset;

	align_info.min_overlap = align_info.min_overlap / 1;
	align_info.mm_ratio = align_info_default->mm_ratio;
	ofstream MultipleAlignment_debug, MultipleAlignment_debug2;
	ofstream o_cleaned_reads,o_chimera;
	size_t read_cnt = 0;
	//ofstream o_overlap_relations("AllOverlaps.txt");
	if (Debug)
	{
		MultipleAlignment_debug.open("MultipleAlignment_debug.txt");
		MultipleAlignment_debug2.open("MultipleAlignment_debug2.txt");
	}

	o_cleaned_reads.open("CleanedLongReads.txt");
	o_chimera.open("ChimeraIdx.txt");
	vector< vector<int> > Contig2LongRead_corr;

	Contig2LongRead_corr.resize(reads_info->Contig2LongRead.size());
	int numReads = reads_info->LenVec.size() - 1;
	reads_info->contig2contig_cov.clear();

	int64_t n_alignments = 0;
	int64_t seq_len_sum = 0, seq_len_sum_compressed=0;


	for (int i = 1; i < numReads; ++i)
	{
		int current_read = i;// (reads_info->LengthRank[i]) - reads_info->LenVec.begin();

	
		if (i % 1000000 == 0&&i>0)
		{
			cout << i << " sequences aligned." << endl;
			cout << "Avg alignment size: " << seq_len_sum / 2 / n_alignments << endl;
			cout << "total alignments: " << n_alignments << endl;
		}
	
		if (reads_info->LongReadIndexVec[current_read].size() == 1)
		{
		//	reads_info->contained[current_read] = 1;
			int contig_no = abs(reads_info->LongReadIndexVec[current_read][0].contig_no);
			Contig2LongRead_corr[contig_no].push_back(current_read);
			vector<Coord_CTG_Cov> read_layout_corrected;
			read_layout_corrected = reads_info->LongReadIndexVec[current_read];
			if (read_layout_corrected.size() > 0)
			{
				read_cnt++;
				o_cleaned_reads  << reads_info->tag_vec[current_read] << endl;
				o_cleaned_reads << reads_info->LenVec[current_read] << endl;

				for (int r = 0; r < read_layout_corrected.size(); ++r)
				{
					o_cleaned_reads << read_layout_corrected[r].coord << ", ";
					o_cleaned_reads << read_layout_corrected[r].contig_no << ", ";
					o_cleaned_reads << read_layout_corrected[r].cov << ", ";
					o_cleaned_reads << read_layout_corrected[r].coord2 << endl;

				}
			}
			
			continue;
		}

		vector<Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[abs(current_read)];

		map<int, int> candidate_reads;
		map<int, vector<int> > overlap_index;
		map<int, vector<Coord_CTG_Cov> > candidate_reads_layouts;
		map<int, vector<int> > sorted_candidate_reads;

		for (int c1 = 0; c1 < reads_info->LongReadIndexVec[abs(current_read)].size(); ++c1)
		{

			int ctg_no = reads_info->LongReadIndexVec[abs(current_read)][c1].contig_no;
			if (reads_info->Contig2LongRead[abs(ctg_no)].size()<reads_info->max_reads)
			for (int r = 0; r < reads_info->Contig2LongRead[abs(ctg_no)].size(); ++r)
			{

				int read_idx = reads_info->Contig2LongRead[abs(ctg_no)][r];
				if (read_idx == current_read)
				{
					continue;
				}
				if (reads_info->LongReadIndexVec[read_idx].size()>1)
				if (abs(reads_info->LongReadIndexVec[read_idx][0].coord2 - reads_info->LongReadIndexVec[read_idx][reads_info->LongReadIndexVec[read_idx].size() - 1].coord2) < reads_info->MinLen)
				{
					//continue;
				}

				//candidate_reads[read_idx]++;
				if (reads_info->LongReadIndexVec[abs(read_idx)].size()>1)//non-contained
				{
					if (align_info.matching_method == 1)
					{
						candidate_reads[read_idx] += contigs_info->contig_sz_vt[abs(ctg_no)];

					}
					if (align_info.matching_method == 2 || align_info.matching_method == 3)
					{
						candidate_reads[read_idx] += reads_info->LongReadIndexVec[abs(current_read)][c1].cov;

					}
				}
				

			}

		}
		//cout << "";

		for (map<int, int>::iterator it1 = candidate_reads.begin(); it1 != candidate_reads.end(); ++it1)
		{
			sorted_candidate_reads[it1->second].push_back(it1->first);
		}

		int n_reads = 0;
		for (map<int, vector<int> >::reverse_iterator it2 = sorted_candidate_reads.rbegin(); it2 != sorted_candidate_reads.rend(); ++it2)
		{
			if (it2->first < align_info.min_overlap)
			{
				break;//stop if the overlap is too small.
			}
			//vector<int> tmp_vec = it2->second;
			for (int jj = 0; jj < it2->second.size(); ++jj)
			{
				if (n_reads>2000)//reads_info->TopOverlaps)
				{
					break;//don't break here
				}
				candidate_reads_layouts[it2->second[jj]] = reads_info->LongReadIndexVec[it2->second[jj]];// candidate_reads_layouts2[it2->second[jj]];
				n_reads++;

			}
		}

		//candidates sorted
		map<int, vector< vector<int> > >  aligned_qry_vec;
		map<int, vector<string> >  aligned_qry_vec_tag, aligned_qry_vec_tag2;
		
		
		map<int, int> matching_cov;
		vector<vector<int> > match_vec;
		match_vec.reserve(2000);

		for (map<int, vector<Coord_CTG_Cov> >::iterator it1 = candidate_reads_layouts.begin(); it1 != candidate_reads_layouts.end(); ++it1)
		{
			int read_idx = it1->first;
			if (read_idx == 5612)
			{
				cout << "";
			}
			vector<Coord_CTG_Cov> Coord_CTG_Cov_Vec = it1->second;
			int result = 0;
			
			int orientation = 0;
			int match_ctgs1 = 0, match_ctgs2 = 0;
			if (align_info.fix_orientation)
			{
				map<int, int> local_index, ref_index, qry_index;
				int ref_sz = current_read_layout.size(), qry_sz = Coord_CTG_Cov_Vec.size();

				for (int i = 0; i < ref_sz; ++i)
				{
					local_index[current_read_layout[i].contig_no] = current_read_layout[i].coord;
				}

				for (int i = 0; i < qry_sz; ++i)
				{
					if (local_index.count(Coord_CTG_Cov_Vec[i].contig_no))
					{
						match_ctgs1++;
					}
					if (local_index.count(-Coord_CTG_Cov_Vec[i].contig_no))
					{
						match_ctgs2++;
					}
				}
			}
			


			
			for (int orientation = 0; orientation <= 1; ++orientation)
			{
				if (align_info.fix_orientation)
				{
					if (orientation == 0)
					{
						align_info.force_flip = 0;

						if (match_ctgs1 == 0)
						{
							continue;
						}

					}
					else
					{
						align_info.force_flip = 1;

						if (match_ctgs2 == 0)
						{
							continue;
						}
					}

				}
				else
				{
					if (orientation > 0)
					{
						break;
					}
				}

	
				align_info.ref_len = reads_info->LenVec[abs(current_read)];
				align_info.qry_len = reads_info->LenVec[abs(read_idx)];

				n_alignments++;
				seq_len_sum += current_read_layout.size();
				seq_len_sum += Coord_CTG_Cov_Vec.size();

				//
				if (align_info.SparseAlign)
				{
					sparse_semi_global_align(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_matrices, &align_info);
				}
				else
				{
					approximate_semi_global_align(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_matrices, &align_info);
				}

				
				seq_len_sum_compressed += align_info.ref_len_compressed;
				seq_len_sum_compressed += align_info.qry_len_compressed;
				//align_info.band_width = max(current_read_layout.size(),Coord_CTG_Cov_Vec.size());
				//approximate_semi_global_align_banded(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_maps, &align_info);

				result = classify_alignment(&align_info);

				// approximate_semi_local_align_banded(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_maps, &align_info);

				//align2ref(contigs_info, reads_info, current_read_layout, &align_maps, &align_info);

				if (result >= 1 && result <= 4)
				{
					int matching_ctgs = 0;
					for (int q = 0; q < align_info.qry_aligned.size(); ++q)
					{
						if (align_info.ref_aligned[q] != 0 && align_info.qry_aligned[q] != 0)
						{
							matching_ctgs++;
						}

					}
					if (matching_ctgs <= 1)
					{
						continue;
					}


					int overlaped_read = abs(read_idx);
					if (align_info.flip)
					{
						overlaped_read = -overlaped_read;
					}
					//o_overlap_relations << current_read << " " << overlaped_read << " "<<align_info.offset << endl;



					int matching_pos = 0, min_match = -1, max_match = -1;
					vector<int> tmp_vec;
					for (int q = 0; q < align_info.qry_aligned.size(); ++q)
					{
						if (align_info.ref_aligned[q] != 0)
						{

							if (align_info.qry_aligned[q] != 0)
							{
								if (min_match < 0)
								{
									min_match = q;
								}
								if (max_match < q)
								{
									max_match = q;
								}
								matching_cov[matching_pos]++;
							}

							matching_pos++;

							//MultipleAlignment_debug << align_info.qry_aligned[q] << ", ";
							tmp_vec.push_back(align_info.qry_aligned[q]);
						}


					}
					match_vec.push_back(tmp_vec);

					if (Debug)
					{
						stringstream itoa_str;
						string temp_str1, temp_str2;
						itoa_str << read_idx;
						temp_str1 = itoa_str.str();
						stringstream itoa_str2;
						itoa_str2 << result;
						temp_str2 = itoa_str2.str();
						string qry_tag = ">Qry_read_" + temp_str1 + "_type:" + temp_str2;
						aligned_qry_vec_tag[min_match].push_back(qry_tag);
						//aligned_qry_vec_tag2[min_match].push_back(reads_info->tag_vec[abs(read_idx)]);
						aligned_qry_vec[min_match].push_back(tmp_vec);

						MultipleAlignment_debug2 << ">Ref_read_" << current_read << endl;

						MultipleAlignment_debug2 << reads_info->tag_vec[abs(current_read)] << endl;
						for (int q = 0; q < align_info.ref_aligned.size(); ++q)
						{
							MultipleAlignment_debug2 << align_info.ref_aligned[q] << ", ";
						}
						MultipleAlignment_debug2 << endl;

						MultipleAlignment_debug2 << ">Qry_read_" << read_idx << "_type:" << result << endl;
						MultipleAlignment_debug2 << reads_info->tag_vec[abs(read_idx)] << endl;
						for (int q = 0; q < align_info.qry_aligned.size(); ++q)
						{
							MultipleAlignment_debug2 << align_info.qry_aligned[q] << ", ";
						}
						MultipleAlignment_debug2 << endl;


					}


				}

			}



		}


		vector<Coord_CTG_Cov> read_layout_corrected;
		vector<vector<int > > match_vec_corr;
		vector<int> cov_corr;
		int match_reads = match_vec.size();
		match_vec_corr.resize(match_reads);
		for (int r = 0; r < current_read_layout.size(); ++r)
		{
			if (matching_cov[r] >= contigs_info->ContigCovTh)
			{
				read_layout_corrected.push_back(current_read_layout[r]);
				cov_corr.push_back(matching_cov[r]);
				for (int rr = 0; rr < match_reads; ++rr)
				{
					match_vec_corr[rr].push_back(match_vec[rr][r]);
				}
			}
			
		}


		if (Debug)
		{
			MultipleAlignment_debug << ">Ref_read_" << current_read << endl;
			//MultipleAlignment_debug << reads_info->tag_vec[abs(current_read)] << endl;

			for (int r = 0; r < current_read_layout.size(); ++r)
			{
				MultipleAlignment_debug << current_read_layout[r].contig_no << ", ";
			}
			MultipleAlignment_debug << endl;

			map<int, vector<vector<int> > >::iterator tmp_it;
			map<int, vector<string> > ::iterator tmp_it2;
			tmp_it2 = aligned_qry_vec_tag.begin();
			for (tmp_it = aligned_qry_vec.begin(); tmp_it != aligned_qry_vec.end(); ++tmp_it)
			{
				for (int r = 0; r < tmp_it->second.size(); ++r)
				{
					MultipleAlignment_debug << tmp_it2->second[r] << endl;

					for (int c = 0; c < tmp_it->second[r].size(); ++c)
					{
						MultipleAlignment_debug << tmp_it->second[r][c] << ", ";
					}
					MultipleAlignment_debug << endl;
				}
				tmp_it2++;

			}

			for (int r = 0; r < current_read_layout.size(); ++r)
			{
				MultipleAlignment_debug << matching_cov[r] << ", ";
			}

			MultipleAlignment_debug << endl;
		
		}

		
		
		vector<int> cov_vec;
		if (read_layout_corrected.size()>1)
		{
			cov_vec.resize(read_layout_corrected.size()-1);
		}
		for (int r = 0; r < match_vec_corr.size(); ++r)
		{
			int min_hit = -1, max_hit = -2;
			for (int rr = 0; rr < match_vec_corr[r].size(); ++rr)
			{
				if (match_vec_corr[r][rr]!=0)
				{
					if (min_hit < 0)
					{
						min_hit = rr;
					}
					max_hit = rr;
				}
			}

			for (int rr = min_hit; rr < max_hit; ++rr)
			{
				cov_vec[rr]++;
			}
		}
	
		for (int r = 0; r < cov_vec.size(); ++r)
		{
			if (cov_vec[r] < reads_info->ChimeraTh)
			{
				reads_info->chimeric[abs(current_read)] = 1;

				//a trick here that changes performance a bit
				if (reads_info->RemoveChimera)
				{					
					reads_info->contained[abs(current_read)] = 1;
					
				}
			}
			//cout << cov_vec[r] << ", ";
		}
		Coord_CTG_Cov Coord_CTG_Cov_best;
		Coord_CTG_Cov_best.cov = 0;
		if (read_layout_corrected.size() <=1 && reads_info->chimeric[abs(current_read)] == 0)
		{
			

			if(1)
			{
				//could be used for gap closing
				reads_info->chimeric[abs(current_read)] = 1;
				reads_info->contained[abs(current_read)] = 1;
				read_layout_corrected = current_read_layout;
				cov_vec.clear();
				cov_corr.clear();
				cov_vec.resize(current_read_layout.size());
				cov_corr.resize(current_read_layout.size());
			}
			else
			{
				//pick the best match for the contained reads

				for (int l = 0; l < current_read_layout.size(); ++l)
				{
					if (current_read_layout[l].cov>Coord_CTG_Cov_best.cov)
					{
						int AdaptiveTh = 3 * reads_info->AdaptiveTh*contigs_info->contig_sz_vt[abs(current_read_layout[l].contig_no)];
						int KmerCovTh = reads_info->KmerCovTh * 3;
						if ((current_read_layout[l].cov > KmerCovTh) && (current_read_layout[l].cov > AdaptiveTh))
						{
							Coord_CTG_Cov_best = current_read_layout[l];
						}

					}
				}
			}
			
		}
		if (Coord_CTG_Cov_best.cov > 0)
		{
			read_layout_corrected.push_back(Coord_CTG_Cov_best);
			cov_vec.resize(1);
			cov_vec[0] = Coord_CTG_Cov_best.cov;
			cov_corr.resize(1);
			cov_corr[0] = Coord_CTG_Cov_best.cov;
		}
		
		//it is safer but worse if we don't use the following part
		if (0)//reads_info->RemoveChimera&&reads_info->chimeric[abs(current_read)])
		{
			int begin_cnt = 0, rear_cnt = 0, begin_len = 0, rear_len = 0;
			for (int r = 0; r < cov_vec.size(); ++r)
			{
				if (cov_vec[r] == 0)
				{
					begin_cnt = r;
					begin_len = abs(read_layout_corrected[begin_cnt].coord - read_layout_corrected[0].coord) + contigs_info->contig_sz_vt[abs(read_layout_corrected[0].contig_no)] / 2 + contigs_info->contig_sz_vt[abs(read_layout_corrected[begin_cnt].contig_no)] / 2;
					break;
				}
			}
			for (int r = (int)cov_vec.size()-1; r >= 0; --r)
			{
				if (cov_vec[r] == 0)
				{
					rear_cnt =  r+1;
					rear_len = abs(read_layout_corrected[cov_vec.size()].coord - read_layout_corrected[rear_cnt].coord) + contigs_info->contig_sz_vt[abs(read_layout_corrected[cov_vec.size()].contig_no)] / 2 + contigs_info->contig_sz_vt[abs(read_layout_corrected[rear_cnt].contig_no)] / 2;
					break;
				}
			}
			if (begin_len >= rear_len)
			{
				vector<Coord_CTG_Cov> tmp_vec;
				for (int r = 0; r <= begin_cnt; ++r)
				{
					tmp_vec.push_back(read_layout_corrected[r]);
				}
				read_layout_corrected = tmp_vec;
			}
			else
			{
				vector<Coord_CTG_Cov> tmp_vec;
				for (int r = rear_cnt; r <= cov_vec.size(); ++r)
				{
					tmp_vec.push_back(read_layout_corrected[r]);
				}
				read_layout_corrected = tmp_vec;
			}
			reads_info->chimeric[abs(current_read)] = 0;
			reads_info->contained[abs(current_read)] = 0;
		}
		//cout << endl;
		
		reads_info->contig2contig_cov[abs(current_read)] = cov_vec;
		
		
		//latest update!!!
		
		reads_info->LongReadIndexVec[abs(current_read)] = read_layout_corrected;

		
		//check above
		
		if (1)
		{
			if (read_layout_corrected.size() > 0)
			{
				read_cnt++;
				o_cleaned_reads<< reads_info->tag_vec[current_read] << endl;
				o_cleaned_reads << reads_info->LenVec[current_read] << endl;

				for (int r = 0; r < read_layout_corrected.size(); ++r)
				{
					o_cleaned_reads << read_layout_corrected[r].coord << ", ";
					o_cleaned_reads << read_layout_corrected[r].contig_no << ", ";
					o_cleaned_reads << read_layout_corrected[r].cov << ", ";
					o_cleaned_reads << read_layout_corrected[r].coord2 << endl;

				}
				if (reads_info->chimeric[current_read])
				{
					o_chimera << read_cnt << endl;
				}
			}

		}



	}
	cout << "Avg alignment size: " << seq_len_sum / 2 / n_alignments << endl;
	if (align_info.SparseAlign)
	cout << "Avg sparse alignment size: " << seq_len_sum_compressed / 2 / n_alignments << endl;
	cout << "total alignments: " << n_alignments << endl;


	reads_info->Contig2LongRead = Contig2LongRead_corr;

}



void RecoverFalseNegatives(reads_info *reads_info, contigs_info *contigs_info, align_info *align_info_default)
{
	cout << "Recover false negative contigs." << endl;

	ofstream debug_out;
	ofstream debug_out2;


	align_info align_info;
	align_info.flip = 0;
	align_info.fix_orientation = 0;
	align_info.force_flip = 0;
	align_info = *align_info_default;
	align_matrices align_matrices;
	bool Debug = reads_info->Debug;
	Debug = 0;
	if (Debug)
	{
		debug_out.open("RecoverFalseNegatives_debug.txt"); debug_out2.open("RecoverFalseNegativesCorr_debug.txt");
	}

//	align_info.mm_penalty = 0.3;
//	align_info.mm_ratio = 1.0;
//	align_info.max_mismatch_FR = align_info_default->min_overlap * 3;
//	align_info.min_overlap = align_info_default->min_overlap / 2;
//	align_info.min_extension_FR = align_info_default->min_overlap * 0;
	/*
	for each read if it is longer than a threshold, then check all the reads that potentially overlap it,
	if the overlap is larger than a threshold and it looks like a prefix-suffix overlap, then build a link.
	Construct the graph from the longest informative reads.
	*/

	int numReads = reads_info->LenVec.size() - 1;



	int64_t n_alignments = 0;
	int64_t seq_len_sum = 0, seq_len_sum_compressed=0;

	cout << "Calculating reads overlaps." << endl;
	time_t beg_time, end_time;
	time(&beg_time);

	align_info.CalculateOffset = 0;

	for (int i = 0; i < numReads; ++i)
	{

		if ((i % 1000000 == 0) && (n_alignments>0))
		{
			cout << i << " reads aligned." << endl;
			cout << "Avg alignment size: " << seq_len_sum / 2 / n_alignments << endl;
			cout << "total alignments: " << n_alignments << endl;
		}
		int current_read = i;// (reads_info->LengthRank[i]) - reads_info->LenVec.begin();

		if (reads_info->contained[current_read])//contained reads are skipped in any round, which may lead to errors
		{
			continue;
		}

		vector<Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[abs(current_read)];

		map<int, int> candidate_reads;//, left_ext_map, right_ext_map;
		map<int, vector<Coord_CTG_Cov> > candidate_reads_layouts;
		map<int, vector<int> > sorted_candidate_reads;
		vector<vector<int> > collected_alignments_ref,collected_alignments_qry;
		vector<map<int, int> > contig_score_vec;
		contig_score_vec.resize(reads_info->LongReadIndexVec[abs(current_read)].size());
		for (int c1 = 0; c1 < reads_info->LongReadIndexVec[abs(current_read)].size(); ++c1)
		{
			int ctg_no = reads_info->LongReadIndexVec[abs(current_read)][c1].contig_no;

			if (reads_info->Contig2LongRead[abs(ctg_no)].size() < reads_info->max_reads)
			for (int r = 0; r < reads_info->Contig2LongRead[abs(ctg_no)].size(); ++r)
			{

				int read_idx = reads_info->Contig2LongRead[abs(ctg_no)][r];

				if (read_idx == current_read)
				{
					continue;
				}
				if (align_info.matching_method == 1)
				{
					candidate_reads[read_idx] += contigs_info->contig_sz_vt[abs(ctg_no)];
				}
				if (align_info.matching_method == 2 || align_info.matching_method == 3)
				{
					candidate_reads[read_idx] += reads_info->LongReadIndexVec[abs(current_read)][c1].cov;
				}
			}

		}

		for (map<int, int>::iterator it1 = candidate_reads.begin(); it1 != candidate_reads.end(); ++it1)
		{
			if (it1->second < align_info.min_overlap)
			{
				continue;//stop if the overlap is too small.
			}
			sorted_candidate_reads[it1->second].push_back(it1->first);
		}

		int n_reads = 0;
		for (map<int, vector<int> >::reverse_iterator it2 = sorted_candidate_reads.rbegin(); it2 != sorted_candidate_reads.rend(); ++it2)
		{

			//vector<int> tmp_vec = it2->second;
			for (int jj = 0; jj < it2->second.size(); ++jj)
			{
				if (n_reads>10000)
				{
					break;
				}
				candidate_reads_layouts[it2->second[jj]] = reads_info->LongReadIndexVec[it2->second[jj]];// candidate_reads_layouts2[it2->second[jj]];
				n_reads++;

			}
		}

		//map<int, int> sorted_left_overlaps, sorted_right_overlaps;
		map<int, vector<int> > sorted_left_overlaps, sorted_right_overlaps;

		for (map<int, vector<Coord_CTG_Cov> >::iterator it1 = candidate_reads_layouts.begin(); it1 != candidate_reads_layouts.end(); ++it1)
		{
			int read_idx = it1->first;
			vector<Coord_CTG_Cov> Coord_CTG_Cov_Vec = it1->second;

			n_alignments++;
			seq_len_sum += current_read_layout.size();
			seq_len_sum += Coord_CTG_Cov_Vec.size();
	

			if (align_info.SparseAlign)
			{
				sparse_semi_global_align(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_matrices, &align_info);
				seq_len_sum_compressed += align_info.ref_len_compressed;
				seq_len_sum_compressed += align_info.qry_len_compressed;
			}
			else
			{
				approximate_semi_global_align(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_matrices, &align_info);

			}
			



			int result = classify_alignment(&align_info);


			if (result >= 1 && result <= 4)
			{
				


				int matching_ctgs = 0;
				for (int q = 0; q < align_info.qry_aligned.size(); ++q)
				{
					if (align_info.ref_aligned[q] != 0 && align_info.qry_aligned[q] != 0)
					{
						matching_ctgs++;
					}

				}
				if (matching_ctgs <= 1)
				{
					continue;
				}
				
				collected_alignments_ref.push_back(align_info.ref_aligned);
				collected_alignments_qry.push_back(align_info.qry_aligned);


			}

		}


		for (int n = 0; n < collected_alignments_qry.size(); ++n)
		{

			int matching_pos = 0, min_match = -1, max_match = -1;
			vector<int> tmp_vec;
			for (int q = 0; q < collected_alignments_qry[n].size(); ++q)
			{
				if (collected_alignments_ref[n][q] != 0)
				{
					if (matching_pos>0)
					{
						for (int ii = 0; ii < tmp_vec.size(); ++ii)
						{
							contig_score_vec[matching_pos - 1][tmp_vec[ii]]++;//score at each location
						}

						tmp_vec.clear();
					}

					matching_pos++;

				}
				else
				{
					if (matching_pos > 0)
					{
						if (collected_alignments_qry[n][q] != 0)
						{
							tmp_vec.push_back(collected_alignments_qry[n][q]);
						}

					}
				}


			}
		}

		vector<vector<int> > false_nagative_path;
		vector<int> max_score_vec;
		max_score_vec.resize(reads_info->LongReadIndexVec[abs(current_read)].size());
		false_nagative_path.resize(reads_info->LongReadIndexVec[abs(current_read)].size()); 

		for (int n = 0; n < collected_alignments_qry.size(); ++n)
		{

			int matching_pos = 0, min_match = -1, max_match = -1;
			vector<int> tmp_vec;
			int tmp_score = 0;
			for (int q = 0; q < collected_alignments_qry[n].size(); ++q)
			{
				if (collected_alignments_ref[n][q] != 0)
				{
					if (matching_pos>0)
					{
						for (int ii = 0; ii < tmp_vec.size(); ++ii)
						{
							tmp_score+=contig_score_vec[matching_pos - 1][tmp_vec[ii]];//score at each location
						}
						if (tmp_score>max_score_vec[matching_pos - 1])
						{
							max_score_vec[matching_pos - 1] = tmp_score;
							false_nagative_path[matching_pos - 1] = tmp_vec;
						}

						tmp_vec.clear();
						tmp_score = 0;
					}

					matching_pos++;
					
				}
				else
				{
					if (matching_pos > 0)
					{
						if (collected_alignments_qry[n][q] != 0)
						{
							tmp_vec.push_back(collected_alignments_qry[n][q]);
						}

					}
				}


			}
		}

		if (Debug)
		{
			debug_out <<"read: "<< i << endl;
			for (int n = 0; n <collected_alignments_qry.size(); ++n)
			{
				for (int q = 0; q < collected_alignments_ref[n].size(); ++q)
				{
					debug_out << collected_alignments_ref[n][q] << ", ";

				}
				debug_out << endl;
				for (int q = 0; q < collected_alignments_qry[n].size(); ++q)
				{
					debug_out << collected_alignments_qry[n][q] << ", ";

				}
				debug_out << endl;
			}

			for (int n = 0; n+1 < max_score_vec.size(); ++n)
			{
				debug_out << "Position: " << n << ", score: "<<max_score_vec[n]<<endl;
				for (int q = 0; q < false_nagative_path[n].size(); ++q)
				{
					debug_out << false_nagative_path[n][q] << ", ";
				}
				debug_out << endl;
			}

		}

		vector<Coord_CTG_Cov> corr_vec;
		if (reads_info->LongReadIndexVec[abs(current_read)].size()>0)
		{
			corr_vec.push_back(reads_info->LongReadIndexVec[abs(current_read)][0]);
		}
		for (int j = 1; j < reads_info->LongReadIndexVec[abs(current_read)].size(); ++j)
		{



			int Ns = abs(reads_info->LongReadIndexVec[i][j].coord - reads_info->LongReadIndexVec[i][j - 1].coord);
			int dist = (contigs_info->contig_sz_vt[abs(reads_info->LongReadIndexVec[i][j].contig_no)] + contigs_info->contig_sz_vt[abs(reads_info->LongReadIndexVec[i][j - 1].contig_no)]) / 2;
			Ns -= dist;

			if (Ns > 1000)
			{
				//try fill
				if (max_score_vec[j - 1] > 0)
				{
					int mass = 0;
					for (int q = 0; q < false_nagative_path[j - 1].size(); ++q)
					{
						mass += contigs_info->contig_sz_vt[abs(false_nagative_path[j - 1][q])];
					}
					int center = reads_info->LongReadIndexVec[i][j - 1].coord + (contigs_info->contig_sz_vt[abs(reads_info->LongReadIndexVec[i][j - 1].contig_no)] + Ns) / 2;
					int start = center - mass / 2;
					for (int q = 0; q < false_nagative_path[j - 1].size(); ++q)
					{
						struct Coord_CTG_Cov Coord_CTG_Cov_tmp;
						Coord_CTG_Cov_tmp.coord = start + contigs_info->contig_sz_vt[abs(reads_info->LongReadIndexVec[i][j - 1].contig_no)]/ 2;
						Coord_CTG_Cov_tmp.contig_no = false_nagative_path[j - 1][q];
						Coord_CTG_Cov_tmp.cov = 1;
						Coord_CTG_Cov_tmp.coord2 = Coord_CTG_Cov_tmp.coord;
						start += contigs_info->contig_sz_vt[abs(reads_info->LongReadIndexVec[i][j - 1].contig_no)];
						corr_vec.push_back(Coord_CTG_Cov_tmp);
					}

				}

			}
			corr_vec.push_back(reads_info->LongReadIndexVec[abs(current_read)][j]);

		}
		if (corr_vec.size() > reads_info->LongReadIndexVec[abs(current_read)].size())
		{
			if (Debug)
			{
				debug_out2 << "orig: ";
				for (int j = 0; j < reads_info->LongReadIndexVec[abs(current_read)].size(); ++j)
				{
					debug_out2 << reads_info->LongReadIndexVec[abs(current_read)][j].contig_no << ", ";
				}
				debug_out2 << endl;
				
				debug_out2 << "corr: ";
				for (int j = 0; j < corr_vec.size(); ++j)
				{
					debug_out2 << corr_vec[j].contig_no << ", ";
				}
				debug_out2 << endl;

			}
			reads_info->LongReadIndexVec[abs(current_read)] = corr_vec;
		}


	}
	time(&end_time);
	cout << "Avg alignment size: " << seq_len_sum / 2 / n_alignments << endl;
	if (align_info.SparseAlign)
	cout << "Avg sparse alignment size: " << seq_len_sum_compressed / 2 / n_alignments << endl;
	cout << n_alignments << " alignments calculated." << endl;
	cout << difftime(end_time, beg_time) << " secs." << endl;


}



void LoadNonContainedRawReads(reads_info *reads_info)
{
	cout << "Loading non-contained sequences." << endl;
	int n_loaded = 0;
	reads_info->selected_long_reads_seq.clear();
	reads_info->selected_long_reads_seq.resize(reads_info->LenVec.size());
	map<string, int> tag2id;
	//cout <<"total: "<< reads_info->LenVec.size() << endl;
	for (int i = 1; i < reads_info->LenVec.size(); ++i)
	{
		//cout << i << endl;
		if (reads_info->contained[abs(i)] == 0 || reads_info->chimeric[abs(i)] == 1 || reads_info->LeftOverlaps[i].size()>0 || reads_info->RightOverlaps[i].size()>0)
		{
			reads_info->tag_vec[i][0] = '>';
			tag2id[reads_info->tag_vec[i]] = i;
		}
	}
	//cout << "loading.." << endl;
	for (int i = 0; i < reads_info->input_files.size(); ++i)
	{
		//cout << i << endl;
		bool fq_flag = 0;
		string tag, n_tag, read_str, quality_score;

		ifstream LongReads_in(reads_info->input_files[i].c_str());

		getline(LongReads_in, tag);
		if (tag[0] == '@')
		{
			fq_flag = 1;
		}
		LongReads_in.close();

		LongReads_in.open(reads_info->input_files[i].c_str());

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

			tag[0] = '>';
			if (tag2id.count(tag))
			{
				//deal with 'N's in a stupid way
				reads_info->selected_long_reads_seq[tag2id[tag]] = read_str;
				n_loaded++;
			}


		}

	}
	cout << n_loaded << " loaded." << endl;
}

/*
After this function, 
Left/RightOverlaps stores all overlaps with scores, and is used for all later score look ups
Left/RightOverlapsWithTies store the best left/right overlaps, but with no scores.
*/
void BuildOverlapGraph(reads_info *reads_info, contigs_info *contigs_info, align_info *align_info_default)
{


	align_info align_info;
	align_info = *align_info_default;

	align_matrices align_matrices;

	//align_info.fix_orientation = 1;
	align_info.force_flip = 0;
	align_info.flip = 0;
	/*
	for each read if it is longer than a threshold, then check all the reads that potentially overlap it,
	if the overlap is larger than a threshold and it looks like a prefix-suffix overlap, then build a link.
	Construct the graph from the longest informative reads.
	*/
	bool Debug =reads_info->Debug;
	ofstream Alignment_debug;
	if (Debug)
	{
		Alignment_debug.open("Alignment_debug.txt");
	}
	
	map<int, int> classification_summary;
	int numReads = reads_info->LenVec.size();
	

	for (int i = 0; i < numReads; ++i)
	{
		reads_info->contained[i] = 0;

	}
	

	cout << numReads-1 << " reads." << endl;

	reads_info->LeftOverlaps.clear();
	reads_info->RightOverlaps.clear();
	reads_info->LeftOverlaps.resize(numReads + 1);
	reads_info->RightOverlaps.resize(numReads + 1);
	int initial_round = 1;
	if (reads_info->AdaptiveTh < 0.0)
	{
		//initial_round = 2;
	}
	for (int round = initial_round; round <= 2; ++round)
	{
		int64_t n_alignments = 0;
		int64_t seq_len_sum = 0,seq_len_sum_compressed=0;

		cout << "Calculating reads overlaps, round " << round << endl;
		time_t beg_time, end_time ;
		time(&beg_time);
		if (reads_info->MSA&&round == 1) //&& contigs_info->ContigCovTh>0)
		{
			align_info.CalculateOffset = 1;
			for (int cc = 1; cc <=1; ++cc)
			{

				cout << "Multiple alignment for error correction." << endl;
				//RecoverFalseNegatives(reads_info, contigs_info, &align_info);
				//reads_info->Debug = 0;
				
				MultipleAlignmentErrorCorrection(reads_info, contigs_info, &align_info);
				
				string filename = "CleanedLongReads.txt";
				vector<string> filename_vec;
				filename_vec.push_back(filename);
				reads_info->Clean = 0;
				LoadLongReadIndex(filename_vec, reads_info, contigs_info);//for long reads
				ifstream in_chimera("ChimeraIdx.txt");
				int chimera_idx;
				while (in_chimera >> chimera_idx)
				{
					reads_info->chimeric[chimera_idx] = 1;

					reads_info->contained[chimera_idx] = 1;
				}
				if (reads_info->RecoverFalseNegatives)
				{
					RecoverFalseNegatives(reads_info, contigs_info, &align_info);
				}


			}
			//align_info.scoring_method = 3;
			align_info.CalculateOffset = 0;
			if (reads_info->RemoveChimera)
			{
				ofstream o_chimeras("Chimeras.txt");

				for (int i = 0; i < reads_info->chimeric.size(); ++i)
				{
					if (reads_info->chimeric[i])
					{
						o_chimeras << reads_info->tag_vec[i] << endl;
					}
				}
			}

			cout << "Done." << endl;
			time(&end_time);
			cout << "MSA time: " << difftime(end_time, beg_time) << " secs." << endl;

		}



		
		if (round == 2 )
		{
			if (reads_info->RecoverFalseNegatives)
			{
				//reads_info->Debug = 1;
				//RecoverFalseNegatives(reads_info, contigs_info, &align_info);
			}
			
			if (reads_info->AdaptiveTh>0.0)
			{

				ofstream o_non_contained_reads("CompressedNonContainedReads.txt");

				for (int i = 1; i < reads_info->tag_vec.size(); ++i)
				{
					if (reads_info->chimeric[i] == 0 && reads_info->contained[i] == 0)
					{

						o_non_contained_reads <<">"<< i << endl;
						o_non_contained_reads << reads_info->LenVec[i] << endl;
						for (int v = 0; v != reads_info->LongReadIndexVec[i].size(); ++v)
						{
							struct Coord_CTG_Cov Coord_CTG_Cov_tmp = reads_info->LongReadIndexVec[i][v];
							
							o_non_contained_reads << Coord_CTG_Cov_tmp.coord << ", ";
							o_non_contained_reads << Coord_CTG_Cov_tmp.contig_no << ", ";
							o_non_contained_reads << Coord_CTG_Cov_tmp.cov << ", ";
							o_non_contained_reads << Coord_CTG_Cov_tmp.coord2;

							o_non_contained_reads << endl;

						}
					}
					else
					{
						o_non_contained_reads << ">" << i << endl;
						o_non_contained_reads << reads_info->LenVec[i] << endl;
					}

				}

				if (0)
				{
					ofstream o_non_contained_reads("NonContainedReadNames.txt");

					for (int i = 0; i < reads_info->tag_vec.size(); ++i)
					{
						if (reads_info->chimeric[i] == 0 && reads_info->contained[i] == 0)
						{
							o_non_contained_reads << reads_info->tag_vec[i] << endl;
						}

					}
				}
				
				//if (reads_info->Debug)
	
				

			}
			

		
		}

		numReads = reads_info->LenVec.size() ;

		for (int i = 0; i < numReads; ++i)
		{
			if (round == 2)
			{
				//cout << i << endl;
			}
			if ((i % 1000000 == 0) && (n_alignments>0))
			{
				cout << i << " reads aligned."<<endl;
				cout << "Avg alignment size: " << seq_len_sum / 2 / n_alignments << endl;
				cout << "total alignments: " << n_alignments << endl;
			}
			int current_read = i;// (reads_info->LengthRank[i]) - reads_info->LenVec.begin();

			if (reads_info->contained[current_read])//contained reads are skipped in any round, which may lead to errors
			{
				continue;
			}
		
			vector<Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[abs(current_read)];

			map<int, int> candidate_reads;// left_ext_map, right_ext_map;
			map<int, vector<Coord_CTG_Cov> > candidate_reads_layouts;
			map<int, vector<int> > sorted_candidate_reads;

			for (int c1 = 0; c1 < reads_info->LongReadIndexVec[abs(current_read)].size(); ++c1)
			{
				int ctg_no = reads_info->LongReadIndexVec[abs(current_read)][c1].contig_no;
				
				if (reads_info->Contig2LongRead[abs(ctg_no)].size()<reads_info->max_reads)
				for (int r = 0; r < reads_info->Contig2LongRead[abs(ctg_no)].size(); ++r)
				{

					int read_idx = reads_info->Contig2LongRead[abs(ctg_no)][r];
				
					if (read_idx == current_read)
					{
						continue;
					}
					
					if (reads_info->contained[read_idx]&&(round==1||align_info.FAST))
					{
						continue;//else we use a slower strategy to avoid some errors
					}

					//candidate_reads[read_idx]++;
					candidate_reads[read_idx] += contigs_info->contig_sz_vt[abs(ctg_no)];

					/*
					if (candidate_reads_layouts.count(read_idx) == 0)
					{
						candidate_reads_layouts[read_idx] = reads_info->LongReadIndexVec[read_idx];
					}
					*/
				}

			}
			
			for (map<int, int>::iterator it1 = candidate_reads.begin(); it1 != candidate_reads.end(); ++it1)
			{
				if (it1->second < align_info.min_overlap)
				{
					continue;//stop if the overlap is too small.
				}
				sorted_candidate_reads[it1->second].push_back(it1->first);
			}
			//map < int, vector<Coord_CTG_Cov> > candidate_reads_layouts2 = candidate_reads_layouts;

			int n_reads = 0;
			for (map<int, vector<int> >::reverse_iterator it2 = sorted_candidate_reads.rbegin(); it2 != sorted_candidate_reads.rend(); ++it2)
			{
				
				//vector<int> tmp_vec = it2->second;
				for (int jj = 0; jj < it2->second.size(); ++jj)
				{
					if (n_reads>reads_info->TopOverlaps)
					{
						break;
					}
					candidate_reads_layouts[it2->second[jj]] = reads_info->LongReadIndexVec[it2->second[jj]];// candidate_reads_layouts2[it2->second[jj]];
					n_reads++;
					
				}
			}

			//map<int, int> sorted_left_overlaps, sorted_right_overlaps;
			map<int, vector<int> > sorted_left_overlaps, sorted_right_overlaps;

			for (map<int, vector<Coord_CTG_Cov> >::iterator it1 = candidate_reads_layouts.begin(); it1 != candidate_reads_layouts.end(); ++it1)
			{
				int read_idx = it1->first;
				vector<Coord_CTG_Cov> Coord_CTG_Cov_Vec = it1->second;

				if (reads_info->contained[current_read])
				{
					//cout << ".";
					break;// case 4 in the following happened
				}
				//align_maps align_maps;
				//approximate_semi_global_align_banded(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_maps, &align_info);
				//int result = classify_alignment(&align_info);
				//approximate_semi_local_align_banded(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_maps, &align_info);

				int orientation=0;
				int match_ctgs1 = 0, match_ctgs2 = 0;

				if (align_info.fix_orientation)
				{
					map<int, int> local_index, ref_index, qry_index;
					int ref_sz = current_read_layout.size(), qry_sz = Coord_CTG_Cov_Vec.size();

					for (int i = 0; i < ref_sz; ++i)
					{
						local_index[current_read_layout[i].contig_no] = current_read_layout[i].coord;
					}

					
					for (int i = 0; i < qry_sz; ++i)
					{
						if (local_index.count(Coord_CTG_Cov_Vec[i].contig_no))
						{
							match_ctgs1++;
						}
						if (local_index.count(-Coord_CTG_Cov_Vec[i].contig_no))
						{
							match_ctgs2++;
						}
					}
					/*
					if (match_ctgs1 > match_ctgs2)
					{
						orientation = 0;
					}
					else
					{
						orientation = 1;
					}
					if (match_ctgs1 == match_ctgs2)
					{
						continue;
					}
					*/
				}
				

				for (int orientation = 0; orientation <= 1; ++orientation)
				{
					if (align_info.fix_orientation)
					{

						if (orientation == 0)
						{
							align_info.force_flip = 0;

							if (match_ctgs1 == 0)
							{
								continue;
							}

						}
						else
						{
							align_info.force_flip = 1;

							if (match_ctgs2 == 0)
							{
								continue;
							}
						}
					}
					else
					{
						if (orientation > 0)
						{
							break;
						}
					}

					n_alignments++;
					seq_len_sum += current_read_layout.size();
					seq_len_sum += Coord_CTG_Cov_Vec.size();

					align_info.ref_len = reads_info->LenVec[abs(current_read)];
					align_info.qry_len = reads_info->LenVec[abs(read_idx)];

					if (align_info.SparseAlign)
					{
						sparse_semi_global_align(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_matrices, &align_info);
						seq_len_sum_compressed += align_info.ref_len_compressed;
						seq_len_sum_compressed += align_info.qry_len_compressed;
					}
					else
					{
						approximate_semi_global_align(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_matrices, &align_info);
					}
					int result = classify_alignment(&align_info);

					
					if (Debug)
					{
						//classification_summary[result]++;
						Alignment_debug << "Input pair: " << endl;
						for (int r = 0; r < current_read_layout.size(); ++r)
						{
							Alignment_debug << current_read_layout[r].contig_no << ", ";
						}
						Alignment_debug << endl;
						for (int r = 0; r < Coord_CTG_Cov_Vec.size(); ++r)
						{
							Alignment_debug << Coord_CTG_Cov_Vec[r].contig_no << ", ";
						}
						Alignment_debug << endl;
						Alignment_debug << ">Ref_read_" << current_read << endl;

						for (int r = 0; r < align_info.ref_aligned.size(); ++r)
						{
							Alignment_debug << align_info.ref_aligned[r] << ", ";
						}
						Alignment_debug << endl;
						//Alignment_debug << ">Read_" << current_read << endl;
						for (int r = 0; r < align_info.ref_aligned.size(); ++r)
						{
							Alignment_debug << contigs_info->contig_sz_vt[abs(align_info.ref_aligned[r])] << ", ";
						}
						Alignment_debug << endl;

						Alignment_debug << ">Qry_read_" << read_idx << endl;
						for (int r = 0; r < align_info.qry_aligned.size(); ++r)
						{
							Alignment_debug << align_info.qry_aligned[r] << ", ";
						}
						Alignment_debug << endl;
						//Alignment_debug << ">Read_" << current_read << endl;
						for (int r = 0; r < align_info.qry_aligned.size(); ++r)
						{
							Alignment_debug << contigs_info->contig_sz_vt[abs(align_info.qry_aligned[r])] << ", ";
						}
						Alignment_debug << endl;
						Alignment_debug << "Case: " << result << ", score: " << align_info.max_score << endl;

					}

					//result = classify_alignment(&align_info);




					if (round == 1)
					{

						if (result == 3)
						{
							if (reads_info->contained[current_read] == 0)
							{
								reads_info->contained[abs(read_idx)] = 1;

							}

						}

						if (result == 4)
						{
							if (reads_info->contained[abs(read_idx)] == 0)
							{
								reads_info->contained[abs(current_read)] = 1;

								break;//!!!!!!!!!!

							}
						}
					}
					if (round == 2)
					{

						///newly added.
						if (result == 1)
						{
							int dist1 = align_info.frontal_mismatch_len[0];
							int dist2 = align_info.rear_mismatch_len[1];
							int next_read = read_idx;
							int score = align_info.max_score;
							//int score = align_info.max_len;
							if (align_info.flip)
							{
								if (reads_info->RightOverlaps[next_read].count(-current_read) == 0)
								{
									reads_info->RightOverlaps[next_read][-current_read] = score;// (dist2);
								}
								next_read = -next_read;

							}
							else
							{
								if (reads_info->LeftOverlaps[next_read].count(current_read) == 0)
								{
									reads_info->LeftOverlaps[next_read][current_read] = score;//(dist2);
								}
							}

							if (reads_info->RightOverlaps[current_read].count(next_read) == 0)
							{
								reads_info->RightOverlaps[current_read][next_read] = score;// (dist1);
							}

							//
							//sorted_right_overlaps[score] = next_read;

						}
						if (result == 2)
						{
							int dist1 = align_info.frontal_mismatch_len[1];
							int dist2 = align_info.rear_mismatch_len[0];
							int next_read = read_idx;
							int score = align_info.max_score;
							//int score = align_info.max_len;
							if (align_info.flip)
							{
								if (reads_info->LeftOverlaps[next_read].count(-current_read) == 0)
								{
									reads_info->LeftOverlaps[next_read][-current_read] = score;//(dist1);
								}
								next_read = -next_read;
							}
							else
							{
								if (reads_info->RightOverlaps[next_read].count(current_read) == 0)
								{
									reads_info->RightOverlaps[next_read][current_read] = score;//(dist1);
								}
							}

							if (reads_info->LeftOverlaps[current_read].count(next_read) == 0)
							{
								reads_info->LeftOverlaps[current_read][next_read] = score;// (dist2);
							}

						}


						int next_read = read_idx;
						if (align_info.flip)
						{
							next_read = -next_read;
						}
						if (result == 1)
						{
							int score = align_info.max_score; //align_info.max_len;// 
							//int score = align_info.max_len;// 
							sorted_right_overlaps[score].push_back(next_read);
						}
						if (result == 2)
						{
							int score = align_info.max_score;//align_info.max_len; //
							//int score = align_info.max_len;// 
							sorted_left_overlaps[score].push_back(next_read);

						}

					}

				}

			}
			if (sorted_left_overlaps.size() > 0)
			{
				map<int, vector<int> >::reverse_iterator rit = sorted_left_overlaps.rbegin();
				reads_info->LeftOverlapsBest[current_read] = rit->second[0];



				for (int c = 0; c < rit->second.size(); ++c)
				{
					int next_read = rit->second[c];
					if (!reads_info->contained[abs(next_read)])
					{
						reads_info->LeftOverlapsBestWithTies[current_read].push_back(next_read);

					}

				}

				if (round == 2 && align_info.FAST == 0)// a little buggy
				{



					int non_contained_neighbors = 0;
					for (int c = 0; c < rit->second.size(); ++c)
					{
						if (!reads_info->contained[abs(rit->second[c])])
						{
							non_contained_neighbors++;
						}
					}

					for (int c = 0; c < rit->second.size(); ++c)
					{
						int next_read = rit->second[c];
						if (round == 2 && align_info.FAST == 0 && reads_info->contained[abs(next_read)])// a little buggy
						{
							if (non_contained_neighbors == 0 && !reads_info->chimeric[abs(next_read)])
							{
								reads_info->contained[abs(next_read)] = 0;
								//cout << next_read << ", ";
								reads_info->LeftOverlapsBestWithTies[current_read].push_back(next_read);
							}
						}
					}
				}


			}
			if (sorted_right_overlaps.size() > 0)
			{
				map<int, vector<int> >::reverse_iterator rit = sorted_right_overlaps.rbegin();
				reads_info->RightOverlapsBest[current_read] = rit->second[0];



				for (int c = 0; c < rit->second.size(); ++c)
				{
					int next_read = rit->second[c];
					if (!reads_info->contained[abs(next_read)])
					{
						reads_info->RightOverlapsBestWithTies[current_read].push_back(next_read);

					}

				}

				if (round == 2 && align_info.FAST == 0)// a little buggy
				{
					int non_contained_neighbors = 0;
					for (int c = 0; c < rit->second.size(); ++c)
					{
						if (!reads_info->contained[abs(rit->second[c])])
						{
							non_contained_neighbors++;
						}
					}

					for (int c = 0; c < rit->second.size(); ++c)
					{
						int next_read = rit->second[c];
						if (round == 2 && align_info.FAST == 0 && reads_info->contained[abs(next_read)])// a little buggy
						{
							if (non_contained_neighbors == 0 && !reads_info->chimeric[abs(next_read)])
							{
								//cout << next_read << ", ";
								reads_info->contained[abs(next_read)] = 0;
								reads_info->RightOverlapsBestWithTies[current_read].push_back(next_read);
							}
						}
					}
				}

			}

		}
		time(&end_time);
		if (n_alignments > 0)
		{
			cout << "Avg alignment size: " << seq_len_sum / 2 / n_alignments << endl;
			if (align_info.SparseAlign)
				cout << "Avg sparse alignment size: " << seq_len_sum_compressed / 2 / n_alignments << endl;

		}
		cout << n_alignments << " alignments calculated." << endl;
		cout << "Round " << round << " takes " << difftime(end_time, beg_time) << " secs." << endl;
	}

	/*
	for (map<int, int>::iterator tmp_it = classification_summary.begin(); tmp_it != classification_summary.end(); ++tmp_it)
	{
		o_classification_summary << tmp_it->first << ", " << tmp_it->second << endl;
	}
	*/

	int n_contained = 0;
	for (int i = 1; i < numReads; ++i)
	{
		if (reads_info->contained[i])
		{
			n_contained++;
		}
	}
	cout << n_contained << " contained out of " << numReads-1 << endl;


}


bool read_match_info_cmp(const struct read_match_info &a, const  struct read_match_info &b)
{
	return (a.score) > (b.score);
}


void CollectConsensusInfo(reads_info *reads_info, contigs_info *contigs_info, align_info *align_info_default)
{


	ofstream Alignment_debug, Alignment_debug2;
	bool Debug = 0;// reads_info->Debug;
	if (Debug)
	{
		Alignment_debug.open("ConsensusAlignment_debug.txt");
		Alignment_debug2.open("ConsensusAlignment_debug2.txt");
		
	}



	align_info align_info;
	align_info.force_flip = 0;
	align_info.fix_orientation = 0;
	align_info = *align_info_default;
	align_matrices align_matrices;

	align_info.mm_penalty = 0.3;
	align_info.mm_ratio = 1.0;
	align_info.max_mismatch_FR_score = align_info_default->min_overlap * 5;
	align_info.min_overlap = align_info_default->min_overlap/5;
	align_info.min_extension_FR = align_info_default->min_overlap * 0;
	/*
	for each read if it is longer than a threshold, then check all the reads that potentially overlap it,
	if the overlap is larger than a threshold and it looks like a prefix-suffix overlap, then build a link.
	Construct the graph from the longest informative reads.
	*/

	int numReads = reads_info->LenVec.size() ;
	
	cout <<"Collecting information for consensus."<<endl<< numReads-1 << " reads." << endl;
	
	
	int64_t n_alignments = 0;
	int64_t seq_len_sum = 0, seq_len_sum_compressed=0;

	cout << "Calculating reads overlaps." << endl;
	time_t beg_time, end_time;
	time(&beg_time);
		
	align_info.CalculateOffset = 0;
		
	for (int i = 0; i < numReads; ++i)
	{

		if ((i % 1000000 == 0) && (n_alignments>0))
		{
			cout << i << " reads aligned." << endl;
			cout << "Avg alignment size: " << seq_len_sum / 2 / n_alignments << endl;
			if (align_info.SparseAlign)
			{
				cout << "Avg sparse alignment size: " << seq_len_sum_compressed / 2 / n_alignments << endl;

			}
			cout << "total alignments: " << n_alignments << endl;
		}
		int current_read = i;// (reads_info->LengthRank[i]) - reads_info->LenVec.begin();

		if (reads_info->contained[current_read])//contained reads are skipped in any round, which may lead to errors
		{
			continue;
		}

		vector<Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[abs(current_read)];

		map<int, int> candidate_reads;//, left_ext_map, right_ext_map;
		map<int, vector<Coord_CTG_Cov> > candidate_reads_layouts;
		map<int, vector<int> > sorted_candidate_reads;

		for (int c1 = 0; c1 < reads_info->LongReadIndexVec[abs(current_read)].size(); ++c1)
		{
			int ctg_no = reads_info->LongReadIndexVec[abs(current_read)][c1].contig_no;

			if (reads_info->Contig2LongRead[abs(ctg_no)].size() < reads_info->max_reads*100)
			for (int r = 0; r < reads_info->Contig2LongRead[abs(ctg_no)].size(); ++r)
			{

				int read_idx = reads_info->Contig2LongRead[abs(ctg_no)][r];

				if (read_idx == current_read)
				{
					continue;
				}
				if (align_info.matching_method == 1)
				{
					candidate_reads[read_idx] += contigs_info->contig_sz_vt[abs(ctg_no)];
				}
				if (align_info.matching_method == 2 || align_info.matching_method == 3)
				{
					candidate_reads[read_idx] += reads_info->LongReadIndexVec[abs(current_read)][c1].cov;
				}
			}

		}

		for (map<int, int>::iterator it1 = candidate_reads.begin(); it1 != candidate_reads.end(); ++it1)
		{
			if (it1->second < align_info.min_overlap)
			{
				continue;//stop if the overlap is too small.
			}
			sorted_candidate_reads[it1->second].push_back(it1->first);
		}

		int n_reads = 0;
		for (map<int, vector<int> >::reverse_iterator it2 = sorted_candidate_reads.rbegin(); it2 != sorted_candidate_reads.rend(); ++it2)
		{

			//vector<int> tmp_vec = it2->second;
			for (int jj = 0; jj < it2->second.size(); ++jj)
			{
				if (n_reads>10000)
				{
					break;
				}
				candidate_reads_layouts[it2->second[jj]] = reads_info->LongReadIndexVec[it2->second[jj]];// candidate_reads_layouts2[it2->second[jj]];
				n_reads++;

			}
		}

		//map<int, int> sorted_left_overlaps, sorted_right_overlaps;
		map<int, vector<int> > sorted_left_overlaps, sorted_right_overlaps;

		for (map<int, vector<Coord_CTG_Cov> >::iterator it1 = candidate_reads_layouts.begin(); it1 != candidate_reads_layouts.end(); ++it1)
		{
			int read_idx = it1->first;
			vector<Coord_CTG_Cov> Coord_CTG_Cov_Vec = it1->second;

			n_alignments++;
			seq_len_sum += current_read_layout.size();
			seq_len_sum += Coord_CTG_Cov_Vec.size();

			align_info.ref_len = reads_info->LenVec[abs(current_read)];
			align_info.qry_len = reads_info->LenVec[abs(read_idx)];

			if (align_info.SparseAlign)
			{
				sparse_semi_global_align(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_matrices, &align_info);
				seq_len_sum_compressed += align_info.ref_len_compressed;
				seq_len_sum_compressed += align_info.qry_len_compressed;
			}
			else
			{
				approximate_semi_global_align(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_matrices, &align_info);

			}
			int result = classify_alignment(&align_info);


			if (result >= 1 && result <= 4)
			{
				
				if (1)//reads_info->Consensus_info_full.count(abs(read_idx)) == 0)//|| (reads_info->Consensus_info_full[abs(read_idx)].score<align_info.max_score)
				{
					//reads_info->Consensus_info[abs(current_read)].push_back(read_idx);
					read_match_info read_match_info;
					read_match_info.read_idx = abs(current_read);
					read_match_info.score = align_info.max_score;
					reads_info->Consensus_info_full[abs(read_idx)].push_back(read_match_info);

				}
			
			}

		
			if (Debug)
			{
				Alignment_debug << ">Ref_read_" << current_read << endl;
				
				for (int r = 0; r < current_read_layout.size(); ++r)
				{
					Alignment_debug << current_read_layout[r].contig_no << ", ";
				}
				Alignment_debug << endl;
				
				for (int r = 0; r < Coord_CTG_Cov_Vec.size(); ++r)
				{
					Alignment_debug << Coord_CTG_Cov_Vec[r].contig_no << ", ";
				}
				Alignment_debug << endl;

			}
		}
	}


	map<int, vector<read_match_info> >::iterator it;
	for (it = reads_info->Consensus_info_full.begin(); it != reads_info->Consensus_info_full.end(); ++it)
	{
		if (it->second.size() > 1)
		{
			sort(it->second.begin(), it->second.end(), read_match_info_cmp);
			for (int ii = 0; ii < reads_info->ConsensusBestN; ++ii)
			{
				if (ii < it->second.size())
				{
					reads_info->Consensus_info[abs((int)it->second[ii].read_idx)].push_back(it->first);
				}
			}
			
		}
	}

	

	time(&end_time);
	if (n_alignments>0)
	{
		cout << "Avg alignment size: " << seq_len_sum / 2 / n_alignments << endl;
		if (align_info.SparseAlign)
		{
			cout << "Avg sparse alignment size: " << seq_len_sum_compressed / 2 / n_alignments << endl;
		}
	}
	cout << n_alignments << " alignments calculated." << endl;
	cout << difftime(end_time, beg_time) << " secs." << endl;

}


void Create_Local_Contig_Kmer_Index(contigs_info *contigs_info, align_info* align_info_default, int block_sz, int K_size, string read1, string read2, string &backbone_raw, vector<Coord_CTG_Cov> compressed_read1, vector<Coord_CTG_Cov> compressed_read2, ref_read_t &long_read)
{
	bool Debug = 0;
	Debug = align_info_default->Debug;
	//ofstream o_debug("Ext_debug.txt",ios_base::app);
	ref_read_t ref;
	align_info align_info0 = *align_info_default;
	align_matrices align_matrices;
	int max_ctg_len = 0, len_sum = 0;
	map<int, int> shared_ctgs, shared_ctgs2;

	map<int, int>::iterator map_it;
	vector<Coord_CTG_Cov> c_read1 = compressed_read1;
	vector<Coord_CTG_Cov> c_read2 = compressed_read2;
	/*
	cout << "read1_raw: " << endl;
	for (int i = 0; i < compressed_read1.size(); ++i)
	{
		cout << "ctg_" << compressed_read1[i].contig_no << ":cov_" << compressed_read1[i].cov << endl;
	}

	cout << "read2_raw: " << endl;
	for (int i = 0; i < compressed_read2.size(); ++i)
	{
		cout << "ctg_"<<compressed_read2[i].contig_no << ":cov_" << compressed_read2[i].cov << endl;
	}
	*/
	for (int i = 0; i < compressed_read1.size(); ++i)
	{
		//shared_ctgs[abs(compressed_read1[i].contig_no)] = 1;
		shared_ctgs[compressed_read1[i].contig_no] = 1;

	}
	for (int i = 0; i < compressed_read2.size(); ++i)
	{
		//if (shared_ctgs.count(abs(compressed_read2[i].contig_no)))
		if (shared_ctgs.count(compressed_read2[i].contig_no))
		{
			//shared_ctgs[abs(compressed_read2[i].contig_no)] = 2;
			shared_ctgs[compressed_read2[i].contig_no] = 2;
		}
	}

	for (map_it = shared_ctgs.begin(); map_it != shared_ctgs.end(); ++map_it)
	{
		if (map_it->second == 2)
		{
			if (contigs_info->contig_sz_vt[abs(map_it->first)] > max_ctg_len)
			{
				max_ctg_len = contigs_info->contig_sz_vt[abs(map_it->first)];
			}
			len_sum += contigs_info->contig_sz_vt[abs(map_it->first)];
		}

	}

	ref.read_bits = (uint64_t *)myMalloc((size_t)(max_ctg_len / 4) + 100);
	ref.alloc_sz = (size_t)(max_ctg_len / 4 + 100);
	struct hashtable ht;
	Init_HT(&ht, len_sum + 1);
	for (map_it = shared_ctgs.begin(); map_it != shared_ctgs.end(); ++map_it)
	{
		if (map_it->second == 2)
		{
			shared_ctgs2[abs(map_it->first)] = 2;
		}
	}
	for (int round = 1; round <= 2; ++round)
	{
		ht.round = round;
		int bucket_count = 0;
		for (map_it = shared_ctgs2.begin(); map_it != shared_ctgs2.end(); ++map_it)
		{
			if (map_it->second == 2)
			{
				string contig_str = contigs_info->contigs_str[abs(map_it->first)];

				Init_Ref_Read(contig_str, ref);
				ref.contig_no = abs(map_it->first);
				Contig_Kmer_Index(&ref, &ht, K_size, &bucket_count);
			}

		}

		if (ht.round == 1)
		{
			SwitchBuckets(&ht);
		}

	}




	Init_Ref_Read(read1, long_read);
	long_read.read_idx = 1;
	struct LongReadContigIndex LongReadContigIndex1, LongReadContigIndex2;
	LongReadContigIndex1.KmerCovTh = contigs_info->KmerCovTh;
	LongReadContigIndex1.nMatches = 0;
	LongReadContigIndex1.BlockSize = block_sz;
	LongReadContigIndex1.FastMap = 0;
	LongReadContigIndex2.KmerCovTh = contigs_info->KmerCovTh;
	LongReadContigIndex2.nMatches = 0;
	LongReadContigIndex2.BlockSize = block_sz;
	LongReadContigIndex2.FastMap = 0;
	ReadCompression(&long_read, &ht, contigs_info, K_size, &LongReadContigIndex1);
	Init_Ref_Read(read2, long_read);
	long_read.read_idx = 2;
	ReadCompression(&long_read, &ht, contigs_info, K_size, &LongReadContigIndex2);




	align_info0.ref_len = read1.size();
	align_info0.qry_len = read2.size();
	align_info0.fix_orientation = 1;
	align_info0.force_flip = 0;
	align_info0.CalculateOffset = 1;
	align_info0.flip = 0;

	compressed_read1.clear();
	compressed_read2.clear();

	map<int, ContigInRead>::iterator it;
	//o_debug << "r1:" << endl;
	for (it = LongReadContigIndex1.layout.begin(); it != LongReadContigIndex1.layout.end(); ++it)
	{
		Coord_CTG_Cov Coord_CTG_Cov_tmp;
		Coord_CTG_Cov_tmp.coord = it->first;
		Coord_CTG_Cov_tmp.contig_no = it->second.ctg_no;
		Coord_CTG_Cov_tmp.cov = it->second.cov;
		Coord_CTG_Cov_tmp.coord2 = it->second.coord2;
		//cout << it->second.ctg_no<<" ";
		//if (shared_ctgs[abs(Coord_CTG_Cov_tmp.contig_no)] == 2)
		if (shared_ctgs[Coord_CTG_Cov_tmp.contig_no] == 2)
		{
			compressed_read1.push_back(Coord_CTG_Cov_tmp);
		//	o_debug << Coord_CTG_Cov_tmp.contig_no << ", ";
		}
		//o_debug << endl;
		
	}
	//cout << endl;
	//o_debug << "r2:" << endl;
	for (it = LongReadContigIndex2.layout.begin(); it != LongReadContigIndex2.layout.end(); ++it)
	{
		Coord_CTG_Cov Coord_CTG_Cov_tmp;
		Coord_CTG_Cov_tmp.coord = it->first;
		Coord_CTG_Cov_tmp.contig_no = it->second.ctg_no;
		Coord_CTG_Cov_tmp.cov = it->second.cov;
		Coord_CTG_Cov_tmp.coord2 = it->second.coord2;
		//if (shared_ctgs[abs(Coord_CTG_Cov_tmp.contig_no)] == 2)
		if (shared_ctgs[Coord_CTG_Cov_tmp.contig_no] == 2)
		{
			compressed_read2.push_back(Coord_CTG_Cov_tmp);
		//	o_debug << Coord_CTG_Cov_tmp.contig_no << ", ";
		}
		//o_debug << endl;
		//cout << it->second.ctg_no << " ";
	}

//	cout << endl;


	/*

	cout << "read1_pp: " << endl;
	for (int i = 0; i < compressed_read1.size(); ++i)
	{
		cout << "ctg_" << compressed_read1[i].contig_no << ":cov_" << compressed_read1[i].cov << endl;
	}

	cout << "read2_pp: " << endl;
	for (int i = 0; i < compressed_read2.size(); ++i)
	{
		cout << "ctg_" << compressed_read2[i].contig_no << ":cov_" << compressed_read2[i].cov << endl;
	}

	*/


	if (align_info0.SparseAlign)
	{
		sparse_semi_global_align(contigs_info, compressed_read1, compressed_read2, &align_matrices, &align_info0);
	}
	else
	{
		approximate_semi_global_align(contigs_info, compressed_read1, compressed_read2, &align_matrices, &align_info0);
	}
	
	int result = classify_alignment(&align_info0);
	if (Debug)
	{
		ofstream o_debug("Problematic_reads.txt");
		o_debug << ">read1" << endl << read1 << endl;
		o_debug << ">read2" << endl << read2 << endl;
		for (map_it = shared_ctgs.begin(); map_it != shared_ctgs.end(); ++map_it)
		{
			if (map_it->second == 2)
			{
				string contig_str = contigs_info->contigs_str[abs(map_it->first)];

				o_debug << ">Contig_" << abs(map_it->first)<< endl << contig_str << endl;
			}

		}

		cout << "result=" << result << endl;
		for (int a = 0; a < align_info0.ref_aligned.size(); ++a)
		{
			cout << align_info0.ref_aligned[a] << " ";
		}
		cout <<endl;
		for (int a = 0; a < align_info0.qry_aligned.size(); ++a)
		{
			cout << align_info0.qry_aligned[a] << " ";
		}
		cout << endl;

		for (int a = 0; a < c_read1.size(); ++a)
		{
			cout << c_read1[a].contig_no << " ";
		}
		cout << endl;
		for (int a = 0; a < c_read2.size(); ++a)
		{
			cout << c_read2[a].contig_no << " ";
		}
		cout << endl;
	}
	result = 1;

	//find the center position matching the contig of interest, and align the contig at that place
	int aligned_contig = 0;
	for (int i = 0; i<compressed_read1.size(); ++i)
	{
		if (compressed_read1[i].coord2 == align_info0.ref_aligned_center)
		{
			aligned_contig = compressed_read1[i].contig_no;
		}
	}
	int ReadCoord1 = LongReadContigIndex1.CTG2LR[aligned_contig][LongReadContigIndex1.CTG2LR[aligned_contig].size() / 2];
	int ContigCoord1 = LongReadContigIndex1.CTG2Coords[aligned_contig][ReadCoord1];
	int mean_offset1;
	if (LongReadContigIndex1.CTG2LR_2[aligned_contig].size() < 30)
	{
		mean_offset1 = LongReadContigIndex1.CTG2LR_2[aligned_contig][LongReadContigIndex1.CTG2LR_2[aligned_contig].size() / 2];
	}
	else
	{
		mean_offset1 = LongReadContigIndex1.CTG2LR_2[aligned_contig][LongReadContigIndex1.CTG2LR_2[aligned_contig].size() *2 /5];
	}
	for (int ii = 0; ii < LongReadContigIndex1.CTG2LR[aligned_contig].size(); ++ii)
	{
		ReadCoord1 = LongReadContigIndex1.CTG2LR[aligned_contig][ii];
		ContigCoord1 = LongReadContigIndex1.CTG2Coords[aligned_contig][ReadCoord1];
		int ContigOffset1 = LongReadContigIndex1.CTG2Offsets[aligned_contig][ReadCoord1];
		if (ContigOffset1 == mean_offset1)
		{
			break;
		}
	}

	int ReadCoord2 = LongReadContigIndex2.CTG2LR[aligned_contig][LongReadContigIndex2.CTG2LR[aligned_contig].size() / 2];
	int ContigCoord2 = LongReadContigIndex2.CTG2Coords[aligned_contig][ReadCoord2];
	int mean_offset2;
	if (LongReadContigIndex1.CTG2LR_2[aligned_contig].size() < 30)
	{
		mean_offset2 = LongReadContigIndex2.CTG2LR_2[aligned_contig][LongReadContigIndex2.CTG2LR_2[aligned_contig].size() / 2];
	}
	else
	{
		mean_offset2 = LongReadContigIndex2.CTG2LR_2[aligned_contig][LongReadContigIndex2.CTG2LR_2[aligned_contig].size() * 3/5];
	}
		
		
	//for (int ii = (int)LongReadContigIndex2.CTG2LR[aligned_contig].size() -1; ii >=0; --ii)
	for (int ii = 0; ii < LongReadContigIndex2.CTG2LR[aligned_contig].size(); ++ii)
	{
		ReadCoord2 = LongReadContigIndex2.CTG2LR[aligned_contig][ii];
		ContigCoord2 = LongReadContigIndex2.CTG2Coords[aligned_contig][ReadCoord2];
		int ContigOffset2 = LongReadContigIndex2.CTG2Offsets[aligned_contig][ReadCoord2];
		if (ContigOffset2 == mean_offset2)
		{
			break;
		}
	}

	string contig_str = contigs_info->contigs_str[abs(aligned_contig)];
	if (aligned_contig < 0)
	{
		reverse_complement_str(contig_str);
		ContigCoord1 = contigs_info->contig_sz_vt[abs(aligned_contig)] - ContigCoord1;
		ContigCoord2 = contigs_info->contig_sz_vt[abs(aligned_contig)] - ContigCoord2;

	}
	
	
	if (ContigCoord2 > ContigCoord1)
	{
		int crop_sz = read1.size() - ReadCoord1;
		backbone_raw = backbone_raw.substr(0, backbone_raw.size() - crop_sz);
		backbone_raw += contig_str.substr(ContigCoord1, ContigCoord2 - ContigCoord1);
		
		backbone_raw += read2.substr(ReadCoord2, read2.size());

	}
	else
	{
		int crop_sz = read1.size() - ReadCoord1 + (ContigCoord1 - ContigCoord2);// this is introducing some stupid errors :(
		if (backbone_raw.size() > crop_sz)
		{
			backbone_raw = backbone_raw.substr(0, backbone_raw.size() - crop_sz);
		}
		else
		{
			backbone_raw.clear();
			cout << "Extension warning." << endl;
		}
		backbone_raw += read2.substr(ReadCoord2, read2.size());
	}
		
	
	

	if (ht.ht_sz > 0)
	{
		FreeSparseKmerGraph(&ht);

	}
}


void ConstructBackbone(reads_info *reads_info, contigs_info *contigs_info, align_info *align_info)
{
	ifstream in_layout("layout_reads.txt");
	ifstream in_layout_compressed("layout_compressed_reads.txt");

	ofstream o_backbone_raw("backbone_raw.fasta");
	//ofstream o_backbone_corr("backbone_corr.fasta");
	ofstream o_layout_full("full_layout.fasta");
	ofstream o_layout_coords("full_layout_coords.txt");
	ofstream ext_log("extension_log.txt");
	ofstream ext_warning("extension_warning.txt");
	struct align_info  align_info0 = *align_info;

	bool Debug = reads_info->Debug;
	Debug = 1;
	
	if (Debug)
	{
	//	ofstream o_debug("Ext_debug.txt");
	//	o_debug.close();

	}

	align_matrices align_matrices;

	ref_read_t long_read;
	size_t max_readlen = 10000000;
	long_read.read_bits = (uint64_t *)myMalloc((size_t)(max_readlen / 4) + 100);
	long_read.alloc_sz = (size_t)(max_readlen / 4 + 100);



	bool success = 1;
	string backbone_tag, tag, n_tag, tag_temp;
	getline(in_layout, backbone_tag);
	getline(in_layout_compressed, backbone_tag);
	if (backbone_tag.size() > 0)
	{
		success = 1;
	}
	else
	{
		success = 0;
	}
	int n_backbone = 0;
	while (success)
	{
		n_backbone++;
		vector<Coord_CTG_Cov> TempIndex;
		vector< vector<Coord_CTG_Cov> > compressed_read_vec;
		vector< int > LenVec;
		vector<string> seq_vec;
		vector<string> seq_tag;
		int n_reads = 0;
		string temp;
		getline(in_layout, temp);
		getline(in_layout_compressed, temp);
		if (temp.size() == 0)
		{
			success = 0;
			break;
		}
		n_reads = atoi(temp.c_str());
		tag.clear();
		n_tag.clear();
		for (int i = 0; i < n_reads; ++i)
		{
			int Len;
			string seq;
			get_a_contig_path(in_layout_compressed, tag_temp, TempIndex, Len, 0, tag_temp);
			get_a_fasta_read(in_layout, tag, seq, n_tag);

			string temp_str, raw_read = seq;
			temp_str.resize(seq.size());
			size_t nn = 0;
			for (int n = 0; n < seq.size(); ++n)
			{
				if (seq[n] != 'N')
				{
					temp_str[nn] = seq[n];
					nn++;
				}
			}
			temp_str.resize(nn);
			seq = temp_str;


			Len = seq.size();
			if (Len == 0)
			{
				cout << "Empty sequence loaded. It looks like you have messed up the data." << endl;
				return;
			}
			compressed_read_vec.push_back(TempIndex);
			seq_tag.push_back(tag);
			LenVec.push_back(Len);
			seq_vec.push_back(seq);
		}
		string backbone_raw, backbone_corr,backbone_raw_new;
		o_layout_coords << ">Backbone_" << n_backbone << endl;

		
		//if (n_backbone>=5084)
		for (int r = 0; r < n_reads; ++r)
		{

			//if (r>=43)
			if (r > 0)
			{

				//cout <<"Read:"<< r << endl;

			
				Create_Local_Contig_Kmer_Index(contigs_info, align_info, 500, reads_info->K_size, seq_vec[r - 1], seq_vec[r], backbone_raw_new, compressed_read_vec[r - 1], compressed_read_vec[r], long_read);

				align_info0.ref_len = LenVec[r - 1];
				align_info0.qry_len = LenVec[r];
				align_info0.CalculateOffset = 1;
				align_info0.fix_orientation = 1;
				align_info0.force_flip = 0;
				if (align_info0.SparseAlign)
				{
					sparse_semi_global_align(contigs_info, compressed_read_vec[r - 1], compressed_read_vec[r], &align_matrices, &align_info0);
				}
				else
				{
					approximate_semi_global_align(contigs_info, compressed_read_vec[r - 1], compressed_read_vec[r], &align_matrices, &align_info0);
				}


				
				int result = classify_alignment(&align_info0);


				
				if (result != 1)
				{
					ext_warning << "Alignment warning: " << result << endl;
					ext_warning << "previous:" << endl;
					for (int l = 0; l < compressed_read_vec[r - 1].size(); ++l)
					{
						ext_warning << compressed_read_vec[r - 1][l].coord << ", ";
						ext_warning << compressed_read_vec[r - 1][l].contig_no << ", ";
						ext_warning << compressed_read_vec[r - 1][l].cov << ", ";
						ext_warning << compressed_read_vec[r - 1][l].coord2 << ", ";
						ext_warning << endl;

					}
					ext_warning << "next:" << endl;
					for (int l = 0; l < compressed_read_vec[r].size(); ++l)
					{
						ext_warning << compressed_read_vec[r][l].coord << ", ";
						ext_warning << compressed_read_vec[r][l].contig_no << ", ";
						ext_warning << compressed_read_vec[r][l].cov << ", ";
						ext_warning << compressed_read_vec[r][l].coord2 << ", ";
						ext_warning << endl;

					}

				}
				result = 1;
				


			}
			else
			{
				//backbone_raw = seq_vec[0];
				backbone_raw_new = seq_vec[0];
				//backbone_corr = seq_vec[0];
				
			}
			o_layout_full << seq_tag[r] << endl;
			o_layout_full << seq_vec[r] << endl;

			if (r + 1 < seq_tag.size())
			{
				string temp_str = seq_tag[r];
				if (temp_str[0] == '>' || temp_str[0] == '@')
				{
					temp_str = temp_str.substr(1, temp_str.size());
				}
				o_layout_coords << temp_str << " ";
			}
			o_layout_coords << endl;
		}

		o_backbone_raw << ">Backbone_" << n_backbone << endl;
		o_backbone_raw << backbone_raw_new << endl;
		//o_backbone_corr << ">Backbone_" << n_backbone << endl;
		//o_backbone_corr << backbone_corr << endl;

	}

}


void CloseGapsWithAllReads(reads_info *reads_info, contigs_info *contigs_info, align_info *align_info_default)
{

	align_info align_info;
	align_info = *align_info_default;

	align_matrices align_matrices;

	align_info.fix_orientation = 1;
	align_info.force_flip = 0;
	align_info.flip = 0;

	bool Debug = reads_info->Debug;
	ofstream Alignment_debug;
	if (Debug)
	{
		Alignment_debug.open("CloseGaps_debug.txt");
	}

	size_t numReads = reads_info->LenVec.size() ;

	

	reads_info->gap_closing.resize(numReads );
	int contig_ends = 0;
	int Chimeras = 0;
	map<int, bool> LeftEnd, RightEnd;
	for (int i = 0; i < numReads; ++i)
	{
		
		reads_info->gap_closing[i] = 0;
		if (reads_info->LeftOverlaps[i].size() != 0 || reads_info->RightOverlaps[i].size() !=0) 
		if (reads_info ->contained[i]==0&&(reads_info->LeftOverlaps[i].size() == 0 || reads_info->RightOverlaps[i].size() == 0))
		{
			
			
			if (reads_info->LeftOverlaps[i].size() == 0)
			{
				LeftEnd[i] = 1;
			}
			if (reads_info->RightOverlaps[i].size() == 0)
			{
				RightEnd[i] = 1;
			}
			reads_info->gap_closing[i] = 1;
			contig_ends++;
		}
		if (reads_info->chimeric[i])
		{
			Chimeras++;
			reads_info->gap_closing[i] = 1;
			reads_info->contained[i] = 0;
			//reads_info->chimeric[i] = 0;
		}
	}




	int initial_round = 1;
	cout << "Closing gaps using " << Chimeras << " reads." << endl;
	cout << "Contig ends " << contig_ends << "." << endl;

	
	for (int round = initial_round; round <= 2; ++round)
	{
		int64_t n_alignments = 0;
		int64_t seq_len_sum = 0, seq_len_sum_compressed=0;

		time_t beg_time, end_time;
		time(&beg_time);
		
		for (int i = 0; i < numReads; ++i)
		{
			
			int current_read = i;// (reads_info->LengthRank[i]) - reads_info->LenVec.begin();

			if (!reads_info->gap_closing[current_read])//contained reads are skipped in any round, which may lead to errors
			{
				continue;
			}

			vector<Coord_CTG_Cov> current_read_layout = reads_info->LongReadIndexVec[abs(current_read)];

			map<int, int> candidate_reads;// left_ext_map, right_ext_map;
			map<int, vector<Coord_CTG_Cov> > candidate_reads_layouts;
			map<int, vector<int> > sorted_candidate_reads;

			for (int c1 = 0; c1 < reads_info->LongReadIndexVec[abs(current_read)].size(); ++c1)
			{
				int ctg_no = reads_info->LongReadIndexVec[abs(current_read)][c1].contig_no;

				if (reads_info->Contig2LongRead[abs(ctg_no)].size()<reads_info->max_reads)
				for (int r = 0; r < reads_info->Contig2LongRead[abs(ctg_no)].size(); ++r)
				{

					int read_idx = reads_info->Contig2LongRead[abs(ctg_no)][r];

					if (read_idx == current_read)
					{
						continue;
					}

					if (!reads_info->gap_closing[read_idx])
					{
						continue;
					}

					//candidate_reads[read_idx]++;
					candidate_reads[read_idx] += contigs_info->contig_sz_vt[abs(ctg_no)];

				}

			}

			for (map<int, int>::iterator it1 = candidate_reads.begin(); it1 != candidate_reads.end(); ++it1)
			{
				if (it1->second < align_info.min_overlap)
				{
					continue;//stop if the overlap is too small.
				}
				sorted_candidate_reads[it1->second].push_back(it1->first);
			}
			//map < int, vector<Coord_CTG_Cov> > candidate_reads_layouts2 = candidate_reads_layouts;

			int n_reads = 0;
			for (map<int, vector<int> >::reverse_iterator it2 = sorted_candidate_reads.rbegin(); it2 != sorted_candidate_reads.rend(); ++it2)
			{

				//vector<int> tmp_vec = it2->second;
				for (int jj = 0; jj < it2->second.size(); ++jj)
				{
					if (n_reads>reads_info->TopOverlaps)
					{
						break;
					}
					candidate_reads_layouts[it2->second[jj]] = reads_info->LongReadIndexVec[it2->second[jj]];// candidate_reads_layouts2[it2->second[jj]];
					n_reads++;

				}
			}

			//map<int, int> sorted_left_overlaps, sorted_right_overlaps;
			map<int, vector<int> > sorted_left_overlaps, sorted_right_overlaps;

			for (map<int, vector<Coord_CTG_Cov> >::iterator it1 = candidate_reads_layouts.begin(); it1 != candidate_reads_layouts.end(); ++it1)
			{
				int read_idx = it1->first;
				vector<Coord_CTG_Cov> Coord_CTG_Cov_Vec = it1->second;

				if (reads_info->contained[current_read])
				{
					//cout << ".";
					break;// case 4 in the following happened
				}
				if (reads_info->contained[read_idx] && (!reads_info->chimeric[read_idx]))
				{
					cout << "bug.";
					break;// case 4 in the following happened
				}
				//align_maps align_maps;
				//approximate_semi_global_align_banded(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_maps, &align_info);
				//int result = classify_alignment(&align_info);
				//approximate_semi_local_align_banded(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_maps, &align_info);

				int orientation = 0;
				int match_ctgs1 = 0, match_ctgs2 = 0;

				if (align_info.fix_orientation)
				{
					map<int, int> local_index, ref_index, qry_index;
					int ref_sz = current_read_layout.size(), qry_sz = Coord_CTG_Cov_Vec.size();

					for (int i = 0; i < ref_sz; ++i)
					{
						local_index[current_read_layout[i].contig_no] = current_read_layout[i].coord;
					}


					for (int i = 0; i < qry_sz; ++i)
					{
						if (local_index.count(Coord_CTG_Cov_Vec[i].contig_no))
						{
							match_ctgs1++;
						}
						if (local_index.count(-Coord_CTG_Cov_Vec[i].contig_no))
						{
							match_ctgs2++;
						}
					}
	
				}

				if (1)
				for (int orientation = 0; orientation <= 1; ++orientation)
				{
					if (align_info.fix_orientation)
					{

						if (orientation == 0)
						{
							align_info.force_flip = 0;

							if (match_ctgs1 == 0)
							{
								continue;
							}

						}
						else
						{
							align_info.force_flip = 1;

							if (match_ctgs2 == 0)
							{
								continue;
							}
						}
					}

					n_alignments++;
					seq_len_sum += current_read_layout.size();
					seq_len_sum += Coord_CTG_Cov_Vec.size();

					align_info.ref_len = reads_info->LenVec[abs(current_read)];
					align_info.qry_len = reads_info->LenVec[abs(read_idx)];


					if (align_info.SparseAlign)
					{
						sparse_semi_global_align(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_matrices, &align_info);
						seq_len_sum_compressed += align_info.ref_len_compressed;
						seq_len_sum_compressed += align_info.qry_len_compressed;
					}
					else
					{
						approximate_semi_global_align(contigs_info, current_read_layout, Coord_CTG_Cov_Vec, &align_matrices, &align_info);
					}

					int result = classify_alignment(&align_info);


					if (Debug)
					{
						//classification_summary[result]++;
						Alignment_debug << "Input pair: " << endl;
						for (int r = 0; r < current_read_layout.size(); ++r)
						{
							Alignment_debug << current_read_layout[r].contig_no << ", ";
						}
						Alignment_debug << endl;
						for (int r = 0; r < Coord_CTG_Cov_Vec.size(); ++r)
						{
							Alignment_debug << Coord_CTG_Cov_Vec[r].contig_no << ", ";
						}
						Alignment_debug << endl;
						Alignment_debug << ">Ref_read_" << current_read << endl;

						for (int r = 0; r < align_info.ref_aligned.size(); ++r)
						{
							Alignment_debug << align_info.ref_aligned[r] << ", ";
						}
						Alignment_debug << endl;
						//Alignment_debug << ">Read_" << current_read << endl;
						for (int r = 0; r < align_info.ref_aligned.size(); ++r)
						{
							Alignment_debug << contigs_info->contig_sz_vt[abs(align_info.ref_aligned[r])] << ", ";
						}
						Alignment_debug << endl;

						Alignment_debug << ">Qry_read_" << read_idx << endl;
						for (int r = 0; r < align_info.qry_aligned.size(); ++r)
						{
							Alignment_debug << align_info.qry_aligned[r] << ", ";
						}
						Alignment_debug << endl;
						//Alignment_debug << ">Read_" << current_read << endl;
						for (int r = 0; r < align_info.qry_aligned.size(); ++r)
						{
							Alignment_debug << contigs_info->contig_sz_vt[abs(align_info.qry_aligned[r])] << ", ";
						}
						Alignment_debug << endl;
						Alignment_debug << "Case: " << result << ", score: " << align_info.max_score << endl;

					}

					//result = classify_alignment(&align_info);

					if (round == 1)
					{

						if (result == 3)
						{
							if (reads_info->contained[current_read] == 0 && reads_info->chimeric[abs(read_idx)])//tolerate the contig ends
							{
								reads_info->contained[abs(read_idx)] = 1;

							}

						}

						if (result == 4 && reads_info->chimeric[abs(current_read)])
						{
							if (reads_info->contained[abs(read_idx)] == 0)
							{
								reads_info->contained[abs(current_read)] = 1;

								break;//!!!!!!!!!!

							}
						}
					}
					
						

					if (round == 2)
					{
						if (result == 3)
						{
							if (!reads_info->chimeric[abs(read_idx)])
							{
								if (align_info.flip)
								{
									if (LeftEnd.count(abs(read_idx)))
									{
										result = 2;
									}
									if (RightEnd.count(abs(read_idx)))
									{
										result = 1;
									}
								}
								else
								{
									if (LeftEnd.count(abs(read_idx)))
									{
										result = 1;
									}
									if (RightEnd.count(abs(read_idx)))
									{
										result = 2;
									}
								}
							}
						}
						if (result == 4)
						{
							if (!reads_info->chimeric[abs(current_read)])
							{

								if (LeftEnd.count(abs(current_read)))
								{
									result = 2;
								}
								if (RightEnd.count(abs(current_read)))
								{
									result = 1;
								}
							}

						}

					
						
					

						//
						if (result == 1)
						{
							int dist1 = align_info.frontal_mismatch_len[0];
							int dist2 = align_info.rear_mismatch_len[1];
							int next_read = read_idx;
							int score = align_info.max_score;
							//int score = align_info.max_len;
							if (align_info.flip)
							{
								if (reads_info->RightOverlaps[next_read].count(-current_read) == 0)
								{
									reads_info->RightOverlaps[next_read][-current_read] = score;// (dist2);
								}
								next_read = -next_read;

							}
							else
							{
								if (reads_info->LeftOverlaps[next_read].count(current_read) == 0)
								{
									reads_info->LeftOverlaps[next_read][current_read] = score;//(dist2);
								}
							}

							if (reads_info->RightOverlaps[current_read].count(next_read) == 0)
							{
								reads_info->RightOverlaps[current_read][next_read] = score;// (dist1);
							}

							//
							//sorted_right_overlaps[score] = next_read;

						}
						if (result == 2)
						{
							int dist1 = align_info.frontal_mismatch_len[1];
							int dist2 = align_info.rear_mismatch_len[0];
							int next_read = read_idx;
							int score = align_info.max_score;
							//int score = align_info.max_len;
							if (align_info.flip)
							{
								if (reads_info->LeftOverlaps[next_read].count(-current_read) == 0)
								{
									reads_info->LeftOverlaps[next_read][-current_read] = score;//(dist1);
								}
								next_read = -next_read;
							}
							else
							{
								if (reads_info->RightOverlaps[next_read].count(current_read) == 0)
								{
									reads_info->RightOverlaps[next_read][current_read] = score;//(dist1);
								}
							}

							if (reads_info->LeftOverlaps[current_read].count(next_read) == 0)
							{
								reads_info->LeftOverlaps[current_read][next_read] = score;// (dist2);
							}

						}


						int next_read = read_idx;
						if (align_info.flip)
						{
							next_read = -next_read;
						}
						

					}

				}

			}


			map<int, int>::iterator temp_it;
			for (temp_it = reads_info->LeftOverlaps[current_read].begin(); temp_it != reads_info->LeftOverlaps[current_read].end(); ++temp_it)
			{
				sorted_left_overlaps[temp_it->second].push_back(temp_it->first);
			}

			for (temp_it = reads_info->RightOverlaps[current_read].begin(); temp_it != reads_info->RightOverlaps[current_read].end(); ++temp_it)
			{
				sorted_right_overlaps[temp_it->second].push_back(temp_it->first);
			}

			if (sorted_left_overlaps.size() > 0)
			{
				map<int, vector<int> >::reverse_iterator rit = sorted_left_overlaps.rbegin();
				reads_info->LeftOverlapsBest[current_read] = rit->second[0];



				for (int c = 0; c < rit->second.size(); ++c)
				{
					int next_read = rit->second[c];
					if (!reads_info->contained[abs(next_read)])
					{
						reads_info->LeftOverlapsBestWithTies[current_read].push_back(next_read);

					}

				}

				if (round == 2 && align_info.FAST == 0)// a little buggy
				{



					int non_contained_neighbors = 0;
					for (int c = 0; c < rit->second.size(); ++c)
					{
						if (!reads_info->contained[abs(rit->second[c])])
						{
							non_contained_neighbors++;
						}
					}

					for (int c = 0; c < rit->second.size(); ++c)
					{
						int next_read = rit->second[c];
						if (round == 2 && align_info.FAST == 0 && reads_info->contained[abs(next_read)])// a little buggy
						{
							if (non_contained_neighbors == 0)
							{
								reads_info->contained[abs(next_read)] = 0;
								//cout << next_read << ", ";
								reads_info->LeftOverlapsBestWithTies[current_read].push_back(next_read);
							}
						}
					}
				}


			}
			if (sorted_right_overlaps.size() > 0)
			{
				map<int, vector<int> >::reverse_iterator rit = sorted_right_overlaps.rbegin();
				reads_info->RightOverlapsBest[current_read] = rit->second[0];



				for (int c = 0; c < rit->second.size(); ++c)
				{
					int next_read = rit->second[c];
					if (!reads_info->contained[abs(next_read)])
					{
						reads_info->RightOverlapsBestWithTies[current_read].push_back(next_read);

					}

				}

				if (round == 2 && align_info.FAST == 0)// a little buggy
				{
					int non_contained_neighbors = 0;
					for (int c = 0; c < rit->second.size(); ++c)
					{
						if (!reads_info->contained[abs(rit->second[c])])
						{
							non_contained_neighbors++;
						}
					}

					for (int c = 0; c < rit->second.size(); ++c)
					{
						int next_read = rit->second[c];
						if (round == 2 && align_info.FAST == 0 && reads_info->contained[abs(next_read)])// a little buggy
						{
							if (non_contained_neighbors == 0 )
							{
								//cout << next_read << ", ";
								reads_info->contained[abs(next_read)] = 0;
								reads_info->RightOverlapsBestWithTies[current_read].push_back(next_read);
							}
						}
					}
				}

			}

		}
		time(&end_time);
		if (n_alignments > 0)
		{
			cout << "Avg alignment size: " << seq_len_sum / 2 / n_alignments << endl;
			if (align_info.SparseAlign)
			{
				cout << "Avg sparse alignment size: " << seq_len_sum_compressed / 2 / n_alignments << endl;

			}

		}
		cout << n_alignments << " alignments calculated." << endl;
		cout << "Round " << round << " takes " << difftime(end_time, beg_time) << " secs." << endl;
	}

	
	for (int i = 0; i < numReads; ++i)
	{
		if (reads_info->LeftOverlapsBestWithTies[i].size() == 0 && reads_info->RightOverlapsBestWithTies[i].size() == 0)
		{
			
			reads_info->contained[i] = 1;
		}
	}
	
	int n_contained = 0;
	for (int i = 1; i < numReads; ++i)
	{
		if (reads_info->contained[i])
		{
			n_contained++;
		}
	}
	cout << n_contained << " contained out of " << numReads << endl;

}

void OutputReadsInfo(reads_info *reads_info)
{
	ofstream o_chimera("ChimeraIdx.txt");
	ofstream o_contained("ContainedIdx.txt");
	for (size_t i = 0; i < reads_info->LenVec.size(); ++i)
	{
		if (reads_info->contained[i])
		{
			o_contained << i << endl;
		}
		if (reads_info->chimeric[i])
		{
			o_chimera << i << endl;
		}
	}
}



void LoadReadsInfo(reads_info *reads_info)
{
	ifstream in_chimera("ChimeraIdx.txt");
	ifstream in_contained("ContainedIdx.txt");

	reads_info->contained.clear();
	reads_info->contained.resize(reads_info->LenVec.size() + 10);
	reads_info->chimeric.clear();
	reads_info->chimeric.resize(reads_info->LenVec.size() + 10);

	size_t idx = 0;
	for (size_t i = 0; i < reads_info->LenVec.size(); ++i)
	{
		reads_info->chimeric[i] = 0;
		reads_info->contained[i] = 0;
	}
	while (in_chimera >> idx)
	{
		reads_info->chimeric[idx] = 1;
	}
	while (in_contained >> idx)
	{
		reads_info->contained[idx] = 1;
	}
}

#endif
