#ifndef __BASIC_DATA_STRUCTURE_H
#define __BASIC_DATA_STRUCTURE_H

#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <sstream>
#include <fstream>
#include "time.h"
#include "getRSS.h"

using namespace std;

void* myMalloc (size_t size) {
	void *p = malloc(size);
	if (NULL==p) {
		cerr << "malloc error: Couldn't allocate memory for size="<<size<<". Exiting with error.\n";
		cerr << "Current usage = "<<getCurrentRSS()<<"\n";
		cerr << "Peak usage = "<<getPeakRSS()<<"\n";
		exit(1);
	}
	return(p);
}

void* myCalloc (size_t num, size_t size) {
	void *p = calloc(num, size);
	if (NULL==p) {
		cerr << "calloc error: Couldn't allocate memory for "<<num<<" elements of size "<<size<<" each. Exiting with error.\n";
		cerr << "Current usage = "<<getCurrentRSS()<<"\n";
		cerr << "Peak usage = "<<getPeakRSS()<<"\n";
		exit(1);
	}
	return(p);
}


// These are the structures to save the k-mers.32 bases, 64 bases, 96 bases, 128 bases.
struct kmer_t
{
	uint64_t kmer:62,cov:2;
};

struct kmer_t2
{
	uint64_t kmer[2];
};

struct kmer_t3
{
	uint64_t kmer[3];
};

struct kmer_t4
{
	uint64_t kmer[4];
};




// edge structure, maximum length 25
struct edge_node
{
	uint64_t edge:50,edge_cov:4,len:6,masked:1,removed:1;
	struct edge_node *nxt_edge;
};
// edge structure, maximum length 64
struct edge_node2
{
	uint64_t edge[2];
	uint16_t edge_cov:8,len:8;
	struct edge_node2 *nxt_edge;
};

// information for a k-mer node 
struct kmer_info
{
	//flags, recording if the node is used in a search, whether the node have branches on its sides...
	//uint8_t used:1,split_left:1,split_right:1,removed:1,flip:1,marked:1,repeat:1,masked:1;
	//pointers to edge links 
	//struct edge_node *left;
	//struct edge_node *right;
	//uint16_t cov1;
	//for building scaffolds
	int32_t contig_no:28,flip:1,cov1:2;
	int32_t cod;
};

//split the above structure, round 1 of graph construction:
struct kmer_info_r1
{
	uint8_t cov1;	
};


//Bloom filter structure
struct BF_info
{
	uint8_t * BF_HT;
	uint64_t m;
	int d;
	bool Bloom;
};

// a complete k-mer node,
struct bucket
{
	struct kmer_t kmer_t;	//32 bp
	struct kmer_info kmer_info;
	bucket *nxt_bucket;
};

struct bucket2
{
	struct kmer_t2 kmer_t2;	//64 bp
	struct kmer_info kmer_info;
	bucket2 *nxt_bucket;
};

struct bucket3
{
	struct kmer_t3 kmer_t3;	
	struct kmer_info kmer_info;
	bucket3 *nxt_bucket;
};

struct bucket4
{
	struct kmer_t4 kmer_t4;	
	struct kmer_info kmer_info;
	bucket4 *nxt_bucket;
};


struct bucket0
{
	uint64_t *kmer_t;	//any kmer size
	struct kmer_info kmer_info;
	bucket0 *nxt_bucket;
};




// a structure recording the information of removed k-mers in the BFS bubble removal

struct bucket_rm
{
	struct kmer_t kmer_t;
	struct kmer_t merged_kmer;
	bool flip;
	bucket_rm *nxt_bucket;
};

struct bucket_rm2
{
	struct kmer_t2 kmer_t2;
	struct kmer_t2 merged_kmer;
	bool flip;
	bucket_rm2 *nxt_bucket;
};
struct bucket_rm3
{
	struct kmer_t3 kmer_t3;
	struct kmer_t3 merged_kmer;
	bool flip;
	bucket_rm3 *nxt_bucket;
};
struct bucket_rm4
{
	struct kmer_t4 kmer_t4;
	struct kmer_t4 merged_kmer;
	bool flip;
	bucket_rm4 *nxt_bucket;
};

struct bucket_rm0
{
	uint64_t *kmer_t;
	uint64_t *merged_kmer;
	bool flip;
	bucket_rm0 *nxt_bucket;
};

//bucket in round 1
struct bucket_r1
{
	struct kmer_t kmer_t;	
	bucket_r1 *nxt_bucket;
};

struct bucket2_r1
{
	struct kmer_t2 kmer_t2;	
	struct kmer_info_r1 kmer_info;
	bucket2_r1 *nxt_bucket;
};

struct bucket3_r1
{
	struct kmer_t3 kmer_t3;	
	struct kmer_info_r1 kmer_info;
	bucket3_r1 *nxt_bucket;
};

struct bucket4_r1
{
	struct kmer_t4 kmer_t4;	
	struct kmer_info_r1 kmer_info;
	bucket4_r1 *nxt_bucket;
};

struct bucket0_r1
{
	uint64_t * kmer_t;	
	struct kmer_info_r1 kmer_info;
	bucket0_r1 *nxt_bucket;
};
//hashtable structure

struct hashtable
{
	struct bucket **store_pos;
	size_t ht_sz;
	int round;
};

struct hashtable2
{
	struct bucket2 **store_pos;
	size_t ht_sz;
};


struct hashtable3
{
	struct bucket3 **store_pos;
	size_t ht_sz;
};

struct hashtable4
{
	struct bucket4 **store_pos;
	size_t ht_sz;
};

struct hashtable0
{
	struct bucket0 **store_pos;
	size_t ht_sz;
};

//read structure

struct read_t
{
	char tag[1000];
	bool error_nt[1000];
	char c_seq[10000];//char representation
	//char *c_seq;//char representation
	
	//uint64_t read_bits[10000];//bit representation
	uint64_t *read_bits;
	//char read[1000];//char representation
	int readLen;// read length
	int read_idx;
};



struct ref_read_t
{
	char tag[1000];
	uint64_t *read_bits;//bit representation
	size_t read_idx;
	int alloc_sz;
	int readLen;// read length
	int contig_no;

};


struct key_table
{
	bool in_use;
	list<uint64_t *> pblocks;
	int BytesPerKey;
	int KeysPerBlock;
	int current_block;
	int current_index;
};






struct read_index
{
	vector< map<uint64_t,bool> > repeat_maps;
	uint64_t repeat_cnt;
	int MaxMatch;
	vector<int> read_len_vt;
};


struct reads_table
{
	bool in_use;
	list<uint64_t *> pblocks;
	int BytesPerBlock;
	int current_block;
	int current_byte;
	int current_read;
	map<int,uint64_t *> read_ptr;
	vector<int> read_len_vt;
};

struct KmerInContig
{
	uint32_t contig_no:31,flip:1;
	int pos;
};

struct ContigInRead
{
	int ctg_no;
	int cov;
	int coord2;

};

struct LongReadContigIndex
{
	map<int, KmerInContig> LR2CTG;
	map<int, vector<int> > CTG2LR, CTG2LR_2;
	map<int, map<int, int> > CTG2Coords;
	map<int, map<int, int> > CTG2Offsets;
	//map<int, ContigInRead > CTG2LR, CTG2LR_2;


	map<int, ContigInRead > layout;
	int KmerCovTh;
	int nMatches;
	int BlockSize;
	bool FastMap;
	//	vector<int> matching_positions_LR;
	//vector<KmerInContig> matching_positions_ctg;
};

struct Coord_CTG_Cov
{
	int coord;//observed coords
	int contig_no;
	int cov;
	int coord2;//contig coord in a read
};

struct raw_overlap_str
{
	int score;
	int index;
	string raw_seq;
};

struct read_match_info
{
	uint32_t read_idx;
	uint32_t score;

};
struct reads_info
{

	vector< vector<Coord_CTG_Cov> > LongReadIndexVec, RefIndexVec, LongReadIndexVecCorr,ShortReadIndexVec;
	vector<string> tag_vec, input_files;
	
	vector< vector<int> > Contig2LongRead,Contig2Ref,Contig2ShortRead;
	vector<int> LenVec, RefLenVec;
	vector < vector<int>::iterator> LengthRank;
	vector<string> selected_long_reads_seq;
	vector<map<int, int > > LeftOverlaps, RightOverlaps;//record all left and right overlaps with scores
	vector<int> LeftOverlapsBest, RightOverlapsBest,degree_vec;//record non-tied single best overlap
	vector<map<int, int> > LeftBestOverlapsTemp, RightBestOverlapsTemp;//record overlaps with scores
	vector< vector<int> > LeftOverlapsBestWithTies, RightOverlapsBestWithTies;//record best overlaps with ties
	map<int, vector<int> > contig2contig_cov;
	map<int, vector<int> > Consensus_info;
	map<int, vector<read_match_info> > Consensus_info_full;
	int ConsensusBestN;
	map<int, map<int, int> > LeftNodes, RightNodes, LeftNodes_ShortReads, RightNodes_ShortReads, LeftNodes_LongReads, RightNodes_LongReads;
	map< vector<int>, vector<raw_overlap_str> > Selected_Overlaps;
	bool RemoveChimera,MSA,RecoverFalseNegatives, Clean,CloseGaps;
	int MinLen;
	int MinOverlap;
	int K_size;
	int KmerCovTh;
	double Redundancy;
	int TopOverlaps;
	int64_t TotalReads;
	int PathCovTh;
	int MaxReadLen,MaxCompressedReadLen;
	int max_reads;
	double HitRate, AdaptiveTh;
	int ContigCovTh,ChimeraTh;
	vector<bool> contained, used, used_vt_left, used_vt_right,both_stand_used,chimeric,gap_closing;
	int n_deleted_edges;
	int mode;
	bool inward;
	bool outward;
	int insert_size;
	size_t numReads;
	vector<vector<int> > aligned_nodes, gap_idx_vec;
	
    bool Debug;
    //bool Debug;
};


//contig graph

struct adjacent_contig_info
{
	int32_t dist_sum;
	int8_t cov;
	string bridge;

};




struct contigs_info
{
	int total_contigs;
	int K_size;
	vector<int> contig_sz_vt,kmer_cnt_vt,comp_vt;
	vector<int> contigs_hp_b,contigs_hp_e;
	vector<string> contigs_str,contig_tag;
	map<int, vector<int> > Length_ID;
	map<int, vector<int> > Cov_Length;
	string ContigsFileName;
	string ContigGraphName;
	//map<int, vector<int> > ctg_in_scf;
	//vector<vector<int> > scaffolds;
	//vector<vector<int> > gaps_in_scaffolds;
	vector < vector<double>::iterator> SortedAstats;
	vector<double> AstatsVec;
	vector < vector<int>::iterator> LengthRank;
	vector<int> cov_vt;
	int ContigCovTh;
	int KmerCovTh;
	bool CallConsensus;
	vector<bool> used;
	//vector<struct c_info> c_info_vt;
	//vector< map<int,struct scaffold_contig_info> > scaffold_adjacency_left,scaffold_adjacency_right;
	vector< map<int,struct adjacent_contig_info> > contig_adjacency_left,contig_adjacency_right;	
};




//path information in the BFS bubble removal

struct BFS_path_info
{
	int cov;
	int depth;
	int len;
	bool BothEndsUsed;
	struct bucket* last_bkt;
	struct edge_node* last_bkt_edge;
};





struct BFS_path_info2
{
	int cov;
	int depth;
	int len;
	struct bucket2* last_bkt;
	struct edge_node* last_bkt_edge;
};

struct BFS_path_info3
{
	int cov;
	int depth;
	int len;
	struct bucket3* last_bkt;
	struct edge_node* last_bkt_edge;
};

struct BFS_path_info4
{
	int cov;
	int depth;
	int len;
	struct bucket4* last_bkt;
	struct edge_node* last_bkt_edge;
};

struct BFS_path_info0
{
	int cov;
	int depth;
	int len;
	struct bucket0* last_bkt;
	struct edge_node* last_bkt_edge;
};


//outdated

struct BFS_path_info_ctg
{
	int cov;
	int depth;
	int len;
	int last_ctg;	
};


//a stack in the BFS
struct stacked_bucket
{
	bucket *bktptr;
	bool RightSearch;
	bool BothEndsUsed;
};
struct stacked_bucket2
{
	bucket2 *bktptr;
	bool RightSearch;
};
struct stacked_bucket3
{
	bucket3 *bktptr;
	bool RightSearch;
};
struct stacked_bucket4
{
	bucket4 *bktptr;
	bool RightSearch;
};

struct stacked_bucket0
{
	bucket0 *bktptr;
	bool RightSearch;
};

bool it_cmp(const vector<int>::iterator &a, const vector<int>::iterator &b)
{
	return (*a) > (*b);
}


bool get_a_fasta_read(ifstream & fasta_in, string &tag, string &str, string & n_tag)
{
	
	ifstream tmp_ifstream;
	string temp;
	if(!getline(fasta_in,temp))
	{return 0;}
	if(temp[temp.size()-1]=='\n'||temp[temp.size()-1]=='\r')
	{temp.resize(temp.size()-1);}

	str.clear();
	if(temp[0]=='>')
	{tag=temp;}
	else
	{
		tag=n_tag;
		str=temp;
	}

	
	while(getline(fasta_in,temp)&&temp.size()>0)
	{
		
		if(temp[temp.size()-1]=='\n'||temp[temp.size()-1]=='\r')
		{temp.resize(temp.size()-1);}

		if((temp.size()>0&&(temp[0]=='>'||temp[0]=='\n'||temp[0]=='\r')))
		{
			n_tag=temp;
			return 1;
		}
		else
		{
			str+=temp;
		}
		
	}
	return 1;
}


bool get_a_fastq_read(ifstream & fastq_in, string &tag, string &seq, string & quality)
{
	
	ifstream tmp_ifstream;
	string temp;
	if(!getline(fastq_in,temp))
	{return 0;}
	seq.clear();
	if(temp[0]=='@')
	{tag=temp;}
	else
	{
		return 0;
	}
	getline(fastq_in,seq);//seq
	if(seq[seq.size()-1]=='\n'||seq[seq.size()-1]=='\r')
	{seq.resize(seq.size()-1);}
	getline(fastq_in,temp);//'+'
	getline(fastq_in,quality);
	if(quality[quality.size()-1]=='\n'||quality[quality.size()-1]=='\r')
	{quality.resize(quality.size()-1);}

	return 1;
}


bool get_a_contig_path(ifstream & cp_in, string &tag, vector<struct Coord_CTG_Cov> &TempIndex, int &Len, int KmerCovTh, string & n_tag)
{
	string tmp1, tmp2, tmp3, str;
	int coord, contig, cov, coord2;
	TempIndex.clear();
	TempIndex.reserve(500);
	
	Coord_CTG_Cov Coord_CTG_Cov;
	

	string temp;
	if (!getline(cp_in, temp))
	{
		return 0;
	}
	if (temp.size() == 0)
	{
		return 0;
	}
	if (temp[temp.size() - 1] == '\n' || temp[temp.size() - 1] == '\r')
	{
		temp.resize(temp.size() - 1);
	}

	str.clear();
	if (temp[0] == '>' || temp[0] == '@')
	{
		tag = temp;
	}
	else
	{
		tag = n_tag;
		str = temp;
		stringstream ss(str);
		tmp1.clear();
		ss >> coord >> tmp1 >> contig >> tmp2 >> cov >> tmp3 >> coord2;
		if (tmp1.size() == 0 || tmp1[0] != ',')
		{
			Len = coord;
		}
		else
		{
			if (cov > KmerCovTh)
			{
				Coord_CTG_Cov.coord = coord;
				Coord_CTG_Cov.contig_no = contig;
				Coord_CTG_Cov.cov = cov;
				Coord_CTG_Cov.coord2 = coord2;

				TempIndex.push_back(Coord_CTG_Cov);
			}
		}
	}


	while (getline(cp_in, temp))
	{
		if (temp.size() == 0)
		{
			break;
		}
		if (temp[temp.size() - 1] == '\n' || temp[temp.size() - 1] == '\r')
		{
			temp.resize(temp.size() - 1);
		}

		if ((temp.size()>0 && (temp[0] == '@'||temp[0] == '>' || temp[0] == '\n' || temp[0] == '\r')))
		{

			n_tag = temp;
			return 1;
		}
		else
		{
			stringstream ss(temp);
			tmp1.clear();
			ss >> coord >> tmp1 >> contig >> tmp2 >> cov >> tmp3 >> coord2;
			if (tmp1.size() == 0 || tmp1[0] != ',')
			{
				Len = coord;
				continue;
			}
			if (cov > KmerCovTh)
			{
				Coord_CTG_Cov.coord = coord;
				Coord_CTG_Cov.contig_no = contig;
				Coord_CTG_Cov.cov = cov;
				Coord_CTG_Cov.coord2 = coord2;

				TempIndex.push_back(Coord_CTG_Cov);
			}
			
		}

	}
	return 1;
}


bool basic_quality_check(string &seq_s)
{
	bool good_read=1;
	int seq_sz=seq_s.size();

	for(int i=0;i<seq_sz;++i)
	{
		if(seq_s[i]!='A'&&seq_s[i]!='C'&&seq_s[i]!='G'&&seq_s[i]!='T'&&seq_s[i]!='N')
		{
			good_read=0;
			seq_s.clear();
			return good_read;
		}
	}

	int nN=seq_sz-1,isN=-1;
	for(int i=0;i<seq_sz;++i)
	{
						
		if(seq_s[i]=='-'||seq_s[i]=='N')
		{
			if(i<=seq_sz/2)
			{
				isN=i;
				continue;
			}
			else
			{
				nN=i-1;
				break;
			}
		}
	}
	int s=0;
	if((nN-isN)<=seq_sz/2)
	{
		good_read=0;
	}
					
	if(good_read==0)
	{
		seq_s.clear();
		return good_read;
	}

	if(isN>=0)
	{
		for(int i=isN+1;i<=nN;++i)
		{
			seq_s[s]=seq_s[i];
			s++;
		}
		seq_s[s]='\0';
		seq_s.resize(s);
	}
					
	return good_read;
}



//left shift and right shift of shift_sz bits of the whole bit array, arr_sz is the array length
static inline void L_shift_NB(uint64_t * bitsarr, int shift_sz,int arr_sz)
{

	uint64_t temp_arr[100];

	/*
	for (int i=0;i<arr_sz;++i)
	{
		temp_arr[i]=0;
	}
	memset(temp_arr,0,sizeof(uint64_t)*arr_sz);
	*/

	int jmp=shift_sz/64;
	int offset=shift_sz%64;

	for (int i=0;i<arr_sz;++i)
	{
		if(i+jmp+1<arr_sz)
		{

			uint64_t tt=0;
			if(offset==0)
			{
				tt=0;
			}
			else
			{
				tt=(bitsarr[i+jmp+1]>>(64-offset));
			}
			temp_arr[i]=((bitsarr[i+jmp]<<offset)|tt);
		}
		if(i+jmp+1==arr_sz)
		{temp_arr[i]=bitsarr[i+jmp]<<offset;}
		if(i+jmp+1>arr_sz)
		{temp_arr[i]=0;}

	}

	memcpy(bitsarr,temp_arr,sizeof(uint64_t)*arr_sz);
/*
	for (int i=0;i<arr_sz;++i)
	{
		bitsarr[i]=temp_arr[i];
	}
	*/
}

static inline void R_shift_NB(uint64_t * bitsarr, int shift_sz,int arr_sz)
{
	uint64_t temp_arr[100];
	/*
	for (int i=0;i<arr_sz;++i)
	{
		temp_arr[i]=0;
	}
	memset(temp_arr,0,sizeof(uint64_t)*arr_sz);
	*/
	int jmp=shift_sz/64;
	int offset=shift_sz%64;

	for (int i=arr_sz-1;i>=0;--i)
	{
		if(i-jmp>0)
		{
			if(offset>0)
			{temp_arr[i]=(bitsarr[i-jmp]>>offset)|(bitsarr[i-jmp-1]<<(64-offset));}
			else
			{temp_arr[i]=bitsarr[i-jmp];}
		}
		if (i-jmp==0)
		{
			if(offset>0)
			{temp_arr[i]=(bitsarr[i-jmp]>>offset);}
			else
			{temp_arr[i]=bitsarr[i-jmp];}
		}
		if (i-jmp<0)
		{temp_arr[i]=0;}

	}
	memcpy(bitsarr,temp_arr,sizeof(uint64_t)*arr_sz);
	/*
	for (int i=0;i<arr_sz;++i)
	{
		bitsarr[i]=temp_arr[i];
	}
	*/

}

// get reverse complement of a k-mer.
static inline uint64_t get_rev_comp_seq(uint64_t seq, int seq_size)
{
	seq =~seq;

	seq = ((seq & 0x3333333333333333 )<< 2) | ((seq & 0xCCCCCCCCCCCCCCCC )>> 2);
	seq = ((seq & 0x0F0F0F0F0F0F0F0F )<< 4) | ((seq & 0xF0F0F0F0F0F0F0F0 )>> 4);
	seq = ((seq & 0x00FF00FF00FF00FF )<< 8) | ((seq & 0xFF00FF00FF00FF00 )>> 8);
	seq = ((seq & 0x0000FFFF0000FFFF )<<16) | ((seq & 0xFFFF0000FFFF0000 )>>16);
	seq = ((seq & 0x00000000FFFFFFFF )<<32) | ((seq & 0xFFFFFFFF00000000 )>>32);

	return seq >> (64 - (seq_size*2));
}

static inline uint64_t* get_rev_comp_seq_arr(uint64_t *seq_arr, int seq_size,int arr_sz)
{


	int tot_bits=arr_sz*64;
	for(int i=0;i<arr_sz;++i)
	{
		seq_arr[i]=~seq_arr[i];
		seq_arr[i] = ((seq_arr[i] & 0x3333333333333333 )<< 2) | ((seq_arr[i] & 0xCCCCCCCCCCCCCCCC )>> 2);
		seq_arr[i] = ((seq_arr[i] & 0x0F0F0F0F0F0F0F0F )<< 4) | ((seq_arr[i] & 0xF0F0F0F0F0F0F0F0 )>> 4);
		seq_arr[i] = ((seq_arr[i] & 0x00FF00FF00FF00FF )<< 8) | ((seq_arr[i] & 0xFF00FF00FF00FF00 )>> 8);
		seq_arr[i] = ((seq_arr[i] & 0x0000FFFF0000FFFF )<<16) | ((seq_arr[i] & 0xFFFF0000FFFF0000 )>>16);
		seq_arr[i] = ((seq_arr[i] & 0x00000000FFFFFFFF )<<32) | ((seq_arr[i] & 0xFFFFFFFF00000000 )>>32);
	}

	int j=0,k=arr_sz-1;
	for (;j<k;++j,--k)
	{
		uint64_t temp;
		temp=seq_arr[j];
		seq_arr[j]=seq_arr[k];
		seq_arr[k]=temp;
	}

	R_shift_NB(seq_arr,tot_bits-(seq_size*2),arr_sz);
	return seq_arr;
	//return seq >> (64 - (seq_size*2));
}

// get sub bit array from a bit array.
 inline void get_sub_arr(uint64_t * bitsarr_in,int bitsarr_len,int begin_pos,int sub_sz,uint64_t * bitsarr_out)
{
	if(bitsarr_len<sub_sz)
	{cout<<"Error! Input kmer too short."<<bitsarr_len <<" "<<sub_sz<<endl;return;}
	int arr_sz_in=bitsarr_len/32+1;
	int rem=bitsarr_len%32;
	if(rem==0)
	{arr_sz_in--;}

	int arr_sz_out=sub_sz/32+1;
	if(sub_sz%32==0)
	{arr_sz_out--;}

	uint64_t temp_arr[10];
	memset(temp_arr,0,sizeof(temp_arr));

	memset(bitsarr_out,0,sizeof(uint64_t)*arr_sz_out);

	int rem2=(32-rem+begin_pos)%32;
	int block_beg=(32-rem+begin_pos)/32;
	if(rem==0)
	{block_beg--;}

	int rem3=(32-rem+begin_pos+sub_sz)%32;
	int block_end=(32-rem+begin_pos+sub_sz)/32;
	if(rem3!=0)
	{rem3=32-rem3;}
	else
	{
		block_end--;
	}
	if(rem==0)
	{block_end--;}

	int orig_sz=(block_end-block_beg+1);
	memcpy(temp_arr,&bitsarr_in[block_beg],orig_sz*sizeof(uint64_t));
	L_shift_NB(temp_arr,rem2*2,orig_sz);
	R_shift_NB(temp_arr,(rem2+rem3)%32*2,arr_sz_out);
	memcpy(bitsarr_out,temp_arr,arr_sz_out*sizeof(uint64_t));


}

 uint64_t* str2bitsarr(const char * c_str,int len, uint64_t* b_str,int arr_sz )
{

	for (int  k=0;k<arr_sz;++k)
	{
		b_str[k]=0;
	}
	int arr_sz_needed=len/32+1;
	int rem=len%32;
	if(rem==0)
	{arr_sz_needed--;}

	int beg_arr_idx=arr_sz-arr_sz_needed;
	if(rem==0&&arr_sz_needed>0)
	{rem=32;}
	for (int k=0;k<len;k++)
	{
		if(rem==0)
		{beg_arr_idx++;rem=32;}


		switch(c_str[k])
		{
		case ('A'):case ('a'):case ('0'):
			b_str[beg_arr_idx]<<=2;//L_shift_NB(b_str,2,arr_sz);
			rem--;
			//b_str<<=2;
			break;


		case ('C'):case ('c'):case ('1'):
			b_str[beg_arr_idx]<<=2;//L_shift_NB(b_str,2,arr_sz);
			++b_str[beg_arr_idx];//++b_str[arr_sz-1];
			rem--;
//			++(b_str<<=2);
			break;


		case 'G':case 'g':case '2':
			b_str[beg_arr_idx]<<=1;//L_shift_NB(b_str,1,arr_sz);
			++b_str[beg_arr_idx];//++b_str[arr_sz-1];
			b_str[beg_arr_idx]<<=1;//L_shift_NB(b_str,1,arr_sz);
			rem--;//(++(b_str<<=1))<<=1;
			break;

		case 'T':case 't':case '3':
			b_str[beg_arr_idx]<<=1;//L_shift_NB(b_str,1,arr_sz);
			++b_str[beg_arr_idx];
			b_str[beg_arr_idx]<<=1;//L_shift_NB(b_str,1,arr_sz);
			++b_str[beg_arr_idx];
			rem--;
			//++((++(b_str<<=1))<<=1);
			break;
		default:
			return b_str;
		}

	//	cout<<b_str<<endl;
	}
	return b_str;
}

 //hash functions
uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed )
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);

	while(data != end)
	{
		uint64_t k = *data++;

		k *= m;
		k ^= k >> r;
		k *= m;

		h ^= k;
		h *= m;
	}

	const unsigned char * data2 = (const unsigned char*)data;

	switch(len & 7)
	{
	case 7: h ^= uint64_t(data2[6]) << 48;
	case 6: h ^= uint64_t(data2[5]) << 40;
	case 5: h ^= uint64_t(data2[4]) << 32;
	case 4: h ^= uint64_t(data2[3]) << 24;
	case 3: h ^= uint64_t(data2[2]) << 16;
	case 2: h ^= uint64_t(data2[1]) << 8;
	case 1: h ^= uint64_t(data2[0]);
	        h *= m;
	};

	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
}

uint64_t MurmurHash64B ( const void * key, int len, unsigned int seed )
{
	const unsigned int m = 0x5bd1e995;
	const int r = 24;

	unsigned int h1 = seed ^ len;
	unsigned int h2 = 0;

	const unsigned int * data = (const unsigned int *)key;

	while(len >= 8)
	{
		unsigned int k1 = *data++;
		k1 *= m; k1 ^= k1 >> r; k1 *= m;
		h1 *= m; h1 ^= k1;
		len -= 4;

		unsigned int k2 = *data++;
		k2 *= m; k2 ^= k2 >> r; k2 *= m;
		h2 *= m; h2 ^= k2;
		len -= 4;
	}

	if(len >= 4)
	{
		unsigned int k1 = *data++;
		k1 *= m; k1 ^= k1 >> r; k1 *= m;
		h1 *= m; h1 ^= k1;
		len -= 4;
	}

	switch(len)
	{
	case 3: h2 ^= ((unsigned char*)data)[2] << 16;
	case 2: h2 ^= ((unsigned char*)data)[1] << 8;
	case 1: h2 ^= ((unsigned char*)data)[0];
			h2 *= m;
	};

	h1 ^= h2 >> 18; h1 *= m;
	h2 ^= h1 >> 22; h2 *= m;
	h1 ^= h2 >> 17; h1 *= m;
	h2 ^= h1 >> 19; h2 *= m;

	uint64_t h = h1;

	h = (h << 32) | h2;

	return h;
}

 //convert a string of nucleotide bases into bit array.
void Init_Read(string &seq,struct read_t & read)
{
	read.readLen=strlen(seq.c_str());
	int Read_arr_sz=read.readLen/32+1;
	int rem=read.readLen%32;
	if(rem==0)
	{Read_arr_sz--;}

	str2bitsarr(seq.c_str(),(int)seq.size(),read.read_bits,Read_arr_sz);

}



void Init_Ref_Read(string &seq, struct ref_read_t & read)
{
	read.readLen = strlen(seq.c_str());
	int Read_arr_sz = read.readLen / 32 + 1;
	int rem = read.readLen % 32;
	if (rem == 0)
	{
		Read_arr_sz--;
	}

	str2bitsarr(seq.c_str(), (int)seq.size(), read.read_bits, Read_arr_sz);

}

//static 


//get the complement of a string of nucleotide bases
void complement_str(string & str)
{
	for (size_t i=0;i!=str.size();++i)
	{
		switch (str[i])
		{
		case 'A':case 'a':
			str[i]='T';
			break;
		case 'C':case 'c':
			str[i]='G';
			break;
		case 'G':case 'g':
			str[i]='C';
			break;
		case 'T':case 't':
			str[i]='A';
			break;
		case 'N':case 'n':
			str[i]='N';
			break;
		case '-':
			str[i]='-';
			break;
		default:
			cout<<"error: complement_str"<<str[i]<<endl;
			return;


		}

	}
}


void reverse_complement_str(string & str)
{
	if (str.size() == 0)
	{
		return;
	}

	reverse(str.begin(), str.end());
	for (size_t i = 0; i != str.size(); ++i)
	{
		switch (str[i])
		{
		case 'A':case 'a':
			str[i] = 'T';
			break;
		case 'C':case 'c':
			str[i] = 'G';
			break;
		case 'G':case 'g':
			str[i] = 'C';
			break;
		case 'T':case 't':
			str[i] = 'A';
			break;
		case 'N':case 'n':
			str[i] = 'N';
			break;
		case '-':
			str[i] = '-';
			break;
		default:
			cout << "error: complement_str" << str[i] << endl;
			return;


		}

	}
}
//convert a bit array into a string
char * bitsarr2str(uint64_t* b_seq, int len,char * c_str,int arr_sz)
{

	int tot_bits=arr_sz*64;
	//char *c_str;
	//c_str=(char*) malloc(sizeof(char)*(len+1));
	//#pragma omp parallel for
	for (int i=0;i<len;++i)
	{
		uint64_t temp,temp2[100];/////////////////////////
		for (int k=0;k<arr_sz;++k)
		{
			temp2[k]=b_seq[k];
		}
		L_shift_NB(temp2,tot_bits-(len-i)*2,arr_sz);
		R_shift_NB(temp2,tot_bits-2,arr_sz);
		//uint64_t temp=(b_seq<<(64-(len-i)*2))>>62;
		temp=temp2[arr_sz-1];
		switch(temp)
		{
		case 0:
			c_str[i]='A';
			break;
		case 1:
			c_str[i]='C';
			break;
		case 2:
			c_str[i]='G';
			break;
		case 3:
			c_str[i]='T';
			break;

		}
	}
	c_str[len]='\0';
	return c_str;
}



void string_shrinkage(string &s)
{
	int s_sz=s.size();
	string t=s;
	s.clear();
	s.reserve(s_sz);
	s.push_back(t[0]);
	for(int i=0;i<s_sz; ++i)
	{
		if(t[i]!=s[s.size()-1])
		{
			s.push_back(t[i]);
		}
	}

}



bool it_cmp_double(const vector<double>::iterator &a, const vector<double>::iterator &b)
{
	return (*a) > (*b);
}



#endif
