#include<fstream>
#include<iostream>
#include<vector>
#include<string>
#include<stdio.h>
#include<map>
#include<stack>
#include<stdlib.h>
#include<string.h>
#include<set>
#include<time.h>
#include<unistd.h>
#include<getopt.h>
#define search_range 4
unsigned int get2bit[] =	//to get each code of character 
{
    (1 << 30) + (1 << 31),
    (1 << 28) + (1 << 29),
    (1 << 26) + (1 << 27),
    (1 << 24) + (1 << 25),
    (1 << 22) + (1 << 23),
    (1 << 20) + (1 << 21),
    (1 << 18) + (1 << 19),
    (1 << 16) + (1 << 17),
	(1 << 14) + (1 << 15),
	(1 << 12) + (1 << 13),
	(1 << 10) + (1 << 11), 
	(1 << 8) + (1 << 9), 
	(1 << 6) + (1 << 7),
	(1 << 4) + (1 << 5),
   	(1 << 2) + (1 << 3), 
	(1 << 0) + (1 << 1),
};

unsigned char get1bit[] = //to get each bit
{
	1 << 7,
	1 << 6,
	1 << 5,
	1 << 4,
	1 << 3,
	1 << 2,
	1 << 1,
	1 << 0
};

struct globalArgs_t
{
	int k;
	int n;
	int lnth0;
	int lnth1;
	int d1;
	int d2;
	int kmercov;
	int readcov;
	int readcov0;
	int readlen;
	int err;
	int paircov;
	int Min;
	int Max;
	char * output;
	char * contig_info;
	char * MUP_info;
	char * single;
	char * left;
	char * right;
} globalArgs;



static const char * optString = "k:n:m:L:d:D:c:C:R:e:o:1:2:s:l:r:p:S:M:";




using namespace std;
typedef struct	//structure to save reads
{
	short size;
	unsigned int * p;
} reads;
typedef reads* ReadTable;	//pointer to structure that save reads

struct from_read; //to indicate which read it is from
typedef struct from_read from_read;

typedef struct					//structure of kmer
{
	unsigned int mkey, nkey;			//hash keys
  //we set isdel 1 which means it is deleted
	bool isdel;			
	unsigned char n;				//record it's neighboor
	from_read * head;		//to indicate which read it is from
	int MUP_ind;			//to indicate which MUP it belong to 
} edge;

struct from_read
{
	unsigned int i;
	from_read* next;
};

typedef map<unsigned int, edge> tree;	//nkey represents the key of map
typedef vector<tree> KmerTable;	//mkey represents the index of vector





unsigned int char_to_key(int c)    //if it contains 'N', we set it as a random value
{
	switch (c)
	{
		case 'A':
			return 0;
		case 'G':
			return 1;
		case 'C':
			return 2;
		case 'T':
			return 3;
		default:
			return rand() % 4;
	}
}


//read reads from file to ReadTable, return 0 if there is no input data
//pairkey save the position for paired reads
//single point to file name if it's single-end
//left and right point to files if it's pairend
//return a pointer which point the stucture save all coded reads
ReadTable Input_Reads(char * single, char * left, char * right, unsigned int &pairkey)	
{
	vector<string> V;
	string s1;
	string s2;
	string s3;
	string s4;
	if (single != 0)
	{
		ifstream in(single);
		pairkey = 0;
		while(getline(in, s1) && getline(in, s2) && getline(in, s3) && getline(in, s4))
			V.push_back(s2);
	}
	else if(left != 0 && right != 0)
	{
		ifstream in1(left);
		ifstream in2(right);
		while(getline(in1, s1) && getline(in1, s2) && getline(in1, s3) && getline(in1, s4))
			V.push_back(s2);
		if ((pairkey = V.size()) == 0)
			return 0;
		while(getline(in2, s1) && getline(in2, s2) && getline(in2, s3) && getline(in2, s4))
			V.push_back(s2);
		if ((pairkey << 1) != V.size())
			return 0;
	}
	ReadTable p = (ReadTable)malloc(sizeof(reads) * (V.size() + 1));
	if (p == 0)
	{
		printf("short of memory!\n");
		return 0;
	}
	int state = 0;
	reads * temp = 0;
	string s;
	for (size_t i = 0; i != V.size(); ++i)
	{
		temp = p + i;
		s = V[i];
		temp->size = (short)s.size();
		int size = temp->size / 16 + 1;
		temp->p = (unsigned int*)malloc(sizeof(unsigned int) * size);
		unsigned int * q = temp->p;
		for (int k = 0; k != size; ++k)
			*q++ = 0;
		q = temp->p;
		int k = 0;
		int l = -1;
		for (int j = 0; j != temp->size; ++j)
		{
			k = j % 16;
			if (k == 0)
				++l;
			unsigned int result = 0;
			result = char_to_key(s.at(j)); 
			*(q + l) += (result << 2 * (15 - k));
		}
	}
	(p + V.size())->size = -1;
	(p + V.size())->p = 0;
	return p;
}



inline int kmer2hash(reads * r, unsigned int &mkey, unsigned int &nkey, int startpos, int k, int n) 
{
	if (startpos + k > r->size)
		return 0;
	mkey = 0;
	nkey = 0;
	int result = startpos / 16;
	int remain = startpos % 16;
	int m = k - n;
	unsigned int temp = 0;
	unsigned int * p = r->p + result;
	if (remain + m <= 16)
	{
		mkey = (*p << 2 * remain) >> (2 * (16 - m));
		if (remain + m + n <= 16)
			nkey = (*p << 2 * (remain + m)) >> 2 * (16 - n);
		else
		{
			unsigned int n_remain = 16 - remain - m;
			if (n_remain == 0)
			{
				++p;
				nkey = *p >> 2 * (16 - n);
			}
			else
			{
				nkey = (*p << 2 * (remain + m)) >> 2 * (16 - n_remain);
				++p;
				nkey = (nkey << 2 * (n - n_remain)) + (*p >> 2 * (16 - n + n_remain));
			}
		}
	}
	else
	{
		unsigned int m_remain = 16 - remain;
		mkey = (*p << 2 * remain) >> 2 * (16 - m_remain);
		++p;
		mkey = (mkey << 2 * (m - m_remain)) + (*p >> 2 * (16 - m + m_remain));
		if (m - m_remain + n <= 16)
			nkey = (*p << 2 * (m - m_remain)) >> (2 * (16 - n));
		else
		{
			unsigned int n_remain = 16 - m + m_remain;
			nkey = (*p << 2 * (m - m_remain)) >> 2 * (16 - n_remain);
			++p;
			nkey = (nkey << 2 * (n - n_remain)) + (*p >> 2 * (16 - n + n_remain));
		}
	}
	return 1;
}

//translate code to its respondence character
inline char key_to_char(int key)	
{
		switch (key)
		{
			case 0:
				return 'A';
			case 1:
				return 'G';
			case 2:
				return 'C';
			case 3:
				return 'T';
			default:
				return 'N';
		}
}


//to get the kmer which hash value is mkey and nkey,
//can't verify if it is a valid hash value!
//return a string which key value are mkey and nkey
string inv_hash_complete(unsigned int mkey, unsigned int nkey, int k, int n)
{
	string s;
	unsigned int t = 0;
	int m = k - n;
	for (int i = 0; i != m; ++i)
	{
		t = (mkey >> 2 * (m - 1 - i)) & 3;
		s += key_to_char(t);
	}
	for (int i = 0; i != n; i++)
	{
		t = (nkey >> 2 * (n - 1 - i)) & 3;
		s += key_to_char((int)t);
	}
	return s;
}


//compute the hash key for kmer
inline void hash_kmer(unsigned int c, unsigned int &mkey, unsigned int &nkey, int k, int n)	
{
	int m = k - n;
	unsigned int t = 0;
	t = nkey >> 2 * (n - 1);
	t &= 3;
	mkey <<= 2;
	if (m < 16)
		mkey = (mkey + t) & ((1 << 2 * m) - 1);
	nkey <<= 2;
	if (n < 16)
	nkey &= ((1 << 2 * n) - 1);
	t = c;
	nkey += t;
	return;
}

inline void init_edge(edge &e)	//initalize the kmer
{
	e.mkey = 0;
	e.nkey = 0;
	e.isdel = 0;
	e.MUP_ind = -1;
	e.head = 0;
	e.n = 0;
}

unsigned int get_edge_cov(KmerTable &T, edge * e, int k, int n)
{
	unsigned int result = 0;
	int m = k - n;
	if (e == 0)
		return result;
	from_read * p = e->head;
	if (p == 0)
	{
		unsigned int mkey = e->mkey;
		unsigned int nkey = e->nkey;
		if (m == 16)
			mkey = ~mkey;
		else
			mkey = (1 << 2 * m) - mkey - 1;
		if (n == 16)
			nkey = ~nkey;
		else
			nkey = (1 << 2 * n) - nkey - 1;
		if (T[mkey].count(nkey) > 0)
		{
			p = T[mkey][nkey].head;
			while (p != 0)
			{
				++result;
				p = p->next;
			}
		}
	}
	else
		while (p != 0)
		{
			++result;
			p = p->next;
		}
	return result;
}

from_read * Creat_from_read_Node(unsigned int item)
{
	from_read * p = (from_read*)malloc(sizeof(from_read));
	p->i = item;
	p->next = NULL;
	return p;
}

from_read * read_index_insert(from_read * head, unsigned int index)
{
	from_read * p = NULL;
	if (head == NULL)
	{
		if (!(p = Creat_from_read_Node(index)))
		{
			printf("out of memory\n");
			return 0;
		}
		return p;
	}
	else
	{
		p = head;
		if (p->i != index)
		{
			if (!(p = Creat_from_read_Node(index)))
			{
				printf("out of memory\n");
				return 0;
			}
			p->next = head;
		}
		return p;
	}
}


//insert kmer with hash value(mkey, nkey) 
//index represents which read it is from
void insert_kmer(unsigned int mkey, unsigned int nkey, unsigned int index, KmerTable &T)	
{ 
	//initalize the kmer
	edge e;
	init_edge(e);
	e.mkey = mkey;
	e.nkey = nkey;
	pair<tree::iterator, bool> result;
	result = T[mkey].insert(make_pair(nkey, e));
	from_read * p = result.first->second.head;
	p = read_index_insert(p, index);
	if (p != 0)
		result.first->second.head = p;
	return;
}


//input reads to generate kmers
void read_to_kmer(reads * s, KmerTable &T, int k, int n, unsigned index)	
{
	if (s->size < k)
		return;
	unsigned int mkey = 0;
	unsigned int nkey = 0;
	int l = (k - 1) / 16;
	int t = 0;
	unsigned int temp = 0;
	kmer2hash(s, mkey, nkey, 0, k, n);
	insert_kmer(mkey, nkey, index, T);
	for (unsigned int i = k; i < s->size; ++i)
	{
		t = i % 16;
		if (t == 0)
			l++;
		temp = *(s->p + l);
		temp &= get2bit[t];
		temp >>= 2 * (15 - t);
		hash_kmer(temp, mkey, nkey, k, n);
		insert_kmer(mkey, nkey, index, T);
	}
	return;
}


//insert kmer with hash value(mkey, nkey) 
//index represents which read it is from
void insert_kmer2(unsigned int mkey, unsigned int nkey, unsigned int index, KmerTable &T)	
{ 
	//initalize the kmer
	edge e;
	init_edge(e);
	e.mkey = mkey;
	e.nkey = nkey;
	pair<tree::iterator, bool> result;
	result = T[mkey].insert(make_pair(nkey, e));
	return;
}

void read_to_kmer2(reads * s, KmerTable &T, int k, int n, unsigned index)	
//input reads to generate kmers
{
	if (s->size < k)
		return;
	unsigned int mkey = 0;
	unsigned int nkey = 0;
	int l = (k - 1) / 16;
	int t = 0;
	unsigned int temp = 0;
	kmer2hash(s, mkey, nkey, 0, k, n);
	insert_kmer2(mkey, nkey, index, T);
	for (unsigned int i = k; i < s->size; ++i)
	{
		t = i % 16;
		if (t == 0)
			l++;
		temp = *(s->p + l);
		temp &= get2bit[t];
		temp >>= 2 * (15 - t);
		hash_kmer(temp, mkey, nkey, k, n);
		insert_kmer2(mkey, nkey, index, T);
	}
	return;
}


void complement(reads * p)
{
	short size = p->size;
	unsigned int * q = p->p;
	int k = size / 16;
	int mod = size % 16;
	while (k > 0)
	{
		*q = ~(*q);
		++q;
		--k;
	}
	if (mod != 0)
	{
		unsigned int result = *q;
		result = (result >> 2 * (16 - mod));
		result = (1 << 2 * mod) - 1 - result;
		result = (result << 2 * (16 - mod));
		*q = result;
	}
	return;
}


//NOTICE:Build_Kmer don't use complement infomation of reads!!!!!!!!!!!!!
void Build_Kmer(KmerTable &T, ReadTable V, int k, int n)	
{
	reads * p = V;
#ifndef RELEASE
	unsigned int count = 0;
	clock_t start = clock();
	clock_t begin = start;
	clock_t end;
	double time = 0;
#endif
	if (p == 0)
	{
		printf("No input file!\n");
		return;
	}
	while (p->size != -1)
	{
		read_to_kmer(p, T, k, n, count);
#ifndef RELEASE
		if (count % 1000000 == 0)
		{
			end = clock();
			time = (double)(end - start) / CLOCKS_PER_SEC;
			printf("%lf\n", (double)count);
			printf("time need to process 10,000,000 reads is %lf\n", time);
			start = end;
		}
#endif
		++p;
		++count;
	}
	
#ifndef RELEASE
	time = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("%lf\n", (double)count);
	printf("Total time needed to Build Kmer graph is %lf\n", time);
#endif
	p = V;
	while (p->size != -1)
	{
		complement(p);
		read_to_kmer2(p, T, k, n, count);
		complement(p);
		++p;
	}
	return;
}


//make isdel in Kmer to 1 if its coverage below threshold
void Del_Low_Kmer(KmerTable &T, unsigned int kmer_threshold, int k, int n)	
{
	KmerTable::iterator itr;
	tree::iterator t;
	KmerTable::iterator end = T.end();
	tree::iterator t_end;
	for (itr = T.begin(); itr != end; ++itr)
	{
		if (itr->size() > 0)
		{
			t_end = itr->end();
			for (t = itr->begin(); t != t_end; ++t)
				if (get_edge_cov(T, &t->second, k, n) < kmer_threshold)
					t->second.isdel = 1;
		}
	}
	return;
}


//return pointer to kmer if it finds; return NULL if it can't find
//find kmer based on hash value
edge * find_kmer(unsigned int mkey, unsigned int nkey, KmerTable &T){
	edge * p = NULL;
	if (T[mkey].count(nkey) > 0)
	{
		p = &T[mkey][nkey];
		return p;
	}
	else
		return 0;
}


//find right neighbourhood of kmer.
//	We don't use the kmer as neighbour if it is deleted!
void find_rnei(edge &e, int k, int n, KmerTable &T)	
{
	int m = k - n;
	unsigned int mkey = 0;
	unsigned int nkey = 0;
	edge * p = NULL;
	int temp = (e.nkey >> 2 * (n - 1)) & 3;
	mkey = (e.mkey << 2) + temp;
	if (m != 16)
		mkey &= (1 << 2 * m) - 1;
	nkey = (e.nkey << 2);
	if (n != 16)
	nkey &= (1 << 2 * n) - 1;
	for (int i = 0; i != 4; i++)
	{
		if (p = find_kmer(mkey, nkey + i, T))
			if (p->isdel == 0)
				e.n += get1bit[i + 4];
	}
	return;
}


//find left neighbourhood of kmer
void find_lnei(edge &e, int k, int n, KmerTable &T)
{
	int m = k - n;
	unsigned int mkey = 0;
	unsigned int nkey = 0;
	edge * p = NULL;
	int temp = e.mkey & 3;
	nkey = (e.nkey >> 2) + (temp << 2 * (n - 1));
	if (n != 16)
	nkey &= (1 << 2 * n) - 1;
	mkey = (e.mkey >> 2);
	if (m != 16)
	   mkey &= (1 << 2 * m) - 1;
	for (int i = 0; i != 4; i++)
	{
		if (p = find_kmer(mkey + (i << 2 * (m - 1)), nkey, T))
			if (p->isdel == 0)
				e.n += get1bit[i];
	}
	return;
}


//to build Debrujin Graph for all kmers
//NOTICE: we only consider kmers that are not deleted!
void Build_DeBrujin(KmerTable &T, int threshold, int k, int n)	
{
	tree::iterator t;
	tree::iterator t_end;
#ifndef RELEASE
	int num = 0;
#endif
	for (size_t i = 0; i != T.size(); ++i)
	{
		if (T[i].size() > 0)
		{
#ifndef RELEASE
			num++;
#endif
			t_end = T[i].end();
			for (t = T[i].begin(); t != t_end; ++t)
				if (!t->second.isdel)
				{
					find_lnei(t->second, k, n, T);
					find_rnei(t->second, k, n, T);
				}
		}
	}
	return;
}

typedef set<int> read_contain;	//it contains reads
//int in map<int, read_contain> represent distance between MUPs which support by reads in read_contain
typedef map<int, read_contain> gap_record;
//int in map<int, gap_record> represent the index of neighbour MUP
typedef map<int, gap_record> rnei_record;

struct MUP;
struct MUP_rnei
{
	struct MUP * p;
	short dist;
	struct MUP_rnei * next;
};
struct MUP_lnei
{
	struct MUP * p;
	struct MUP_lnei * next;
};
typedef struct MUP_rnei MUP_rnei;
typedef struct MUP_lnei MUP_lnei;

struct MUP
{
	edge * start;	//starting kmer for MUP
	edge * end;		//ending kmer for MUP
	int length;		//length of MUP which is equal to the number of kmers
	MUP_lnei * lnei;
	MUP_rnei * rnei;	//structure of MUP's right neighbour
	int index;	//index of long_MUP to which MUP belongs
};
typedef struct MUP MUP;
typedef vector<MUP> MUP_Graph;	//all MUPs are saved to MUP-graph
//NOTICE:MUP_graph[0] is a special MUP represents complex region

struct LongMUP;
struct LongMUP_lnei
{
	struct LongMUP * p;
	struct LongMUP_lnei * next;
};
struct LongMUP_rnei
{
	struct LongMUP * p;
	struct LongMUP_rnei * next;
	int dist;
};
typedef struct LongMUP_lnei LongMUP_lnei;
typedef struct LongMUP_rnei LongMUP_rnei;
struct LongMUP
{
	MUP * start;
	MUP * end;
	int length;
	LongMUP_lnei * lnei;
	LongMUP_rnei * rnei;
	bool isused;
};

typedef struct LongMUP LongMUP;
typedef vector<LongMUP> LongMUP_Graph;
typedef vector<LongMUP*> ContigPath_record;
typedef vector<ContigPath_record> Contig_Pathway;

void init_MUP(MUP &mup, edge * start = 0, edge * end = 0, int length = -1)	//initalize MUP
	//NOTICE:in_contig is assigned to -3 for initialization
{
	mup.start = start;
	mup.end = end;
	mup.length = length;
	mup.index = -1;
	mup.lnei = 0;
	mup.rnei = 0;
	return;
}


inline int get_kmer_nei_num(edge * p, int isleft)
{
	int count = 0;
	unsigned char c = p->n;
	if (isleft)
	{
		for (int i = 0; i != 4; ++i)
//8 is the limit of c
//one k_mer has most 8 neighboor
		{
			if ((c & get1bit[i]) != 0)
				++count;
		}
		return count;
	}
	else
	{
		for (int i = 4; i != 8; ++i)
		{
			if ((c & get1bit[i]) != 0)
				++count;
		}
		return count;
	}
}

edge * get_kmer_nei(edge * p, KmerTable &T, int index)
//index range is 0, 1, 2, 3, 4, 5, 6, 7
{
	int m = globalArgs.k - globalArgs.n;
	int n = globalArgs.n;
	int has_nei = p->n & get1bit[index];
	if (!has_nei)
		return 0;
	unsigned int key = index % 4;
	if (index < 4)
	{
		unsigned int mkey = 0;
		unsigned int nkey = 0;
		unsigned int temp = p->mkey & 3;
		nkey = (p->nkey >> 2) + (temp << 2 * (n - 1));
		if (n != 16)
		nkey &= (1 << 2 * n) - 1;
		mkey = (p->mkey >> 2);
		if (m != 16)
	   		mkey &= (1 << 2 * m) - 1;
		return find_kmer(mkey + (key << 2 * (m - 1)), nkey, T);
	}
	else
	{
		unsigned int mkey = 0;
		unsigned int nkey = 0;
		unsigned int temp = (p->nkey >> 2 * (n - 1)) & 3;
		mkey = (p->mkey << 2) + temp;
		if (m != 16)
			mkey &= (1 << 2 * m) - 1;
		nkey = (p->nkey << 2);
		if (n != 16)
			nkey &= (1 << 2 * n) - 1;
		return find_kmer(mkey, nkey + key, T);
	}
}
	

//form MUP based on kmer pointed by pointer p
//In other words, find MUP that pass kmer pointed by p.
//we define MUP which length below mup_lnth0 as complex region. complex region is region with many branches.
//We make complex region as MUP_graph[0] before use Make_MUP
void Make_MUP(KmerTable &T, MUP_Graph &S, edge * p, int mup_lnth0)	
{
	if (p == 0 || p->isdel == 1)
		return;
	if (p->MUP_ind != -1)
		return;
	MUP mup;
	init_MUP(mup);
	int index = S.size();
	int lnth = 1;
	edge * q = 0;
	edge * rq = 0;
	q = p;
	int rnum = 0;
	int lrnum = 0;
	rnum = get_kmer_nei_num(q, 0);
	if (rnum == 1)
	{
		for (int i = 4; i != 8; ++i)
			if (rq = get_kmer_nei(q, T, i))
				break;
		lrnum = get_kmer_nei_num(rq, 1);
	}
	while ( rnum == 1 && lrnum == 1 && rq->MUP_ind == -1)
//keep the path unique and kmer is not used
	{
		q->MUP_ind = index;
		q = rq;
		lnth++;
		rnum = get_kmer_nei_num(q, 0);
		if (rnum != 1)
			break;
		for (int i = 4; i != 8; ++i)
			if (rq = get_kmer_nei(q, T, i))
				break;
		lrnum = get_kmer_nei_num(rq, 1);
	} 
	q->MUP_ind = index;
	mup.end = q;
	q = p;
	q->MUP_ind = -1;
	edge * lq = 0;
	int lnum = 0;
	int rlnum = 0;
	lnum = get_kmer_nei_num(q, 1);
	if (lnum == 1)
	{
		for (int i = 0; i != 4; ++i)
			if (lq = get_kmer_nei(q, T, i))
				break;
		rlnum = get_kmer_nei_num(lq, 0);
	}
	while (lnum == 1 && rlnum == 1 && lq->MUP_ind == -1)
//keep the path unique and kmer is not used
	{
		q->MUP_ind = index;
		q = lq;
		lnth++;
		lnum = get_kmer_nei_num(q, 1);
		if (lnum != 1)
			break;
		for (int i = 0; i != 4; ++i)
			if (lq = get_kmer_nei(q, T, i))
				break;
		rlnum = get_kmer_nei_num(lq, 0);
	} 
	q->MUP_ind = index;
	if (lnth < mup_lnth0)	//for short mup, we see it as bubble
		//bubble region is as complex region, so we make MUP_ind zero and make these kmers belong to complex region
	{
		while (q != mup.end)
		{
			q->MUP_ind = 0;
			for (int i = 4; i != 8; ++i)
				if (rq = get_kmer_nei(q, T, i))
					break;
			q = rq;
		}
		q->MUP_ind = 0;
	}
	else
	{
		mup.start = q;
		mup.length = lnth;	//lnth represents edges the MUP passes
		S.push_back(mup);
	}
	return;
}


//save all MUPs to MUP_graph
void Form_MUP(KmerTable &T, MUP_Graph &S, int lnth0)	
{
	//generate the complex region as MUP_graph[0]
	MUP mup0;
	init_MUP(mup0);
	S.push_back(mup0);
	//generate the complex region as MUP_graph[0]
	tree::iterator itr_end;
	for (size_t i = 0; i != T.size(); i++)
	{
		if (T[i].size() > 0)
		{
			itr_end = T[i].end();
			for (tree::iterator itr = T[i].begin(); itr != itr_end; itr++)
			{
				if (itr->second.MUP_ind == -1 && itr->second.isdel != 1)	
//to check if the kmer has already used or has been deleted
					Make_MUP(T, S, &itr->second, lnth0);
			}
		}
	}
	return;
}

//will be delete
string Print_MUP(MUP * p, KmerTable &T)
{
	edge * e = 0;
	edge * np = 0;
	unsigned int key = 0;
	char c = 0;
	string s;
	e = p->start;
	while (e != p->end)
	{
		key = e->nkey & 3;
		s += key_to_char(key);
		for (int i = 4; i != 8; ++i)
		{
			if (np = get_kmer_nei(e, T, i))
				break;
		}
		e = np;
	}
	key = p->end->nkey & 3;
	s += key_to_char(key);
	return s;
}
//will be delete!


//to get the kmer which hash value is mkey and nkey,
//	and save it to kmer which is a c_string
//It can't verify if it is a valid hash value!
string inv_hash(unsigned int mkey, unsigned int nkey, int k, int n)
{
	string s;
	unsigned int t = 0;
	int m = k - n;
	for (int i = 0; i != m; ++i)
	{
		t = (mkey >> 2 * (m - 1 - i)) & 3;
		s += key_to_char(t);
	}
	for (int i = 0; i != n - 1; i++)
	{
		t = (nkey >> 2 * (n - 1 - i)) & 3;
		s += key_to_char(t);
	}
	return s;
}

string Print_TotalMUP(MUP * p, KmerTable &T, int k, int n)
{
	string s;
	edge * e = p->start;
	s = inv_hash(e->mkey, e->nkey, k, n);
	s += Print_MUP(p, T);
	return s;
}


//return -1 if get_startpos could not find the position
int get_startpos(reads * s, unsigned int mkey, unsigned int nkey, int k, int n)
{
	if (s->size < k)
		return 10000;
	unsigned int mkey1 = 0;
	unsigned int nkey1 = 0;
	int l = (k - 1) / 16;
	int t = 0;
	unsigned int temp = 0;
	unsigned int startpos = 0;
	kmer2hash(s, mkey1, nkey1, startpos, k, n);
	if (mkey1 == mkey && nkey1 == nkey)
		return startpos;
	for (int i = k; i != s->size; ++i)
	{
		++startpos;
		t = i % 16;
		if (t == 0)
			l++;
		temp = *(s->p + l);
		temp &= get2bit[t];
		temp >>= 2 * (15 - t);
		hash_kmer(temp, mkey1, nkey1, k, n);
		if (mkey1 == mkey && nkey1 == nkey)
			return startpos;
	}
	return -1;
}

inline edge * complement_kmer(KmerTable &T, edge * p)
{
	unsigned int mkey = 0;
	unsigned int nkey = 0;
	mkey = p->mkey;
	nkey = p->nkey;
	int m = globalArgs.k - globalArgs.n;
	if (globalArgs.n != 16)
		nkey = (1 << 2 * globalArgs.n) - 1 - nkey;
	else
		nkey = ~nkey;
	if (m != 16)
		mkey = (1 << 2 * m) - 1 - mkey;
	else
		mkey = ~mkey;
	return &T[mkey][nkey];
}


//build right neighbour for complex region
void mup0_rneiBuild(KmerTable &T, ReadTable V, MUP_Graph &S, int k, int n, int dist, unsigned int read_cov0)
{
	tree::iterator itr;
	tree::iterator itr_end;
	map<int, set<int> > record;
	if (T.size() == 0)
		return;
	for (size_t i = 0; i != T.size(); i++)
	{
		if (T[i].size() > 0)
		{
			itr_end = T[i].end();
			for (itr = T[i].begin(); itr != itr_end; itr++)
				if (itr->second.isdel == 0 && itr->second.MUP_ind == 0)
				{
					edge * e = &itr->second;
					from_read * f = e->head;
					unsigned int mkey = 0;
					unsigned int nkey = 0;
					if (f == 0)
					{
						edge * ce = complement_kmer(T, e);
						f = ce->head;
						while (f != 0)
						{
							unsigned int read_ind = f->i;
							reads * q = V + read_ind;
							complement(q);
							int startpos = get_startpos(q, e->mkey, e->nkey, k, n);
							if (startpos == -1)
							{
								complement(q);
								f = f->next;
								continue;
							}
							if (kmer2hash(q, mkey, nkey, startpos + dist, k, n)) 
//success find am if Get_AM return 1
							{	
								edge * p = &T[mkey][nkey];
								if (p->MUP_ind == 0 || p->isdel != 0)
								{
									complement(q);
									f = f->next;
									continue;
								}
								record[p->MUP_ind].insert(read_ind);
							}
							complement(q);
							f = f->next;
						}
					}
					else
						while (f != NULL)
						{
							unsigned int read_ind = f->i;
							reads * q = V + read_ind;
							int startpos = get_startpos(q, e->mkey, e->nkey, k, n);
							if (startpos == -1)
							{
								f = f->next;
								continue;
							}
							if (kmer2hash(q, mkey, nkey, startpos + dist, k, n)) 
//success find am if Get_AM return 1
							{
								edge * p = &T[mkey][nkey];
								if (p->MUP_ind == 0 || p->isdel != 0)
								{
									f = f->next;
									continue;
								}
								record[p->MUP_ind].insert(read_ind);
							}
							f = f->next;
						}	
					}
		}
	}
	map<int, set<int> >::iterator it;
	map<int, set<int> >::iterator it_end;
	it_end = record.end();
	for (it = record.begin(); it != it_end; it++)
		if (it->second.size() >= read_cov0)    
//it is as a valid neighbour of complex region if its coverage is larger than read_cov0
		{
			MUP_lnei * p = S[it->first].lnei;
			if (p == 0)
			{
				p = (MUP_lnei*)malloc(sizeof(MUP_lnei));
				if (p == 0)
				{
					printf("short of memory!\n");
					return;
				}
				p->p = &S[0];
				p->next = 0;
				S[it->first].lnei = p;
			}
		}
	return;
}


//compute the length between two MUPs
//one is from kmer "edge * e", the other is from "edge * am";
//dist is the distance of AM
//MUP included "edge *am" is on the right from MUP included "edge * e"
//we don't computer the distance if two kmer are from the same MUP or one of them is from complex region; in that case we return -1
//here the distance means move steps from one end to other's start, that is to say if two MUPs's distance is one, they are natural connected.
int Get_gaplnth(KmerTable &T, MUP_Graph &S, edge * e, edge * am, int dist)	
{
	if (e->MUP_ind == am->MUP_ind)
		return -1;
	if (e->MUP_ind == 0 || am->MUP_ind == 0)
		return -1;
	int d1 = 0;
	int d2 = 0;
	edge * e1 = S[e->MUP_ind].end;
	edge * e2 = S[am->MUP_ind].start;
	edge * next = 0;
	while (e != e1)
	{
		d1++;
		for (int i = 4; i != 8; ++i)
			if (next = get_kmer_nei(e, T, i))
				break;
		e = next;
	}
	while (am != e2)
	{
		d2++;
		for (int i = 0; i != 4; ++i)
			if (next = get_kmer_nei(am, T, i))
				break;
		am = next;
	}
	return dist - d1 - d2;
}


//isleft is to decide whether the insertion should happen in left(isleft = 1) or right(isleft = 0)
int mup_nei_insert(MUP * s, MUP * des, int isleft, short dist = 1)
{
	MUP * p = s;
	if (isleft)
	{
		if (p->lnei == 0)
		{
			p->lnei = (MUP_lnei*)malloc(sizeof(MUP_lnei));
			if (p->lnei == 0)
			{
				printf("short of memory\n");
				return 0;
			}
			p->lnei->p = des;
			p->lnei->next = 0;
			return 1;
		}
		else
		{
			MUP_lnei * q = p->lnei;
			while (q->next != 0)
			{
				if (q->p == des)
					return 0;
				q = q->next;
			}
			if (q->p == des)
				return 0;
			q->next = (MUP_lnei*)malloc(sizeof(MUP_lnei));
			if (q->next == 0)
			{
				printf("short of memory\n");
				return 0;
			}
			q->next->p = des;
			q->next->next = 0;
			return 1;
		}
	}
	else
	{
		if (p->rnei == 0)
		{
			p->rnei = (MUP_rnei*)malloc(sizeof(MUP_rnei));
			if (p->rnei == 0)
			{
				printf("short of memory\n");
				return 0;
			}
			p->rnei->p = des;
			p->rnei->dist = dist;
			p->rnei->next = 0;
			return 1;
		}
		else
		{
			MUP_rnei * q = p->rnei;
			while (q->next != 0)
			{
				if (q->p == des)
					return 0;
				q = q->next;
			}
			if (q->p == des)
				return 0;
			q->next = (MUP_rnei*)malloc(sizeof(MUP_rnei));
			if (q->next == 0)
			{
				printf("short of memory\n");
				return 0;
			}
			q->next->p = des;
			q->next->dist = dist;
			q->next->next = 0;
			return 1;
		}
	}
}


//pair.first record number of reads support dist 1 between MUPs
//pair.second record number of reads support dist 2 to 1 + dist_error between MUPs
typedef pair<int, int> r_cov;


//check links between MUP based on read_cov, int cov is the coverage of lnei links to RBTR;
//read_cov is the number of supported reads for two related MUPs
//cov is the number of supported reads for link between MUP and complex region
//"int index" is the index of MUP
void check_MUP(MUP_Graph &S, MUP * index, map<int, r_cov> &record, int read_cov, int read_cov0, int cov)
{
	MUP * p = index;
	if (cov >= read_cov0)
	{
		mup_nei_insert(p, &S[0], 0, -1);
	}
	if (record.size() == 0)
		return;
	map<int, r_cov>::iterator itr;
	for (itr = record.begin(); itr != record.end(); ++itr)
	{
		if (itr->second.first >= read_cov)
		{
			mup_nei_insert(p, &S[itr->first], 0);
			MUP * q = &S[itr->first];
			mup_nei_insert(q, index, 1);
		}
		else if (itr->second.first == read_cov - 1)
		{
			if (itr->second.second > 0)
			{
				mup_nei_insert(p, &S[itr->first], 0);
				MUP * q = &S[itr->first];
				mup_nei_insert(q, index, 1);
			}
		}
	}
	return;
}


//build MUP_DeBrujin Graph
void Form_MUP_Graph(KmerTable &T, ReadTable V, MUP_Graph &S, int k, int n, int dist, int read_cov, int read_cov0, int dist_error)	
{
//find right neighbour for complex region
	mup0_rneiBuild(T, V, S, k, n, dist, read_cov0);
	for (size_t j = 1; j != S.size(); j++)
	{
		edge * lp = 0;
		map<int, r_cov> record;
		edge * e = S[j].end;
//search_range has been defined in #define search_range
		int t = search_range;	//to indicate the range of kmer being searched to form MUP
		int cov0 = 0;
		set<unsigned int> cov;
		while (t > 0)
		{
			if (e == 0)
				return;
			from_read * f = e->head;
			unsigned int mkey = 0;
			unsigned int nkey = 0;
			if (f == 0)
			{
				edge * ce = complement_kmer(T, e);
				f = ce->head;
				while (f != 0)
				{
					unsigned int read_ind = f->i;
					reads * q = V + read_ind;
					complement(q);
					int startpos = get_startpos(q, e->mkey, e->nkey, k, n);
					if (startpos == -1)
					{
						complement(q);
						f = f->next;
						continue;
					}
					if (kmer2hash(q, mkey, nkey, startpos + dist, k, n)) 
//success find am if Kmer2hash return 1
					{
						edge * am = &T[mkey][nkey];
						int mup_ind = am->MUP_ind;
						if (am->isdel == 0)
						{
							if (mup_ind == 0)	
//we don't compute the distance between complex region
								cov.insert(read_ind);
							else
							{
								int length = Get_gaplnth(T, S, e, am, dist);
								if (length > 0)
								{	
									if (length == 1)
										record[mup_ind].first++;
									else if (length <= 1 + dist_error && length > 1)
										record[mup_ind].second++;
								}
							}
						}
					}
					complement(q);
					f = f->next;
				}
			}
			else
			{
				while (f != 0)
				{   
					unsigned int read_ind = f->i;
					reads * q = V + read_ind;
					int startpos = get_startpos(q, e->mkey, e->nkey, k, n);
					if (startpos == -1)
					{
						f = f->next;
						continue;
					}
					if (kmer2hash(q, mkey, nkey, startpos + dist, k, n)) 
//success find am if Kmer2hash return 1
					{
						edge * am = &T[mkey][nkey];
						int mup_ind = am->MUP_ind;
						if (am->isdel == 0)
						{
							if (mup_ind == 0)	
//we don't compute the distance between complex region
								cov.insert(read_ind);
							else
							{
								int length = Get_gaplnth(T, S, e, am, dist);
								if (length > 0)
								{
									if (length == 1)
										record[mup_ind].first++;
									else if (length <= 1 + dist_error && length > 1)
										record[mup_ind].second++;
								}
							}
						}		
					}
					f = f->next;
				}
			}
			t--;
			for (int i = 0; i != 4; ++i)
				if (lp = get_kmer_nei(e, T, i))
					break;
			e = lp;
			cov0 = cov.size();
		}
		check_MUP(S, &S[j], record, read_cov, read_cov0, cov0);
	}
	return;
}


//for test purpose
int find_MUP_index(MUP_Graph &S, MUP * p)
{
	for (size_t i = 0; i != S.size(); i++)
		if (p == &S[i])
			return i;
	return -1;
}

int get_LongMUP_id(LongMUP * p, LongMUP_Graph &L)
{
	for (size_t i = 0; i != L.size(); ++i)
	{
		if (p == &L[i])
			return i;
	}
	return -1;
}
//for test purpose

typedef map<int, map<int, set<int> > > MUPlink_Record;

struct Best_MUPlink
{
	MUP * self;
	int des;
	unsigned int freq;
	int length;
};
typedef struct Best_MUPlink Best_MUPlink;


inline MUP_lnei * get_MUP_lnei0(MUP * p, const MUP * mup0)
{
	MUP_lnei * lp = p->lnei;
	if (lp == 0)
		return lp;
	while (lp != 0)
	{
		if (lp->p == mup0)
			return lp;
		lp = lp->next;
	}
	return 0;
}

inline MUP_rnei * get_MUP_rnei0(MUP * p, const MUP * mup0)
{
	MUP_rnei * rp = p->rnei;
	if (rp == 0)
		return rp;
	while (rp != 0)
	{
		if (rp->p == mup0)
			return rp;
		rp = rp->next;
	}
	return 0;
}


void del_MUP_connect(MUP * p, MUP_rnei * rp)
{
	if (rp == 0)
		return;
	MUP_rnei * np = p->rnei;
	if (np == 0)
		return;
	if (np == rp)
	{
		p->rnei = rp->next;
		free(rp);
		rp = 0;
	}
	else
	{
		while (np->next != rp)
		{
			if (np->next == 0)
				return;
			np = np->next;
		}
		np->next = rp->next;
		free(rp);
		rp = 0;
	}
	return;
}

void del_MUP_connect(MUP * p, MUP_lnei * lp)
{
	if (lp == 0)
		return;
	MUP_lnei * np = p->lnei;
	if (np == 0)
		return;
	if (np == lp)
	{
		p->lnei = lp->next;
		free(lp);
		lp = 0;
	}
	else
	{
		while (np->next != lp)
		{
			if (np->next == 0)
				return;
			np = np->next;
		}
		np->next = lp->next;
		free(lp);
		lp = 0;
	}
	return;
}

void check_bubble(KmerTable &T, MUP_Graph &S, ReadTable V, MUP * bubble, MUPlink_Record &record, int dist)
{
	int k = globalArgs.k;
	int n = globalArgs.n;
	edge * e = bubble->end;
//search_range has been defined in #define search_range
	int t = search_range;	
//to indicate the range of kmer being searched to form MUP
	while (t > 0)
	{
		if (e == 0)
			return;
		from_read * f = e->head;
		unsigned int mkey = 0;
		unsigned int nkey = 0;
		if (f == 0)
		{
			edge * ce = complement_kmer(T, e);
			f = ce->head;
			while (f != 0)
			{
				unsigned int read_ind = f->i;
				reads * q = V + read_ind;
				complement(q);
				int startpos = get_startpos(q, e->mkey, e->nkey, k, n);
				if (startpos == -1)
				{
					complement(q);
					f = f->next;
					continue;
				}
				if (kmer2hash(q, mkey, nkey, startpos + dist, k, n)) 
//success find am if Kmer2hash return 1
				{
					edge * am = &T[mkey][nkey];
					int mup_ind = am->MUP_ind;
					MUP_lnei * lp = 0;
					if (am->isdel == 0)
					{
						if (mup_ind != 0 && mup_ind != e->MUP_ind)
							if (get_MUP_lnei0(&S[mup_ind], &S[0]))
//we don't compute the distance between complex region
                               {
								int length = Get_gaplnth(T, S, e, am, dist);
								if (length > 1)
									record[mup_ind][length].insert(read_ind);
							   }
					}
				}
				complement(q);
				f = f->next;
			}
		}
		else
		{
			while (f != 0)
			{   
				unsigned int read_ind = f->i;
				reads * q = V + read_ind;
				int startpos = get_startpos(q, e->mkey, e->nkey, k, n);
				if (startpos == -1)
				{
					f = f->next;
					continue;
				}
				if (kmer2hash(q, mkey, nkey, startpos + dist, k, n)) 
//success find am if Kmer2hash return 1
				{
					edge * am = &T[mkey][nkey];
					int mup_ind = am->MUP_ind;
					if (am->isdel == 0)
					{
						if (mup_ind != 0 && mup_ind != e->MUP_ind)	
							if (get_MUP_lnei0(&S[mup_ind], &S[0]))
//we don't compute the distance between complex region
							{
								int length = Get_gaplnth(T, S, e, am, dist);
								if (length > 1)
									record[mup_ind][length].insert(read_ind);
							}
					}
				}
				f = f->next;
			}
		}
		t--;
		edge * le = 0;
		for (int i = 0; i != 4; ++i)
			if (le = get_kmer_nei(e, T, i))
				break;
		e = le;
	}
	return;
}

int likely_dist_MUPlink(MUPlink_Record &record, MUP * p, Best_MUPlink &result)
{
	if (record.size() == 0)
		return 0;
	MUPlink_Record::iterator itr;
	map<int, set<int> >::iterator it;
	result.freq = 0;
	result.des = 0;
	result.length = 0;
	result.self = p;
	for (itr = record.begin(); itr != record.end(); ++itr)
	{
		if (itr->second.size() > 0)
		{
			for (it = itr->second.begin(); it != itr->second.end(); ++it)
			{
				if (it->second.size() > result.freq)
				{
					result.freq = it->second.size();
					result.des = itr->first;
					result.length = it->first;
				}
			}
		}
	}
	return result.length;
}

void deal_bubble(MUP_Graph &S, MUPlink_Record &record, Best_MUPlink &result, unsigned int complex_th)
{
	if (result.freq < complex_th)
		return;
	MUP * p = result.self;
	MUP * q = &S[result.des];
	MUP_rnei * rp = p->rnei;
	MUP_lnei * lp = 0;
	if (rp == 0)
	   return;
	if (!(rp = get_MUP_rnei0(p, &S[0])))
		return;
	if (S[result.des].lnei == 0 || !(lp = get_MUP_lnei0(q, &S[0])))
		return;
	MUP_rnei * r = p->rnei;
	while (r != 0)
	{
		if (r->p == q)
		{
			del_MUP_connect(p, rp);
			break;
		}
		r = r->next;
	}
	r = rp;
	MUP_lnei * l = p->lnei;
	while (l != 0)
	{
		if (l->p == p)
		{
			del_MUP_connect(q, lp);
			break;
		}
		l = l->next;
	}
	l = lp;
	if (r == 0 || l == 0)
		return;
	r->p = q;
	r->dist = result.length;
	l->p = p;
}

void deal_complexregion(KmerTable &T, MUP_Graph &S, ReadTable pV, unsigned int complex_th, int dist)
{
	int k = globalArgs.k;
	int n = globalArgs.n;
	for (size_t i = 0; i != S.size(); ++i)
	{
		MUP * p = &S[i];
		if (p->rnei != 0 && get_MUP_rnei0(p, &S[0]))
		{
			MUPlink_Record record;
			Best_MUPlink result;
			int length = 0;
			check_bubble(T, S, pV, &S[i], record, dist);
			length = likely_dist_MUPlink(record, p, result);
			if (length <= 0 || length > dist)
				continue;
            deal_bubble(S, record, result, complex_th);
		}
	}
	return;
}


void LongMUP_init(LongMUP * p, MUP * start = 0, MUP * end = 0, int length = -1, LongMUP_lnei * lnei = 0, LongMUP_rnei * rnei = 0)
{
	p->start = start;
	p->end = end;
	p->length = length;
	p->lnei = lnei;
	p->rnei = rnei;
	p->isused = 0;
}

int get_MUP_nei_num(MUP * p, int isleft)
{
	if (p == 0)
		return -1;
	int count = 0;
	if (isleft)
	{
		MUP_lnei * q = 0;
		q = p->lnei;
		while (q != 0)
		{
			++count;
			q = q->next;
		}
		return count;
	}
	else
	{
		MUP_rnei * q = 0;
		q = p->rnei;
		while (q != 0)
		{
			++count;
			q = q->next;
		}
		return count;
	}
}

int left_MUP_dist(MUP * left, MUP * right)
{
	if (left == 0 || right == 0)
		return -1;
	MUP_rnei * np = left->rnei;
	if (np == 0)
		return -1;
	while (np != 0)
	{
		if (np->p == right)
			return np->dist;
		np = np->next;
	}
	return -1;
}

void Preprocess1LongMUP(MUP * p, const MUP * terminal)
{
	if (p == 0)
		return;
	MUP_rnei * rp = p->rnei;
	int num = get_MUP_nei_num(p, 0);
	if (num >= 2)
	{
		while (rp != 0)
		{
			if (rp->p == terminal)
			{
				del_MUP_connect(p, rp);
				break;
			}
			rp = rp->next;
		}
	}
	MUP_lnei * lp = p->lnei;
	num = get_MUP_nei_num(p, 1);
	if (num >= 2)
	{
		while (lp != 0)
		{
			if (lp->p == terminal)
			{
				del_MUP_connect(p, lp);
				break;
			}
			lp = lp->next;
		}
	}
	return;
}


void Preprocess_LongMUP(MUP_Graph &S)
{
	for (int i = 1; i != (int)S.size(); ++i)
		Preprocess1LongMUP(&S[i], &S[0]);
	return;
}


void extend_MUP(LongMUP_Graph &L, MUP * pmup, const MUP * terminal)
//ter is &S[0], it refers to complex region
{
	if (pmup == 0 && pmup->index != -1)
		return;
	if (pmup == terminal)
		return;
	LongMUP lmup;
	LongMUP_init(&lmup);
	MUP * p = pmup;
	MUP_rnei * q = pmup->rnei;
	int length = p->length;
	int rnum = 0;
	int lnum = 0;
	if (q != 0)
	{
		rnum = get_MUP_nei_num(p, 0);
		lnum = get_MUP_nei_num(q->p, 1);
		while (rnum == 1 && lnum == 1 && q->p->index == -1)
		{
			if (q->p != terminal)
			{
				p->index = L.size();
				length += q->dist - 1;
				p = q->p;
				q = p->rnei;
				length += p->length;
				if (q == 0)
					break;
				rnum = get_MUP_nei_num(p, 0);
				lnum = get_MUP_nei_num(q->p, 1);
			}
		}
	}
	p->index = L.size();
	lmup.end = p;
	p = pmup;
	p->index = -1;
	MUP_lnei * q1 = p->lnei;
	if (q1 != 0)
	{
		lnum = get_MUP_nei_num(p, 1);
		rnum = get_MUP_nei_num(q1->p, 0);
		while (rnum == 1 && lnum == 1 && q1->p->index == -1)
		{
			if (q1->p != terminal)
			{
				int dist = 0;
				p->index = L.size();
				if ((dist = left_MUP_dist(q1->p, p)) <= 0)
					return;
				length += dist - 1;
				p = q1->p;
				q1 = p->lnei;
				length += p->length;
			}
			if (q1 == 0)
				break;
			lnum = get_MUP_nei_num(p, 1);
			rnum = get_MUP_nei_num(q1->p, 0);
		}
	}
	p->index = L.size();
	lmup.start = p;
	lmup.length = length;
	L.push_back(lmup);
	return;
}

void Form_LongMUP(LongMUP_Graph &L, MUP_Graph &S)
{
	LongMUP lmup0;
	LongMUP_init(&lmup0);
	L.push_back(lmup0);
	if (S.size() == 1)
	{
		printf("No MUP found!\n");
		return;
	}
	for (unsigned int i = 1; i != S.size(); i++)
	{
		if (S[i].index != -1)
			continue;
		MUP * p = &S[i];
		extend_MUP(L, p, &S[0]);
	}
	return;
}

int dist_LongMUP(LongMUP * p1, LongMUP * p2)
{
	if (p1 == 0 || p2 == 0)
		return -1;
	MUP * p = p1->end;
	MUP_rnei * np = p->rnei;
	if (np == 0)
		return -1;
	while (np != 0)
	{
		if (np->p == p2->start)
			return np->dist;
		np = np->next;
	}
	return -1;
}

int Creat_LongMUP_nei(LongMUP * p, LongMUP * p1, int isleft)
{
	if (isleft)
	{
		LongMUP_lnei * q = p->lnei;
		if (q == 0)
		{
			q = (LongMUP_lnei*)malloc(sizeof(LongMUP_lnei));
			q->p = p1;
			q->next = 0;
			p->lnei = q;
			return 1;
		}
		else
		{
			if (q->p == p1)
				return 0;
			while (q->next != 0)
			{
				if (q->next->p == p1)
					return 0;
				q = q->next;
			}
			q->next = (LongMUP_lnei *)malloc(sizeof(LongMUP_lnei));
			q->next->p = p1;
			q->next->next = 0;
			return 1;
		}
	}
	else
	{
		LongMUP_rnei * q = p->rnei;
		if (q == 0)
		{
			q = (LongMUP_rnei *)malloc(sizeof(LongMUP_rnei));
			q->p = p1;
			q->next = 0;
			q->dist = dist_LongMUP(p, p1);
			p->rnei = q;
			return 1;
		}
		else
		{
			if (q->p == p1)
				return 0;
			while (q->next != 0)
			{
	        	if (q->p == p1)
					return 0;
				q = q->next;
			}
			q->next = (LongMUP_rnei *)malloc(sizeof(LongMUP_rnei));
			q->next->p = p1;
			q->next->dist = dist_LongMUP(p, p1);
			q->next->next = 0;
			return 1;
		}
	}
}


//LongMUP 1, 2 can link if 1's end and 2's start are connected as MUP.
//np is the neighboor of p->end
//here we only consider simple region connection. In other word,
//we remove connection structure in MUP_Graph to form LongMUP_Graph
void Link_LongMUP(LongMUP * p, MUP_rnei * np, LongMUP_Graph &L, MUP_Graph &S)
{
	LongMUP * p1 = 0;
	if (np == 0 || p == 0)
		return;
	while (np != 0)
	{
		int i = np->p->index;
		if (i >0)
		{
			if (L[i].start != np->p || p == &L[i])
			{
				np = np->next;
				continue;
			}
			p1 = &L[i];
			if (Creat_LongMUP_nei(p, p1, 0))
				Creat_LongMUP_nei(p1, p, 1);
			np = np->next;
		}
		return;
	}
}


void Form_LongMUP_Graph(LongMUP_Graph &L, MUP_Graph &S)
{
	if (L.size() == 0)
	{
		printf("No Long MUP exist!\n");
		return;
	}
	LongMUP * p = 0;
	MUP * pmup = 0;
	MUP_rnei * head = 0;
	for (size_t i = 1; i != L.size(); ++i)
	{
		p = &L[i];
		pmup = p->end;
		head = pmup->rnei;
		Link_LongMUP(p, head, L, S);
	}
	return;
}

int count_LongMUP_nei(LongMUP * p, int isleft)
{
	if (p == 0)
		return -1;
	int count = 0;
	if (isleft)
	{
		LongMUP_lnei * q = 0;
		q = p->lnei;
		while (q != 0)
		{
			++count;
			q = q->next;
		}
		return count;
	}
	else
	{
		LongMUP_rnei * q = 0;
		q = p->rnei;
		while (q != 0)
		{
			++count;
			q = q->next;
		}
		return count;
	}
}

void LongMUP_clear(LongMUP_Graph &L)
{
	for (int i = 0; i != L.size(); ++i)
		L[i].isused = 0;
	return;
}

void find_LongMUP_beginpoint(LongMUP_Graph &L, set<int> &start)
{
	if (L.size() == 0)
		return;
	LongMUP * p = 0;
	int num = 0;
//L[0] represents the complex region in LongMUP structure
	for (size_t i = 1; i != L.size(); ++i)
	{
		p = &L[i];
		if (p->lnei == 0)
		{
			start.insert(i);
			continue;
		}
		LongMUP_lnei * q = p->lnei;
		num = count_LongMUP_nei(q->p, 0);
		if (num > 1)
		{
			start.insert(i);
			break;
		}
	}
	return;
}

string fill_MUPgap(KmerTable &T, ReadTable V, MUP * p, MUP_rnei * rp)
{
	int k = globalArgs.k;
	int n = globalArgs.n;
	unsigned int key = 0;
	string s;
	if (rp->dist <= 1)
		return s;
	if (p == 0)
		return s;
	edge * e = p->end;
	edge * lp = 0;
	int t = search_range;
	while (t > 0)
	{
		if (e == 0)
			return s;
		from_read * f = e->head;
		unsigned int mkey = 0;
		unsigned int nkey = 0;
		if (f == 0)
		{
				edge * ce = complement_kmer(T, e);
				f = ce->head;
				while (f != 0)
				{
					unsigned int read_ind = f->i;
					reads * q = V + read_ind;
					complement(q);
					int startpos = get_startpos(q, e->mkey, e->nkey, k, n);
					if (startpos == -1)
					{
						complement(q);
						f = f->next;
						continue;
					}
					if (kmer2hash(q, mkey, nkey, startpos + rp->dist, k, n)) 
//success find am if Kmer2hash return 1
					{
						edge * am = &T[mkey][nkey];
						if (am == rp->p->start)
						{
							int i = 1;
							while (i < rp->dist)
							{
								if (kmer2hash(q, mkey, nkey, startpos + i, k, n))
								{
									if (T[mkey][nkey].isdel != 0)
										break;
									s += key_to_char(nkey & 3);
									++i;
								}
							}
							if (i == rp->dist)
							{
								complement(q);
								return s;
							}
						}
					}
					complement(q);
					f = f->next;
				}
			}
			else
			{
				while (f != 0)
				{   
					unsigned int read_ind = f->i;
					reads * q = V + read_ind;
					int startpos = get_startpos(q, e->mkey, e->nkey, k, n);
					if (startpos == -1)
					{
						f = f->next;
						continue;
					}
					if (kmer2hash(q, mkey, nkey, startpos + rp->dist, k, n)) 
//success find am if Kmer2hash return 1
					{
						edge * am = &T[mkey][nkey];
						if (am == rp->p->start)
						{
							int i = 1;
							while (i < rp->dist)
							{
								if (kmer2hash(q, mkey, nkey, startpos + i, k, n))
								{
									if (T[mkey][nkey].isdel != 0)
										break;
									s += key_to_char(nkey & 3);
									++i;
								}
							}
							if (i == rp->dist)
								return s;
						}
					}
					f = f->next;
				}
			}
			t--;
			for (int i = 0; i != 4; ++i)
				if (lp = get_kmer_nei(e, T, i))
					break;
			e = lp;
	}
	s.clear();
	int i = 1;
	while (i != rp->dist)
		s += 'N';
	return s;
}

string Print_LongMUP(KmerTable &T, ReadTable &V, MUP_Graph &S, LongMUP * target)
{
	MUP * p = target->start;
	string s;
	unsigned int key = 0;
	while (p != target->end)
	{
		s += Print_MUP(p, T);
		s += fill_MUPgap(T, V, p, p->rnei);
		p = p->rnei->p;
	}
	s += Print_MUP(target->end, T);
	return s;
}

inline int LongMUPrnei_dist(LongMUP * p1, LongMUP * p2)
{
	LongMUP_rnei * p = p1->rnei;
	while (p != 0)
	{
		if (p->p == p2)
			return p->dist;
		p = p->next;
	}
	return -1;
}


//ContigPath_record records the current path
//sign is the stop point in the road
//pair_distance compute the distance between sign and record.back
//which don't count the length of sign.
//dist is the distance between to LongMUP
int pair_distance(ContigPath_record &record, LongMUP * sign, int dist)
{
	if (sign == 0)
		return -1;
	int length = 0;
	int key = -1;
	LongMUP * p = sign;
	for (int i = 0; i != (int)record.size(); ++i)
	{
		if (p == record[i])
		{
			key = i;
			break;
		}
	}
	if (key == -1)
		return -1;
	for (int i = key; i + 1 != (int)record.size(); ++i)
	{
		length += LongMUPrnei_dist(record[key], record[key + 1]) - 1;
		length += record[key + 1]->length;
	}
	length += dist - 1;
	return length;
}


int check_short(LongMUP * p1, LongMUP * p2, KmerTable &T, ReadTable pV, MUP_Graph &S, LongMUP_Graph &L, int length, int readcov)
{
	int k = globalArgs.k;
	int n = globalArgs.n;
	MUP * p = p1->end;
	edge * e = p->end;
	int t = search_range;
	set<int> cov;
	while (t > 0)
	{
		if (e == 0)
			return 0;
		from_read * f = e->head;
		unsigned int mkey = 0;
		unsigned int nkey = 0;
		if (f == 0)
		{
				edge * ce = complement_kmer(T, e);
				f = ce->head;
				while (f != 0)
				{
					unsigned int read_ind = f->i;
					if (cov.count(read_ind) > 0)
					{
						f = f->next;
						continue;
					}
					reads * q = pV + read_ind;
					complement(q);
					int startpos = get_startpos(q, e->mkey, e->nkey, k, n);
					if (startpos == -1)
					{
						complement(q);
						f = f->next;
						continue;
					}
					if (kmer2hash(q, mkey, nkey, startpos + length + search_range, k, n)) 
//success find am if Kmer2hash return 1
					{
						edge * edge_am = &T[mkey][nkey];
						MUP * MUP_am = &S[edge_am->MUP_ind];
						LongMUP * LongMUP_am = &L[MUP_am->index];
						if (LongMUP_am == p2)
							cov.insert(read_ind);
					}
					complement(q);
					f = f->next;
				}
			}
			else
			{
				while (f != 0)
				{   
					unsigned int read_ind = f->i;
					if (cov.count(read_ind) > 0)
					{
						f = f->next;
						continue;
					}
					reads * q = pV + read_ind;
					int startpos = get_startpos(q, e->mkey, e->nkey, k, n);
					if (startpos == -1)
					{
						f = f->next;
						continue;
					}
					if (kmer2hash(q, mkey, nkey, startpos + length + search_range, k, n)) 
//success find am if Kmer2hash return 1
					{
						edge * edge_am = &T[mkey][nkey];
						MUP * MUP_am = &S[edge_am->MUP_ind];
						LongMUP * LongMUP_am = &L[MUP_am->index];
						if (LongMUP_am == p2)
							cov.insert(read_ind);
					}
					f = f->next;
				}
			}
			--t;
			edge * lp = 0;
			for (int i = 0; i != 4; ++i)
				if (lp = get_kmer_nei(e, T, i))
					break;
			e = lp;
	}
	if (cov.size() >= (unsigned int)readcov)
		return 1;
	return 0;
}


int check_long(LongMUP * p1, LongMUP * p2, KmerTable &T, ReadTable pV, MUP_Graph &S, LongMUP_Graph &L, int length, unsigned int pairkey, int paircov)
{
	if (pairkey == 0)
		return 0;
	int k = globalArgs.k;
	int n = globalArgs.n;
	MUP * p = p1->end;
	edge * e = p->end;
	int t = search_range;
	set<int> cov;
	while (t > 0)
	{
		if (e == 0)
			return 0;
		from_read * f = e->head;
		unsigned int mkey = 0;
		unsigned int nkey = 0;
		if (f == 0)
		{
			edge * ce = complement_kmer(T, e);
			f = ce->head;
			while (f != 0)
			{
				unsigned int read_ind = f->i;
				if (cov.count(read_ind) > 0 || read_ind >= pairkey)
				{
					f = f->next;
					continue;
				}
				reads * q = pV + read_ind;
				complement(q);
				for (int i = 0; i != q->size - k + 1; ++i)
				{
					if (kmer2hash(q, mkey, nkey, i, k, n))
					{
						int mup_id = T[mkey][nkey].MUP_ind;
						int Lmup_id = S[mup_id].index;
						if (p2 == &L[Lmup_id])
						{
							cov.insert(read_ind);
							break;
						}
					}
				}
				complement(q);
				f = f->next;
			}
		}
		else
		{
			while (f != 0)
			{   
				unsigned int read_ind = f->i;
				if (cov.count(read_ind) > 0 || read_ind >= pairkey)
				{
					f = f->next;
					continue;
				}
				reads * q = pV + read_ind + pairkey;
				for (int i = 0; i != q->size - k + 1; ++i)
				{
					if (kmer2hash(q, mkey, nkey, i, k, n))
					{
						int mup_id = T[mkey][nkey].MUP_ind;
						int Lmup_id = S[mup_id].index;
						if (p2 == &L[Lmup_id])
						{
							cov.insert(read_ind);
							break;
						}
					}
					f = f->next;
				}
			}
		}
		--t;
		edge * lp = 0;
		for (int i = 0; i != 4; ++i)
			if (lp = get_kmer_nei(e, T, i))
				break;
		e = lp;
	}
	if (cov.size() >= (unsigned int)paircov)
		return 1;
	return 0;
}

int counter = 0;


//ContigPath_record records current path of contig composed by LongMUPs
//V records alternative path's positon which is from left
void ContigPath_Select(ContigPath_record record, vector<LongMUP*> V, KmerTable &T, ReadTable pV, MUP_Graph &S, LongMUP_Graph &L, \
	Contig_Pathway &w, unsigned int pairkey, int Max, int Min, int Readlen, int Readcov, int paircov)
{
	int k = globalArgs.k;
	int n = globalArgs.n;
	int num = 0;
	int length = 0;
	if (record.size() == 0)
		return;
	LongMUP_rnei * rp = 0;
	LongMUP * p = record.back();
	p->isused = 1;
	num = count_LongMUP_nei(p, 0);
	while (num == 1)
	{
		if (!p->rnei->p->isused)
		{
			V.push_back(p);
			p = p->rnei->p;
			record.push_back(p);
			p->isused = 1;
			num = count_LongMUP_nei(p, 0);
		}
		else
		{
			w.push_back(record);
			return;
		}
	}
	if (num <= 0)
	{
		w.push_back(record);
		return;
	}
	if (V.size() == 0)
	{
		rp = p->rnei;
		while (rp != 0)
		{
			record.push_back(rp->p);
			ContigPath_Select(record, V, T, pV, S, L, w, pairkey, Max, Min, Readlen, Readcov, paircov);
			++counter;
			record.pop_back();
			rp = rp->next;
		}
		return;
	}
	rp = p->rnei;
	LongMUP * pp = V.back();
	while (rp != 0)
	{
		length = pair_distance(record, pp, rp->dist);
		if (length == -1)
			return;
		if (length + k + 1 > Readlen && length + k + 1 < Min)
			w.push_back(record);
		else if (length + k + 1 <= Readlen)
		{
			if (check_short(pp, rp->p, T, pV, S, L, length, Readcov))
			{
				record.push_back(rp->p);
				V.pop_back();
				ContigPath_Select(record, V, T, pV, S, L, w, pairkey, Max, Min, Readlen, Readcov, paircov);
				record.pop_back();
				++counter;
			}
		}
		else if (length + k + 1 >= Min && length + k + 1 <= Max)
		{
			if (check_long(pp, rp->p, T, pV, S, L, length, pairkey, paircov))
			{
				record.push_back(rp->p);
				V.pop_back();
				ContigPath_Select(record, V, T, pV, S, L, w, pairkey, Max, Min, Readlen, Readcov, paircov);
				record.pop_back();
				++counter;
			}
		}
		else
		{
			w.push_back(record);
		}
		rp = rp->next;
	}
	return;
}

void check_contig_start(set<int> &startpoint, LongMUP_Graph &L, Contig_Pathway &w)
{
	if (w.size() == 0)
		return;
	for (size_t i = 0; i != w.size(); ++i)
	{
		for (size_t j = 0; j != w[i].size(); ++j)
		{
			int key = get_LongMUP_id(w[i][j], L);
			if (startpoint.count(key) > 0)
				startpoint.erase(key);
		}
	}
	return;
}

MUP_rnei * get_LongMUP_nei(LongMUP * lp, LongMUP * rp)
{
	MUP * p = lp->end;
	MUP_rnei * np = p->rnei;
	if (np == 0)
		return 0;
	while (np != 0)
	{
		if (np->p == rp->start)
			return np;
		np = np->next;
	}
	return 0;
}


string print1contig(KmerTable &T, ReadTable pV, MUP_Graph &S, LongMUP_Graph &L, ContigPath_record &r)
{
	int k = globalArgs.k;
	int n = globalArgs.n;
	string s;
	MUP_rnei * rp;
	if (r.size() == 0)
		return s;
	LongMUP * p = r[0];
	MUP * p0 = p->start;
	edge * e = p0->start;
	s = inv_hash(e->mkey, e->nkey, k, n);
	for (size_t i = 0; i != r.size(); ++i)
	{
		p = r[i];
		s += Print_LongMUP(T, pV, S, p);
		if (i + 1 < r.size())
		{
			rp = get_LongMUP_nei(p, r[i + 1]);
			if (rp == 0)
			{
				printf("error for print contig!\n");
				s.clear();
				return s;
			}
			s += fill_MUPgap(T, pV, p->end, rp);
		}
	}
	return s;
}

void Print_Contig(KmerTable &T, ReadTable pV, MUP_Graph &S, LongMUP_Graph &L, unsigned int pairkey, int Max, int Min, int Readlen, int Readcov, int paircov, ofstream &out)
{
	ofstream out1("temp1");
	string s;
	int num = 0;
	set<int> startpoint;
	find_LongMUP_beginpoint(L, startpoint);
	set<int> sub = startpoint;
	if (startpoint.size() == 0)
		return;
	set<int>::iterator itr;
	set<int>::iterator itr_end = startpoint.end();
	for (itr = startpoint.begin(); itr != itr_end; ++itr)
	{
		if (sub.count(*itr) > 0)
		{
			LongMUP_clear(L);
			counter = 0;
			out1 << *itr << "\t";
			Contig_Pathway P;
			ContigPath_record record;
			record.push_back(&L[*itr]);
			vector<LongMUP*> V;
			ContigPath_Select(record, V, T, pV, S, L, P, pairkey, Max, Min, Readlen, Readcov, paircov);
			if (P.size() > 0)
			{
				for (size_t i = 0; i != P.size(); ++i)
				{
					++num;
					s = print1contig(T, pV, S, L, P[i]);
					out << ">" << num << "\t" << s.size() << "\t";
					for (size_t j = 0; j != P[i].size(); ++j)
						out << P[i][j] << ":";
					out << endl;
					LongMUP_lnei * lp = 0;
					LongMUP * item = 0;
					lp = P[i][0]->lnei;
					out << "lnei" << "\t";
					while (lp != 0)
					{
						item = lp->p;
						out << get_LongMUP_id(item, L) << ":";
						lp = lp->next;
					}
					out << endl;
					LongMUP_rnei * rp = 0;
					rp = P[i].back()->rnei;
					out << "rnei" << "\t";
					while (rp != 0)
					{
						item = rp->p;
						out << get_LongMUP_id(item, L) << ":";
						rp = rp->next;
					}
					out << endl;
					out << s << endl;
				}
			}
			check_contig_start(sub, L, P);
			out1 << counter << endl;
		}
	}
	return;
}



int main(int argc, char * argv[])
{
	globalArgs.k = 27;
	globalArgs.n = 14;
	globalArgs.lnth0 = 5;
	globalArgs.lnth1 = 100;
	globalArgs.d1 = 9;
	globalArgs.d2 = 19;
	globalArgs.kmercov = 3;
	globalArgs.readcov = 3;
	globalArgs.readcov0 = 4;
	globalArgs.err = 2;
	globalArgs.paircov = 10;
	globalArgs.Min = 200;
	globalArgs.Max = 700;
	globalArgs.contig_info = NULL;
	globalArgs.output = "contig.txt";
	globalArgs.MUP_info = NULL;
	globalArgs.single = NULL;
	globalArgs.left = NULL;
	globalArgs.right = NULL;
	globalArgs.readlen = 90;
	char c = 0;
	static const struct option longOpts[]={
		{"K_mer_length", 1, NULL, 'k'},
		{"n_length", 1, NULL, 'n'},
		{"Min_MUP_length", 1, NULL, 'm'},
		{"Min_Contig_length", 1, NULL, 'L'},
		{"AM_short_dist", 1, NULL, 'd'},
		{"AM_long_dist", 1, NULL, 'D'},
		{"K_mer_cov", 1, NULL, 'c'},
		{"AM_cov", 1, NULL, 'C'},
		{"AM_cov_complex", 1, NULL, 'R'},
		{"alignment_err", 1, NULL, 'e'},
		{"output", 1, NULL, 'o'},
		{"Contig_info", 1, NULL, '1'},
		{"MUP_info", 1, NULL, '2'},
		{"single", 1, NULL, 's'},
		{"left", 1, NULL, 'l'},
		{"right", 1, NULL, 'r'},
		{"pair_cov", 1, NULL, 'p'},
		{"Min_fregment", 1, NULL, 'S'},
		{"Max_fregment", 1, NULL, 'M'},
		{"Read_length", 1, NULL, 't'}
	};
	while ((c = getopt_long(argc, argv, optString, longOpts, NULL)) >= 0)
	{
		switch (c)
		{
			case 'k': globalArgs.k = atoi(optarg); break;
			case 'n': globalArgs.n = atoi(optarg); break;
			case 'c': globalArgs.kmercov = atoi(optarg); break;
			case 'C': globalArgs.readcov = atoi(optarg); break;
			case 'R': globalArgs.readcov0 = atoi(optarg); break;
			case 'd': globalArgs.d1 = atoi(optarg); break;
			case 'D': globalArgs.d2 = atoi(optarg); break;
			case 'e': globalArgs.err = atoi(optarg); break;
			case 'm': globalArgs.lnth0 = atoi(optarg); break;
			case 'L': globalArgs.lnth1 = atoi(optarg); break;
			case 'o': globalArgs.output = optarg; break;
			case '1': globalArgs.contig_info = optarg; break;
			case '2': globalArgs.MUP_info = optarg; break;
			case 's': globalArgs.single = optarg; break;
			case 'l': globalArgs.left = optarg; break;
			case 'r': globalArgs.right = optarg; break;
			case 'p': globalArgs.paircov = atoi(optarg); break;
			case 'S': globalArgs.Min = atoi(optarg); break;
			case 'M': globalArgs.Max = atoi(optarg); break;
			case 't': globalArgs.readlen = atoi(optarg); break;
		}
	}
	int m = globalArgs.k - globalArgs.n;
	unsigned int pairkey = 0;
	clock_t start;
	clock_t end;
	ofstream out(globalArgs.output);
	start = clock();
	ReadTable V = Input_Reads(globalArgs.single, globalArgs.left, globalArgs.right, pairkey);
	end = clock();
	printf("time to read reads %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
	reads * q = V;
	if (V == 0)
	{
		printf("error input files!\n");
		return 0;
	}
	KmerTable T(1 << 2 * m);
	Build_Kmer(T, V, globalArgs.k, globalArgs.n);
	printf("start to build Kmer-DeBrujin!\n");
	start = clock();
	Del_Low_Kmer(T, globalArgs.kmercov, globalArgs.k, globalArgs.n);
	Build_DeBrujin(T, globalArgs.kmercov, globalArgs.k, globalArgs.n);	
//to build Debrujin Graph for all kmers
	end = clock();
	printf("time to Build_DeBrujin %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
	MUP_Graph S;
	start = clock();
    Form_MUP(T, S, globalArgs.lnth0);	//save all MUPs to MUP_graph
	end = clock();
	printf("time to build MUP %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
	start = clock();
	for (int i = 1; i != 6; ++i)
   	Form_MUP_Graph(T, V, S, globalArgs.k, globalArgs.n, i, globalArgs.readcov, globalArgs.readcov0, 1);
//1 is the tolerated err for distance computed by AM
	end = clock();
	printf("time to build MUP graph %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
//	Print_MUP(T, S, out);
	printf("start to deal complex region\n");
	start = clock();
	deal_complexregion(T, S, V, globalArgs.readcov0, globalArgs.d1);
	deal_complexregion(T, S, V, globalArgs.readcov0, globalArgs.d2);
	end = clock();
	printf("time to finish solve bubble %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
	Preprocess_LongMUP(S);
	LongMUP_Graph L;
	printf("start to build LongMUP Graph\n");
	start = clock();
	Form_LongMUP(L, S);
	Form_LongMUP_Graph(L, S);
	end = clock();
	printf("time to build LongMUP Graph %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
	ofstream out1("temp");
	start = clock();
	printf("time to print contig!\n");
	Print_Contig(T, V, S, L, pairkey, globalArgs.Max, globalArgs.Min,\
			globalArgs.readlen, globalArgs.readcov, globalArgs.paircov, out1);
	end = clock();
	printf("finish printing %lf\n", (double)(end - start) / CLOCKS_PER_SEC);

	return 0;
}
