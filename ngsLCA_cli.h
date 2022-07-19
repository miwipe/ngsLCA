#pragma once
#include <zlib.h>
#include <htslib/sam.h>
#include <cstring>
#include <vector>
#include <map>
#include <sys/stat.h>
#include <pthread.h>

struct cmp_str
{
   bool operator()(char const *a, char const *b) const
   {
      return std::strcmp(a, b) < 0;
   }
};

typedef struct{
  int up;
  std::vector<int> down;
  int dist2root;
}node;

typedef std::map<int,char *> int2char;
typedef std::map<int,int> int2int;
typedef std::map<int,node> int2node;
typedef std::map<char *,int,cmp_str> char2int;

typedef struct{
  //filenames
  char *htsfile;//bam,cram,sam
  char *nodesfile;
  char *namesfile;
  char *acc2taxfile;
  //hts strucutures
  samFile *hts;
  sam_hdr_t *header;
  //parameters for filtering reads
  double simscoreLow;
  double simscoreHigh;
  int editdistMin;
  int editdistMax;
  char *outnames;
  FILE *fp1;
  FILE *fp2;
  FILE *fp3;
  int minmapq;
  int discard; //or bitoperation with the flag of the read
  int minlength;
  char2int *charref2taxid;
}pars;

pars *get_pars(int argc,char **argv);
void print_pars(FILE *fp,pars *p);
void pars_free(pars *p);
int fexists(const char* str);
int fexists2(const char*str1,const char* str2);
int fexists3(const char*str1,const char* str2,const char *str3);
BGZF *getbgzf(const char*str1,const char *mode,int nthreads);
BGZF *getbgzf2(const char*str1,const char *str2,const char *mode,int nthreads);
BGZF *getbgzf3(const char*str1,const char *str2,const char *str3,const char *mode,int nthreads);

