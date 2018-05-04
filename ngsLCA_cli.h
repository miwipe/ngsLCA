#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <htslib/sam.h>
  
typedef struct{
  //filenames
  char *htsfile;//bam,cram,sam
  char *nodesfile;
  char *namesfile;
  char *acc2taxfile;
  //hts strucutures
  samFile *hts;
  bam_hdr_t *header;
  //parameters for filtering reads
  double simscore;
  int editdist;
  char *outnames;
  FILE *fp1;
  FILE *fp2;
  FILE *fp3;
  int minmapq;
  int discard; //or bitoperation with the flag of the read
}pars;


pars *get_pars(int argc,char **argv);
void print_pars(FILE *fp,pars *p);
void pars_free(pars *p);
  
#ifdef __cplusplus
}
#endif
  
