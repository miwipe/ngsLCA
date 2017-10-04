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
}pars;


pars *get_pars(int argc,char **argv);
void print_pars(FILE *fp,pars *p);

  
#ifdef __cplusplus
}
#endif
  
