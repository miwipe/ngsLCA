#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct{
  char *htsfile;//bam,cram,sam
  char *nodesfile;
  char *namesfile;
  char *acc2taxfile;
}pars;


pars *get_pars(int argc,char **argv);
void print_pars(FILE *fp,pars *p);

  
#ifdef __cplusplus
}
#endif
  
