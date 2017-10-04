#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "ngsLCA_cli.h"

pars *pars_init(){
  pars *p =(pars*) calloc(1,sizeof(pars));
  p->htsfile = "CHL_155_12485.sort.bam";
  p->acc2taxfile="nucl_gb.accession2taxid.gz";
  p->namesfile = "names.dmp.gz";
  p->nodesfile= "nodes.dmp.gz";
  p->hts=NULL;
  p->header=NULL;
  p->editdist=0;
  p->simscore=1;
  return p;
}


pars *get_pars(int argc,char **argv){
  pars *p = pars_init();
  if(argc % 2){
    fprintf(stderr,"\t-> Must supply arguments in the form -pattern value\n");
    free(p);
    return NULL;
  }
  while(*argv){
    char *key=*argv;
    char *val=*(++argv);
    
    if(!strcasecmp("-bam",key)) p->htsfile=strdup(val);
    else if(!strcasecmp("-names",key)) p->namesfile=strdup(val);
    else if(!strcasecmp("-nodes",key)) p->nodesfile=strdup(val);
    else if(!strcasecmp("-acc2tax",key)) p->acc2taxfile=strdup(val);
    else if(!strcasecmp("-editdist",key)) p->editdist=atoi(val);
    else if(!strcasecmp("-simscore",key)) p->simscore=atof(val);
    else{
      fprintf(stderr,"\t Unknown parameter key:%s val:%s\n",key,val);
      free(p);
      return NULL;
    }
    
    ++argv;
  }
  p->hts = hts_open(p->htsfile,"r");
  p->header = sam_hdr_read(p->hts);
  assert(p->header);
  return p;
}

void print_pars(FILE *fp,pars *p){
  fprintf(fp,"\t-> -bam  \t%s\n",p->htsfile);
  fprintf(fp,"\t-> -names\t%s\n",p->namesfile);
  fprintf(fp,"\t-> -nodes\t%s\n",p->nodesfile);
  fprintf(fp,"\t-> -acc2tax\t%s\n",p->acc2taxfile);
  fprintf(fp,"\t-> -simscore\t%f\n",p->simscore);
  fprintf(fp,"\t-> -editdist\t%d\n",p->editdist);

}


#ifdef __WITH_MAIN__
int main(int argc,char**argv){
  pars *p=get_pars(--argc,++argv);
  assert(p);
  
  print_pars(stdout,p);
}
#endif
