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
  p->editdistMin=0;
  p->editdistMax=10;
  p->simscoreLow=0;
  p->simscoreHigh=1;
  p->fp1=p->fp2=p->fp3=NULL;
  p->outnames="outnames";
  p->minmapq=0;
  p->discard=516;//discard unmapped and read fail
  p->gz_sam = Z_NULL;
  p->minlength=-1;
  return p;
}

void pars_free(pars *p){
  fclose(p->fp1);
  //fclose(p->fp2);
  fclose(p->fp3);
  free(p);
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
    else if(!strcasecmp("-sam_noheader",key)) p->gz_sam=gzopen(val,"r");
    else if(!strcasecmp("-nodes",key)) p->nodesfile=strdup(val);
    else if(!strcasecmp("-acc2tax",key)) p->acc2taxfile=strdup(val);
    else if(!strcasecmp("-editdistMin",key)) p->editdistMin=atoi(val);
    else if(!strcasecmp("-editdistMax",key)) p->editdistMax=atoi(val);
    else if(!strcasecmp("-minmapq",key)) p->minmapq=atoi(val);
    else if(!strcasecmp("-minlength",key)) p->minlength=atoi(val);
    else if(!strcasecmp("-simscoreLow",key)) p->simscoreLow=atof(val);
    else if(!strcasecmp("-simscoreHigh",key)) p->simscoreHigh=atof(val);
    else if(!strcasecmp("-outnames",key)) p->outnames=strdup(val);
    else if(!strcasecmp("-out",key)) p->outnames=strdup(val);
    else if(!strcasecmp("-discard",key)) p->discard=atoi(val);
    else{
      fprintf(stderr,"\t Unknown parameter key:%s val:%s\n",key,val);
      free(p);
      return NULL;
    }
    
    ++argv;
  }
  if(p->gz_sam==Z_NULL){
    fprintf(stderr,"\t-> Will read header\n");
    p->hts = hts_open(p->htsfile,"r");
    p->header = sam_hdr_read(p->hts);
    assert(p->header);
    fprintf(stderr,"\t-> Done reading header\n");
  }
  char buf[1024];
  snprintf(buf,1024,"%s.lca",p->outnames);
  fprintf(stderr,"\t-> Will output lca results in file:\t\t\'%s\'\n",buf);
  p->fp1 = fopen(buf,"wb");

  snprintf(buf,1024,"%s.wlca",p->outnames);
  fprintf(stderr,"\t-> Will output lca weight in file:\t\t\'%s\'\n",buf);
  //  p->fp2 = fopen(buf,"wb");

  snprintf(buf,1024,"%s.log",p->outnames);
  fprintf(stderr,"\t-> Will output log info (problems) in file:\t\'%s\'\n",buf);
  p->fp3 = fopen(buf,"wb");

  return p;
}

void print_pars(FILE *fp,pars *p){
  fprintf(fp,"\t-> -bam  \t%s\n",p->htsfile);
  fprintf(fp,"\t-> -names\t%s\n",p->namesfile);
  fprintf(fp,"\t-> -nodes\t%s\n",p->nodesfile);
  fprintf(fp,"\t-> -acc2tax\t%s\n",p->acc2taxfile);
  fprintf(fp,"\t-> -simscoreLow\t%f\n",p->simscoreLow);
  fprintf(fp,"\t-> -simscoreHigh\t%f\n",p->simscoreHigh);
  fprintf(fp,"\t-> -editdistMin\t%d\n",p->editdistMin);
  fprintf(fp,"\t-> -editdistMax\t%d\n",p->editdistMax);
  fprintf(fp,"\t-> -outnames\t%s\n",p->outnames);
  fprintf(fp,"\t-> -minmapq\t%d\n",p->minmapq);

}


#ifdef __WITH_MAIN__
int main(int argc,char**argv){
  pars *p=get_pars(--argc,++argv);
  assert(p);
  
  print_pars(stdout,p);
}
#endif
