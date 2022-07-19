#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <cassert>
#include <htslib/bgzf.h>
#include "ngsLCA_cli.h"
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
pars *pars_init(){
  pars *p =(pars*) calloc(1,sizeof(pars));
  p->htsfile = strdup("CHL_155_12485.sort.bam");
  p->acc2taxfile=strdup("nucl_gb.accession2taxid.gz");
  p->namesfile = strdup("names.dmp.gz");
  p->nodesfile= strdup("nodes.dmp.gz");
  p->hts=NULL;
  p->header=NULL;
  p->editdistMin=0;
  p->editdistMax=10;
  p->simscoreLow=0;
  p->simscoreHigh=1;
  p->fp1=p->fp2=p->fp3=NULL;
  p->outnames=strdup("outnames");
  p->minmapq=0;
  p->discard=516;//discard unmapped and read fail
  p->minlength=-1;
  p->charref2taxid = NULL;
  return p;
}

void pars_free(pars *p){
  fclose(p->fp1);
  //fclose(p->fp2);
  fclose(p->fp3);
  if(p->htsfile)
    free(p->htsfile);
  if(p->acc2taxfile)
    free(p->acc2taxfile);
  if(p->nodesfile)
    free(p->nodesfile);
  if(p->header)
    sam_hdr_destroy(p->header);
  if(p->hts)
    sam_close(p->hts);
  free(p);
}

int checkIfSorted(char *str){

   //check if proper header exists
   if(strncmp(str,"@HD",3)!=0){
     fprintf(stderr,"\t-> We require a proper header starting with @HD for metadamage\n");
     fprintf(stderr,"\t-> We observed: \'%.10s\' will exit\n",str);
     return 1;
   }
   //check if SO:coordinate exists
   char *so = strstr(str,"SO:queryname");
   if(so==NULL){
     fprintf(stderr,"\t-> ERROR: We require files to be sorted by readname, will exit\n");
     return 2;
   }
   if(strchr(str,'\n')<so){
     fprintf(stderr,"\t-> We require a SO:queryname tag in the first line of header\n");
     return 3;
   }
   return 0;
 }

int fexists(const char* str){///@param str Filename given as a string.
  fprintf(stderr,"\t-> Checking if exits: \'%s\'\n",str);
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

int fexists2(const char*str1,const char* str2){
  unsigned tmp_l = strlen(str1)+strlen(str2)+5;
  char tmp[tmp_l];
  snprintf(tmp,tmp_l,"%s%s",str1,str2);
  return fexists(tmp);
}


int fexists3(const char*str1,const char* str2,const char *str3){
  unsigned tmp_l = strlen(str1)+strlen(str2)+strlen(str3)+5;
  char tmp[tmp_l];
  snprintf(tmp,tmp_l,"%s%s%s",str1,str2,str3);
  return fexists(tmp);
}



void *read_header_thread(void *ptr){
  time_t t=time(NULL);
  pars *p = (pars *) ptr;
  fprintf(stderr,"\t-> [thread1] Will read header\n");
  p->hts = hts_open(p->htsfile,"r");
  p->header = sam_hdr_read(p->hts);
  assert(p->header);
  fprintf(stderr,"\t-> [thread1] Done reading header: %.2f sec, header contains: %d \n",(float)(time(NULL) - t),p->header->n_targets);
  checkIfSorted(p->header->text);
  pthread_mutex_unlock(&mutex1);
  pthread_exit(0);
}


char2int * ass2bin(const char *fname,int redo){
  const char *CONSTNAME = "delmeme.bin";
  gzFile FP=Z_NULL;
  char2int *cm = new char2int;
  //redo=1;
  //load binary representation
  if(redo==0&&fexists(CONSTNAME)){
    time_t t=time(NULL);
    fprintf(stderr,"\t-> [thread2] reading binary represenation\n");
    FP = gzopen(CONSTNAME,"rb");
    int key_l;
    
    while(sizeof(int)==gzread(FP,&key_l,sizeof(int))) {
      char *key =(char*) calloc(key_l+1,sizeof(char));
      assert(key_l=gzread(FP,key,key_l));
      int val;
      assert(sizeof(int)==gzread(FP,&val,sizeof(int)));
      if(cm->find(key)!=cm->end()){
	fprintf(stderr,"\t-> Duplicate entries found \'%s\'\n",key);
      }else
	(*cm)[key] = val;
    }
    fprintf(stderr,"\t-> [thread2] Done reading binary representation: %.2f sec\n",(float)(time(NULL) - t));

  }else{
    FP=gzopen(CONSTNAME,"wb");
  
    gzFile gz= Z_NULL;
    gz=gzopen(fname,"rb");
    if(gz==Z_NULL){
      fprintf(stderr,"\t-> Problems opening file: \'%s\'\n",fname);
      exit(0);
    }
    char buf[4096];
    int at=0;
    char buf2[4096];
    extern int SIG_COND;
    gzgets(gz,buf,4096);//skip header
    while(SIG_COND&&gzgets(gz,buf,4096)){
      if(!((at++ %100000 ) ))
	if(isatty(fileno(stderr)))
	  fprintf(stderr,"\r\t-> At linenr: %d in \'%s\'      ",at,fname);
      strcpy(buf2,buf);
      strtok(buf,"\t\n ");
      char *key =strtok(NULL,"\t\n ");
      int val = atoi(strtok(NULL,"\t\n "));
      if(FP){
	int key_l = strlen(key);
	gzwrite(FP,&key_l,sizeof(int));
	gzwrite(FP,key,key_l);
	gzwrite(FP,&val,sizeof(int));
	fprintf(stderr,"key_l: %d key:%s val:%d\n",key_l,key,val);
      }
      if(cm->find(key)!=cm->end())
	fprintf(stderr,"\t-> Duplicate entries found \'%s\'\n",key);
      (*cm)[strdup(key)]=val;
    }
  }

  if(FP)
    gzclose(FP);

  fprintf(stderr,"\n");
  fprintf(stderr,"\t-> [%s] Number of entries to use from accesion to taxid: %lu\n",fname,cm->size());
  return cm;
}



void *read_ass2taxid_thread(void *ptr){
  pars *p = (pars *) ptr;
  p->charref2taxid = ass2bin(p->acc2taxfile,0);
  pthread_mutex_unlock(&mutex2);
  pthread_exit(0);
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
    
    if(!strcasecmp("-bam",key)){
      if(p->htsfile)
	free(p->htsfile);
      p->htsfile=strdup(val);
    }
    else if(!strcasecmp("-names",key)){
      if(p->namesfile)
	free(p->namesfile);
      p->namesfile=strdup(val);
    }
    else if(!strcasecmp("-nodes",key)){
      if(p->nodesfile)
	free(p->nodesfile);
      p->nodesfile=strdup(val);
    }
    else if(!strcasecmp("-acc2tax",key)){
      if(p->acc2taxfile)
	free(p->acc2taxfile);
      p->acc2taxfile=strdup(val);
    }
    else if(!strcasecmp("-editdistMin",key)) p->editdistMin=atoi(val);
    else if(!strcasecmp("-editdistMax",key)) p->editdistMax=atoi(val);
    else if(!strcasecmp("-minmapq",key)) p->minmapq=atoi(val);
    else if(!strcasecmp("-minlength",key)) p->minlength=atoi(val);
    else if(!strcasecmp("-simscoreLow",key)) p->simscoreLow=atof(val);
    else if(!strcasecmp("-simscoreHigh",key)) p->simscoreHigh=atof(val);
    else if(!strcasecmp("-outnames",key)){
      if(p->outnames)
	free(p->outnames);
      p->outnames=strdup(val);
    }
    else if(!strcasecmp("-out",key)){
      if(p->outnames)
	free(p->outnames);
      p->outnames=strdup(val);
    }
    else if(!strcasecmp("-discard",key)) p->discard=atoi(val);
    else{
      fprintf(stderr,"\t Unknown parameter key:%s val:%s\n",key,val);
      free(p);
      return NULL;
    }
    
    ++argv;
  }
  pthread_t thread1, thread2;
  pthread_mutex_lock(&mutex1);
  //  pthread_mutex_lock(&mutex2);
  assert(pthread_create( &thread1, NULL, read_header_thread, (void*) p)==0);
  //  assert(pthread_create( &thread2, NULL, read_ass2taxid_thread, (void*) p)==0);
  
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

  pthread_mutex_lock(&mutex1);
  pthread_mutex_lock(&mutex2);
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

BGZF *getbgzf(const char*str1,const char *mode,int nthreads){
  BGZF *fp = NULL;
  fp = bgzf_open(str1,mode);
  fprintf(stderr,"\t-> opening file: \'%s\' mode: \'%s\'\n",str1,mode);
  assert(fp!=NULL);
  if(nthreads>1){
    fprintf(stderr,"\t-> Setting threads to: %d \n",nthreads);
    bgzf_mt(fp,nthreads,64);
  }
  return fp;
}

BGZF *getbgzf2(const char*str1,const char *str2,const char *mode,int nthreads){
  unsigned tmp_l = strlen(str1)+strlen(str2)+5;
  char tmp[tmp_l];
  snprintf(tmp,tmp_l,"%s%s",str1,str2);
  return  getbgzf(tmp,mode,nthreads);
}

BGZF *getbgzf3(const char*str1,const char *str2,const char *str3,const char *mode,int nthreads){
  unsigned tmp_l = strlen(str1)+strlen(str2)+strlen(str3)+5;
  char tmp[tmp_l];
  snprintf(tmp,tmp_l,"%s%s%s",str1,str2,str3);
  return  getbgzf(tmp,mode,nthreads);
}

#ifdef __WITH_MAIN__
int main(int argc,char**argv){
  pars *p=get_pars(--argc,++argv);
  assert(p);
  
  print_pars(stdout,p);
}
#endif
