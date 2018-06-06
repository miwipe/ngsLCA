int mod_in[] =  {1649555 , 1401172 ,1582271, 374764, 242716,1793725 ,292451,38298,  63403  ,67357  ,\
		 163247  ,328623   ,356150 , 502130, 545877,996989  ,996990,1086724,1169024,1576640,\
		 1769757,1802981,1811974 ,1955118};
int mod_out[]=  {1333996 , 1333996 ,1582270,1914213,1917265,1915309 ,263865,2801 ,1916091,285450,1762941,1916091,157727,1932322,376133,1762939 ,1762946,430531, 1169025,1247960,1769758,1708715,1708715	,1925741};

#include <cassert>
#include <cstdio>
#include <zlib.h>
#include <cstring>
#include <map>
#include <cstdlib>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <vector>
#include <pthread.h>
#include <algorithm>
#include <errno.h>
#include "ngsLCA.h"
#include "ngsLCA_cli.h"
#include "ngsLCA_format.h"
int2int i2i; //refid to taxid
int2int specWeight;// Number of reads that map uniquely to a species.


void mod_db(int *in,int *out,int2int &parent, int2char &rank,int2char &name_map){
  for(int i=0;i<24;i++){
    assert(parent.count(out[i])==1);
    parent[in[i]] = parent[out[i]];
    rank[in[i]] = rank[out[i]];
    name_map[in[i]] = strdup("satan");
  }

}
int2int errmap;

//usefull little function to split
char *strpop(char **str,char split){
  char *tok=*str;
  while(**str){
    if(**str!=split)
      (*str)++;
    else{
      **str='\0'; (*str)++;
      break;
    }
  }
  return tok;
}
//usefull little function to remove tab and newlines
void strip(char *line){
  int at=0;
  //  fprintf(stderr,"%s\n",line);
  for(int i=0;i<strlen(line);i++)
    if(line[i]=='\t'||line[i]=='\n')
      continue;
    else
      line[at++]=line[i];
  //  fprintf(stderr,"at:%d line:%p\n",at,line);
  line[at]='\0';
  //fprintf(stderr,"%s\n",line);
}



 
int2int bamRefId2tax(const char *fname,bam_hdr_t *hdr ){
  fprintf(stderr,"\t-> Number of SQ tags:%d \n",hdr->n_targets);
  int2int am;
  gzFile gz= Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problems opening file: \'%s\'\n",fname);
    exit(0);
  }
  char buf[4096];
  int at=0;
  FILE *fp = NULL;
  if(0){
    fp=fopen("dmp.dmp.deleteme","w");
    assert(fp);
  }
  while(gzgets(gz,buf,4096)){
    if(!((at++ %100000 ) ))
      if(isatty(fileno(stderr)))
	fprintf(stderr,"\r\t-> At linenr: %d in \'%s\'      ",at,fname);
    strtok(buf,"\t\n ");
    char *key =strtok(NULL,"\t\n ");
    int val = atoi(strtok(NULL,"\t\n "));

    //check if the key exists in the bamheader, if not then skip this taxid
    int valinbam = bam_name2id(hdr,key);
    if(valinbam==-1)
      continue;
    
    if(am.find(valinbam)!=am.end())
      fprintf(stderr,"\t-> Duplicate entries found \'%s\'\n",key);
    else{
      am[valinbam] = val;
      if(fp)
	fprintf(fp,"val\t%s\t%d\tjunk\n",key,val);
    }

  }
  fprintf(stderr,"\n");
  fprintf(stderr,"\t-> [%s] Number of entries to use from accesion to taxid: %lu\n",fname,am.size());
  if(fp)
    fclose(fp);
  return am;
}

int nodes2root(int taxa,int2int &parent){
  int dist =0;
  while(taxa!=1){
    int2int::iterator it = parent.find(taxa);
    taxa = it->second;
    dist++;
  }
  return dist;
}

int2int dist2root;


int satan(int taxid,int2int &parant){
  int2int::iterator it=dist2root.find(taxid);
  
  return 0;
}



int do_lca(std::vector<int> &taxids,int2int &parent){
  //  fprintf(stderr,"\t-> [%s] with number of taxids: %lu\n",__func__,taxids.size());
  assert(taxids.size()>0);
  if(taxids.size()==1){
    int taxa=taxids[0];
    if(parent.count(taxa)==0){
      fprintf(stderr,"\t-> Problem finding taxaid: %d will skip\n",taxa);
      taxids.clear();
      return -1;
    }

    taxids.clear();
    return taxa;
  }

  int2int counter;
  for(int i=0;i<taxids.size();i++){
    int taxa = taxids[i];
    while(1){
      //      fprintf(stderr,"taxa:%d\n",taxa);
      int2int::iterator it=counter.find(taxa);
      if(it==counter.end()){
	//	fprintf(stderr,"taxa: %d is new will plugin\n",taxa);
	counter[taxa]=1;
      }else
	it->second = it->second+1;
      it = parent.find(taxa);

      if(it==parent.end()){
	int2int::iterator it=errmap.find(taxa);
	if(it==errmap.end()){
	  fprintf(stderr,"\t-> Problem finding parent of :%d\n",taxa);
	  
	  errmap[taxa] = 1;
	}else
	  it->second = it->second +1;
	taxids.clear();
	return -1;
      }
      
      if(taxa==it->second)//<- catch root
	break;
      taxa=it->second;
    }

  }
  //  fprintf(stderr,"counter.size():%lu\n",counter.size());
  //now counts contain how many time a node is traversed to the root
  int2int dist2root;
  for(int2int::iterator it=counter.begin();it!=counter.end();it++)
    if(it->second==taxids.size())
      dist2root[nodes2root(it->first,parent)] = it->first;
  for(int2int::iterator it=dist2root.begin();0&&it!=dist2root.end();it++)
    fprintf(stderr,"%d\t->%d\n",it->first,it->second);
  taxids.clear();
  if(!dist2root.empty())
    return (--dist2root.end())->second;
  else{
    fprintf(stderr,"\t-> Happens\n");
    return -1;
  }
}

void print_chain1(FILE *fp,int taxa,int2char &rank,int2char &name_map){
  int2char::iterator it1=name_map.find(taxa);
  int2char::iterator it2=rank.find(taxa);
  if(it1==name_map.end()){
    fprintf(stderr,"\t-> Problem finding taxaid:%d\n",taxa);
  }
  assert(it2!=rank.end());
  if(it1==name_map.end()||it2==rank.end()){
    fprintf(stderr,"taxa: %d %s doesnt exists will exit\n",taxa,it1->second);
    exit(0);
  }
  fprintf(fp,"\t%d:%s:%s",taxa,it1->second,it2->second);
  
}


void print_chain(FILE *fp,int taxa,int2int &parent,int2char &rank,int2char &name_map){

    while(1){
      print_chain1(fp,taxa,rank,name_map);
      int2int::iterator it = parent.find(taxa);
      assert(it!=parent.end());
      if(taxa==it->second)//<- catch root
	break;
      taxa=it->second;
    }
    fprintf(fp,"\n");
}

int isuniq(std::vector<int> &vec){
  int ret = 1;
  for(int i=1;i<vec.size();i++)
    if(vec[0]!=vec[i])
      return 0;
  return ret;
}

int get_species1(int taxa,int2int &parent, int2char &rank){
  //  fprintf(stderr,"\t-> taxa:%d\n",taxa);
  int2char::iterator it;
  while(1){
    it=rank.find(taxa);
    if(it==rank.end())
      return -1;
    char *val = it->second;
    if(val==NULL)
      fprintf(stderr,"taxa:%d\n",taxa);
    assert(val);
    if(!strcmp(val,"species"))
      return taxa;
    int next = parent[taxa];
    if(next==taxa)
      break;
    taxa=next;
  }
  return taxa;
  
}

int2int get_species(int2int &i2i,int2int &parent, int2char &rank,int2char &names,FILE *fp){
  int2int ret;
  for(int2int::iterator it=i2i.begin();it!=i2i.end();it++){
    //fprintf(stderr,"%d\t%d\n",it->first,it->second);
    int asdf= get_species1(it->second,parent,rank);
    if(asdf==-1){
      fprintf(fp,"\t-> Removing pair(%d,%d) accnumber:%s since doesnt exists in node list\n",it->first,it->second,names[it->second]);
      i2i.erase(it);
    }else
      ret[it->second] = asdf;

  }
  return ret;
}

char *make_seq(bam1_t *aln){
  int len = aln->core.l_qseq;
  char *qseq = new char[len+1];
  uint8_t *q = bam_get_seq(aln);
  for(int i=0; i< len ; i++)
    qseq[i] = seq_nt16_str[bam_seqi(q,i)];
  qseq[len]='\0';
  //  fprintf(stderr,"seq:%s\n",qseq);
  //exit(0);
  return qseq;
}

int printonce=1;

std::vector<int> purge(std::vector<int> &taxids,std::vector<int> &editdist){
  if(printonce==1)
    fprintf(stderr,"\t-> purging taxids oldsize:%lu\n",taxids.size());
  assert(taxids.size()==editdist.size());
  std::vector<int> tmpnewvec;
  int mylow = *std::min_element(editdist.begin(),editdist.end());
  for(int i=0;i<taxids.size();i++)
    if(editdist[i]<=mylow)
      tmpnewvec.push_back(taxids[i]);
  if(printonce--==1)
    fprintf(stderr,"\t-> purging taxids newsize:%lu this info is only printed once\n",tmpnewvec.size());
  return tmpnewvec;
}


void hts(FILE *fp,samFile *fp_in,int2int &i2i,int2int& parent,bam_hdr_t *hdr,int2char &rank, int2char &name_map,FILE *log,int minmapq,int discard,int editMin, int editMax, double scoreLow,double scoreHigh){
  fprintf(stderr,"\t-> editMin:%d editmMax:%d scoreLow:%f scoreHigh:%f\n",editMin,editMax,scoreLow,scoreHigh);
  assert(fp_in!=NULL);
  bam1_t *aln = bam_init1(); //initialize an alignment
  int comp ;

  char *last=NULL;
  char *seq =NULL;
  std::vector<int> taxids;
  std::vector<int> specs;
  std::vector<int> editdist;
  int lca;
  int2int closest_species;
  int skip=0;
  int inc=0;
  while(sam_read1(fp_in,hdr,aln) > 0) {
    char *qname = bam_get_qname(aln);
    int chr = aln->core.tid ; //contig name (chromosome)
    //    fprintf(stderr,"%d %d\n",aln->core.qual,minmapq);
    if(aln->core.qual<minmapq){
      fprintf(stderr,"discarding due to low mapq");
      continue;
    }
    //    fprintf(stderr,"Discard: %d coreflag:%d OR %d\n",discard,aln->core.flag,aln->core.flag&discard);
    if(discard>0&&(aln->core.flag&discard)){
      fprintf(stderr,"Discardding due to core flag\n");
      continue;
    }
    if(last==NULL){
      last=strdup(qname);
      seq=make_seq(aln);
    }
    //change of ref
    if(strcmp(last,qname)!=0) {
      if(taxids.size()>0&&skip==0){
	//	fprintf(stderr,"length of taxids:%lu and other:%lu minedit:%d\n",taxids.size(),editdist.size(),*std::min_element(editdist.begin(),editdist.end()));
	
	int size=taxids.size();
	if(editMin==-1&&editMax==-1)
	  taxids = purge(taxids,editdist);

	lca=do_lca(taxids,parent);
	if(lca!=-1){
	  fprintf(fp,"%s:%s:%lu:%d",last,seq,strlen(seq),size);//fprintf(stderr,"size:%d\n");
	  //	  fprintf(stderr,"adfsadfsafafdad: %d size\n",size);
	  print_chain(fp,lca,parent,rank,name_map);
	  int varisunique = isuniq(specs);
	  //fprintf(stderr,"varisunquieu:%d spec.size():%lu\n",varisunique,specs.size());
	  if(varisunique){
	    int2int::iterator it=specWeight.find(specs[0]);
	    //fprintf(stderr,"specs: %d specs.size:%lu wiehg.szei():%lu\n",specs[0],specs.size(),specWeight.size());
	    if(it==specWeight.end())
	      specWeight[specs[0]] = 1;//specs.size();
	    else
	      it->second = it->second +1;
	    
	    //fprintf(stderr,"specs.size:%lu wiehg.szei():%lu\n",specs.size(),specWeight.size());
	  }


	}
      }
      skip=0;
      specs.clear();
      editdist.clear();
      free(last);
      delete [] seq;
      last=strdup(qname);
      seq=make_seq(aln);
    }
    
    
    //filter by nm
    uint8_t *nm = bam_aux_get(aln,"NM");
    int thiseditdist;
    if(nm!=NULL){
      thiseditdist = (int) bam_aux2i(nm);
      //      fprintf(stderr,"[%d] nm:%d\t",inc++,val);
      if(editMin!=-1&&thiseditdist<editMin){
	skip=1;
	//fprintf(stderr,"skipped1\n");
	continue;
      }else if(editMax!=-1&&thiseditdist>editMax){
	//fprintf(stderr,"continued1\n");
	continue;
      }
      double seqlen=aln->core.l_qseq;
      double myscore=1.0-(((double) thiseditdist)/seqlen);
      //      fprintf(stderr," score:%f\t",myscore);
      if(myscore>scoreHigh){
	//fprintf(stderr,"skipped2\n");
	skip=1;
	continue;
      }
      else if(myscore<scoreLow){	
	//	fprintf(stderr,"continued2\n");
	continue;
      }
    }
    int2int::iterator it = i2i.find(chr);
    //See if cloests speciest exists and plug into closests species
    int dingdong=-1;
    if(it!=i2i.end()) {
      int2int::iterator it2=closest_species.find(chr);
      if(it2!=closest_species.end())
	dingdong=it->second;
      else{
	dingdong=get_species1(it->second,parent,rank);
	if(dingdong!=-1)
	  closest_species[it->second]=dingdong;
      }
      //fprintf(stderr,"\t-> closests size:%lu\n",closest_species.size());
    }

    if(it==i2i.end()||dingdong==-1)
      fprintf(log,"\t-> problem finding chrid:%d chrname:%s\n",chr,hdr->target_name[chr]);
    else{
            
      taxids.push_back(it->second);
      specs.push_back(dingdong);
      editdist.push_back(thiseditdist);
      //fprintf(stderr,"it-.second:%d specs:%d thiseditdist:%d\n",it->second,dingdong,thiseditdist);
      //      fprintf(stderr,"EDIT\t%d\n",thiseditdist);
    }
  }
  if(taxids.size()>0&&skip==0){
    int size=taxids.size();
    if(lca!=-1){
      if(editMin==-1&&editMax==-1)
	taxids = purge(taxids,editdist);
      
      lca=do_lca(taxids,parent);
      if(lca!=-1){
	fprintf(fp,"%s:%s:%d",last,seq,size);fflush(stdout);
	print_chain(fp,lca,parent,rank,name_map);
	if(isuniq(specs)){
	  int2int::iterator it=specWeight.find(specs[0]);
	  if(it==specWeight.end())
	    specWeight[specs[0]] = specs.size();
	  else
	    it->second = it->second +specs.size();
	  
	}
      }
    }
  }
  specs.clear();
  editdist.clear();
  bam_destroy1(aln);
  sam_close(fp_in);
  
  return ;//0;
}

int2char parse_names(const char *fname){

  gzFile gz= Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problems opening file: \'%s\'\n",fname);
    exit(0);
  }
  int2char name_map;
  char buf[4096];
  int at=0;
  char **toks = new char*[5];
  
  while(gzgets(gz,buf,4096)){
    strip(buf);//fprintf(stderr,"buf:%s\n",buf);
    char *saveptr = buf;
    toks[0]=strpop(&saveptr,'|');
    toks[1]= strpop(&saveptr,'|');
    toks[2]= strpop(&saveptr,'|');
    toks[3]= strpop(&saveptr,'|');
    for(int i=0;0&&i<4;i++)
      fprintf(stderr,"%d):\'%s\'\n",i,toks[i]);

    int key=atoi(toks[0]);
    //    fprintf(stderr,"key:%d\n",key);
    if(toks[3]&&strcmp(toks[3],"scientific name")==0){
      int2char::iterator it=name_map.find(key);
      
      if(it!=name_map.end())
	fprintf(stderr,"\t->[%s] duplicate name(column1): %s\n",fname,toks[0]);
      else
	name_map[key]=strdup(toks[1]);

    }
    if(0&&at++>10)
      break;
  }
  //  int2char::iterator it = name_map.find(61564);  assert(it!=name_map.end());
  fprintf(stderr,"\t-> [%s] Number of unique names (column1): %lu with third column 'scientific name'\n",fname,name_map.size());
  return name_map;
}


void parse_nodes(const char *fname,int2char &rank,int2int &parent){

  gzFile gz= Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problems opening file: \'%s\'\n",fname);
    exit(0);
  }
  char buf[4096];
  int at=0;
  char **toks = new char*[5];

  while(gzgets(gz,buf,4096)){
    strip(buf);//fprintf(stderr,"buf:%s\n",buf);
    char *saveptr = buf;
    toks[0]= strpop(&saveptr,'|');
    toks[1]= strpop(&saveptr,'|');
    toks[2]= strpop(&saveptr,'|');
    for(int i=0;0&&i<3;i++)
      fprintf(stderr,"%d):\'%s\'\n",i,toks[i]);

    int2int::iterator it = parent.find(atoi(toks[0]));
    if(it!=parent.end())
      fprintf(stderr,"\t->[%s] duplicate name(column0): %s\n",fname,toks[0]);
    else{
      int key=atoi(toks[0]);
      int val= atoi(toks[1]);
      parent[key]=val;
      rank[key]=strdup(toks[2]);
    }
  }
  fprintf(stderr,"\t-> Number of unique names (column1): %lu from file: %s\n",rank.size(),fname);

}

void print_ref_rank_species(bam_hdr_t *h,int2int &i2i,int2char &names,int2char &rank){
  for(int i=0;i<h->n_targets;i++){
    fprintf(stdout,"%d\t%s\t%s\n",i,names[i2i[i]],rank[i2i[i]]);
    
  }

}

int calc_valens(int2int &i2i, int2int &parent){
  int2int ret;
  for(int2int::iterator it=i2i.begin();it!=i2i.end();it++){
    int2int::iterator pit=parent.find(it->second);
    if(pit==parent.end()||pit->second==1)
      continue;
    //    fprintf(stdout,"%d %d\n",pit->first,pit->second);
    int2int::iterator rit=ret.find(pit->second);
    if(rit==ret.end())
      ret[pit->second]=1;
    else
      rit->second = rit->second +1;
  }

  for(int2int::iterator it=ret.begin();it!=ret.end();it++)
    fprintf(stdout,"%d\t%d\n",it->first,it->second);
  return 0;
}



int calc_dist2root(int2int &i2i, int2int &parent){
  int2int ret;
  for(int2int::iterator it=i2i.begin();it!=i2i.end();it++){
    int2int::iterator pit=parent.find(it->second);
    if(pit==parent.end()||pit->second==1)
      continue;
    //    fprintf(stdout,"%d %d\n",pit->first,pit->second);
    int2int::iterator rit=ret.find(pit->second);
    if(rit==ret.end())
      ret[pit->second]=1;
    else
      rit->second = rit->second +1;
  }

  for(int2int::iterator it=ret.begin();it!=ret.end();it++)
    fprintf(stdout,"%d\t%d\n",it->first,it->second);
  return 0;
}
#if 0
int2node makeNodes(int2int &parent){
  int2node ret;
  for(int2int::iterator it=parent.begin();it!=parent.end();it++){
    int2node::iterator it2=ret.find(it->second);
    if(it2==ret.end()){
      node nd;
      nd.up=it->second;
      ret[it->first] = nd;
    }
    
  }

  
  return ret;
}
#endif




int main(int argc, char **argv){
  if(argc==1){
    fprintf(stderr,"\t-> ngsLCA -names -nodes -acc2tax [-editdist[min/max] -simscore[low/high] -minmapq -discard] -bam \n");
    return 0;
  }
#if 0
  int2int ww = get_weight(argv[1]);
  return 0;
#endif
  time_t t2=time(NULL);
  if(argc>1&&!strcasecmp(argv[1],"format")){
    return ngsLCA_format(--argc,++argv);
  }
  pars *p=get_pars(--argc,++argv);
  print_pars(stderr,p);
  
  //map of bamref ->taxid

  int2int i2i= bamRefId2tax(p->acc2taxfile,p->header);
 
  //map of taxid -> taxid
  int2int parent;
  //map of taxid -> rank
  int2char rank;
  //map of taxid -> name
  int2char name_map = parse_names(p->namesfile);
  parse_nodes(p->nodesfile,rank,parent);
  //  calc_valens(i2i,parent);
  if(0){
    print_ref_rank_species(p->header,i2i,name_map,rank);
    return 0;
  }

  //closes species (direction of root) for a taxid
  //  int2int closest_species=get_species(i2i,parent,rank,name_map,p->fp3);
  //  fprintf(stderr,"\t-> Number of items in closest_species map:%lu\n",closest_species.size());

  
  fprintf(stderr,"\t-> Will add some fixes of the ncbi database due to merged names\n");
  mod_db(mod_in,mod_out,parent,rank,name_map);

  hts(p->fp1,p->hts,i2i,parent,p->header,rank,name_map,p->fp3,p->minmapq,p->discard,p->editdistMin,p->editdistMax,p->simscoreLow,p->simscoreHigh);  
  fprintf(stderr,"\t-> Number of species with reads that map uniquely: %lu\n",specWeight.size());
  
  for(int2int::iterator it=errmap.begin();it!=errmap.end();it++)
    fprintf(p->fp3,"err\t%d\t%d\n",it->first,it->second);

  for(int2int::iterator it=specWeight.begin();it!=specWeight.end();it++)
    fprintf(p->fp2,"%d\t%s\t%d\n",it->first,name_map[it->first],it->second);
  pars_free(p);
  fprintf(stderr, "\t-> [ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  return 0;
}
