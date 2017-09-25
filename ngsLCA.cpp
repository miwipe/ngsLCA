#include <cassert>
#include <cstdio>
#include <zlib.h>
#include <cstring>
#include <map>
#include <cstdlib>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <vector>

struct cmp_str
{
   bool operator()(char const *a, char const *b)
   {
      return std::strcmp(a, b) < 0;
   }
};

typedef std::map<char *, int, cmp_str> char2int;
typedef std::map<int,char *> int2char;
typedef std::map<int,int> int2int;


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




int2int ref2tax(const char *fname,bam_hdr_t *hdr ){
  char2int revmap;
  for(int i=0;i<hdr->n_targets;i++)
    revmap[hdr->target_name[i]] = i;
  fprintf(stderr,"\t-> Number of SQ tags:%d\n",revmap.size());


  
  int2int am;
  gzFile gz= Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problems opening file: \'%s\'\n",fname);
    exit(0);
  }
  char buf[4096];
  int at=0;
  while(gzgets(gz,buf,4096)){
    if(!((at++ %100000 ) ))
       fprintf(stderr,"\r\t-> At linenr: %d in \'%s\'      ",at,fname);
    strtok(buf,"\t\n ");
    char *key =strtok(NULL,"\t\n ");
    int val = atoi(strtok(NULL,"\t\n "));

    //check if the key exists in the bamheader, if not then skip this taxid
    char2int::iterator it=revmap.find(key);
    if(it==revmap.end())
      continue;

    if(am.find(it->second)!=am.end())
      fprintf(stderr,"\t-> Duplicate entries found \'%s\'\n",key);
    else
      am[it->second] = val;

  }
  fprintf(stderr,"\n");
  fprintf(stderr,"\t-> [%s] Number of entries to use from accesion to taxid: %lu\n",fname,am.size());
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



int do_lca(std::vector<int> &taxids,int2int &parent){
  //  fprintf(stderr,"\t-> [%s] with number of taxids: %lu\n",__func__,taxids.size());
  if(taxids.size()<2){
    taxids.clear();
    return 0;
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
      assert(it!=parent.end());
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
}


void hts(const char *fname,int2int &i2i,int2int& parent,bam_hdr_t *hdr){
  samFile *fp_in = hts_open(fname,"r"); //open bam file
  fprintf(stderr,"fp_in:%p\n",fp_in);
  bam1_t *aln = bam_init1(); //initialize an alignment
  bam_hdr_t *bamHdr = sam_hdr_read(fp_in)  ;
  fprintf(stderr,"bamhdfadf:%p\n",bamHdr);
  int comp ;

  char *last=NULL;
  std::vector<int> taxids;
  int lca;
  while(sam_read1(fp_in,bamHdr,aln) > 0){
    char *qname = bam_get_qname(aln);
    int chr = aln->core.tid ; //contig name (chromosome)

    if(last==NULL)
      last=strdup(qname);
    if(strcmp(last,qname)!=0){
      if(taxids.size()>0){
	lca=do_lca(taxids,parent);
	fprintf(stdout,"%s\t%d\n",last,lca);
      }
      free(last);
      last=strdup(qname);
    }
    
    int2int::iterator it = i2i.find(chr);
    if(it==i2i.end())
      fprintf(stderr,"\t-> problem finding chrid:%d chrname:%s\n",chr,hdr->target_name[chr]);
    else
      taxids.push_back(it->second);
  }
  if(taxids.size()>0){
    lca=do_lca(taxids,parent);
    fprintf(stdout,"%s\t%d\n",last,lca);
  }
  bam_destroy1(aln);
  sam_close(fp_in);
  
  return ;//0;
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
    }
    int key=atoi(toks[0]);
    rank[key]=strdup(toks[2]);
  }

  fprintf(stderr,"\t-> Number of unique names (column1): %lu from file: %s\n",rank.size(),fname);
}


int main(int argc, char **argv){
  const char *htsfile="new.bam";
  const char *as2tax="nucl_gb.accession2taxid.gz";
  const char *nodes = "nodes.dmp";

  fprintf(stderr,"\t-> as2fax:%s nodes:%d hts:%s\n",as2tax,nodes,htsfile);
  
  samFile *fp_in = hts_open(htsfile,"r"); //open bam file
  bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
  sam_close(fp_in);

  //map of taxid -> taxid
  int2int parent;
  //map of taxid -> rank
  int2char rank;
  
  //map of bamref ->taxid
  int2int i2i= ref2tax(as2tax,bamHdr);
  parse_nodes(nodes,rank,parent);  





  hts(htsfile,i2i,parent,bamHdr);
  return 0;
}
