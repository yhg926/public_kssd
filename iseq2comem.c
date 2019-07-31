# 1 "iseq2comem.c"
#include "iseq2comem.h"
#include "command_dist.h"
#include "global_basic.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <math.h>
#include <string.h>
#include <err.h>
#include <errno.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>



#define HIBITSET1 0x8000000000000000LLU
#define _64MASK 0xffffffffffffffffLLU
#define H1(K,HASH_SZ) ((K)%(HASH_SZ))
#define H2(K,HASH_SZ) ( 1 + (K) % ( (HASH_SZ) - 1 ) )
#define HASH(K,I,HASH_SZ) ( ( H1(K,HASH_SZ) + I * H2(K,HASH_SZ) ) % HASH_SZ )

static int rand_id;
static int half_ctx_len;
static int half_subctx_len;
static int half_outctx_len;
static int drlevel;

static int comp_bittl;
static int crvsaddmove ;
static llong tupmask;
static int TL ;

static llong domask ;
static llong undomask ;

dim_shuffle_t* dim_shuffle;
static int *dim_shuf_arr ;
static int dim_shuf_arr_len;
static int dim_start;
static int dim_end;

unsigned int hashsize;
unsigned int hashlimit;
int component_num;

void seq2co_global_var_initial(void)
{
 rand_id = dim_shuffle->dim_shuffle_stat.id ;
 half_ctx_len = dim_shuffle->dim_shuffle_stat.k ;
 half_subctx_len = dim_shuffle->dim_shuffle_stat.subk ;
 half_outctx_len = half_ctx_len - half_subctx_len;
  drlevel = dim_shuffle->dim_shuffle_stat.drlevel;
 hashlimit = hashsize * LD_FCTR ;
 printf("rand_id=%d\thalf_ctx_len=%d\thashsize=%d\thashlimit=%d\n",rand_id,half_ctx_len,hashsize,hashlimit);

 component_num = half_ctx_len - drlevel > COMPONENT_SZ ?
         1LU << 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 1 ;

 comp_bittl = 64-4*half_ctx_len;
 crvsaddmove = 4*half_ctx_len-2;
 tupmask = _64MASK >> comp_bittl;
 TL = 2*half_ctx_len ;



 domask = ( (1LLU << ( half_subctx_len*4 ) ) - 1 ) << (2*half_outctx_len);
 undomask = ( (1LLU << (half_outctx_len*2)) - 1 )
         << (2*(half_ctx_len + half_subctx_len));

 dim_shuf_arr = dim_shuffle->shuffled_dim;
 dim_shuf_arr_len = 1LLU << (4*half_subctx_len) ;
 dim_start = 0;
 dim_end = MIN_SUBCTX_DIM_SMP_SZ ;
};


const char gzpipe_cmd[]= "zcat -fc";

llong * fasta2co(char* seqfname, llong *co)
{
 llong tuple = 0LLU, crvstuple = 0LLU,
    unituple, drtuple, pfilter;

 memset(co,0LLU,hashsize*sizeof(llong));
 char seqin_buff[ READSEQ_BUFFSZ + 1 ];
 FILE *infp;
 char fas_fname[PATHLEN];
 sprintf(fas_fname,"%s %s",gzpipe_cmd,seqfname);
 if( (infp=popen(fas_fname,"r")) == NULL ) err(errno,"fasta2co():%s",fas_fname);

 int newLen = fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);
 if(! (newLen >0) ) err(errno,"fastco():eof or fread error file=%s",seqfname);

 llong base = 1; char ch; int basenum;
 unsigned int keycount = 0;

 for(int pos = 0; pos <= newLen; pos++)
 {
  if(pos == newLen){
    newLen = fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);
   if(newLen > 0)
    pos = 0;
   else break;
  };
  ch = seqin_buff[pos];
  basenum = Basemap[(int)ch];

  if(basenum != DEFAULT)
  {
   tuple = ( ( tuple<< 2 ) | (llong)basenum ) & tupmask ;
   crvstuple = ( crvstuple >> 2 ) + (((llong)basenum^3LLU) << crvsaddmove);
   base++;
  }
  else if ( (ch == '\n') || (ch == '\r') ) { continue;}
  else if (isalpha(ch)){ base=1; continue; }
  else if ( ch == '>' )
  {
   while( (pos < newLen ) && ( seqin_buff[pos] != '\n' ) )
   {
    if (pos < newLen - 1)
     pos++;
    else
    {
     newLen = fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);
     if(newLen > 0) pos = -1;
     else err(errno,"fasta2co(): can not find seqences head start from '>' %d",newLen);
    };
   };
      base = 1;
      continue;
  }
  else {

      base=1;
      continue;
    };

  if( base > TL )
  {

   unituple = tuple < crvstuple ? tuple:crvstuple;


   int dim_tup = ((unituple & domask) >> ( (half_outctx_len)*2 ) ) ;
   pfilter = dim_shuf_arr[dim_tup];
   if( ( pfilter >= dim_end) || (pfilter < dim_start ) ) continue;
   pfilter = pfilter - dim_start;
   drtuple = ( ( (unituple & undomask)
       + ( ( unituple & ( ( 1LLU<< ( half_outctx_len*2) ) - 1)) << (TL*2 - half_outctx_len*4) ) )
             >> ( drlevel*4 ) )
             + pfilter ;

   unsigned int i,n ;
     for(i=0;i<hashsize;i++)
     {
       n = HASH(drtuple,i,hashsize);
       if (co[n] == 0)
       {
         co[n] = drtuple;
         keycount++;
         if( keycount > hashlimit)
           err(errno,"the context space is too crowd, try rerun the program using -k%d", half_ctx_len + 1);
         break;
       }
       else if ( co[n] == drtuple )
          break;
   };
  };
 };
  pclose(infp);
  return co;
};



#define LEN 4096
#define CT_BIT 4
#define CT_MAX 0xfLLU

llong * fastq2co(char* seqfname, llong *co, int Q, int M )
{
 if(M >= CT_MAX) err(errno,"fastq2co(): Occurence num should smaller than %d", (int)CT_MAX);

 llong tuple = 0LLU, crvstuple = 0LLU, unituple, drtuple, pfilter;
 memset(co,0LLU,hashsize*sizeof(llong));

 FILE *infp;
 char fq_fname[PATHLEN];
 sprintf(fq_fname,"%s %s",gzpipe_cmd,seqfname);
 if( (infp=popen(fq_fname,"r")) == NULL ) err(errno,"fastq2co():%s",fq_fname);

 char seq[LEN];
 char qual[LEN];

 fgets(seq,LEN,infp); fgets(seq,LEN,infp);
 fgets(qual,LEN,infp); fgets(qual,LEN,infp);

 llong base = 1; char ch ; int basenum,line_num = 0 ;
 unsigned int keycount =0 ;
 for(int pos = 0; pos < strlen(seq); pos++){
  if(seq[pos] == '\n' ){
   fgets(seq,LEN,infp); fgets(seq,LEN,infp);
   fgets(qual,LEN,infp); fgets(qual,LEN,infp);
   line_num+=4;
   if( !feof(infp) ) {
    base = 1;
    pos = -1;
    continue ;
   }
   else break;
  }
  else{
   ch = seq[pos];
   basenum = Basemap[(int)ch];
   if( (basenum != DEFAULT) && ( qual[pos] >= Q ) ){
    tuple = ( ( tuple<< 2 ) | (llong)basenum ) & tupmask ;
    crvstuple = ( crvstuple >> 2 ) + (((llong)basenum^3LLU) << crvsaddmove);
    base++;
   }

   else {




    base = 1;
    continue;
   };
  };

  if( base > TL ){
   unituple = tuple < crvstuple ? tuple:crvstuple;
   int dim_tup = ((unituple & domask) >> ( (half_outctx_len)*2 ) ) ;
      pfilter = dim_shuf_arr[dim_tup];
      if( ( pfilter >= dim_end) || (pfilter < dim_start ) ) continue;
   pfilter = pfilter - dim_start;
      drtuple = ( ( (unituple & undomask)
              + ( ( unituple & ( ( 1LLU<< ( half_outctx_len*2) ) - 1)) << (TL*2 - half_outctx_len*4) ) )
              >> ( drlevel*4 ) )
              + pfilter ;

   unsigned int i,n ;
   for(i=0;i<hashsize;i++) {
    n = HASH(drtuple, i, hashsize);
    if (co[n] == 0LLU){
     co[n] = (drtuple << CT_BIT) + 1LLU ;
     keycount++;
          if( keycount > hashlimit)
            err(errno,"the context space is too crowd, try rerun the program using -k%d", half_ctx_len + 1);
     break;
    } else if ( ( co[n] >> CT_BIT ) == drtuple ) {

      if( (co[n] & CT_MAX) < M )
       co[n]+=1LLU;
      else
       co[n]|= CT_MAX ;
     break ;
        };
   };
  };
 }
 pclose(infp);
 return co;
};





llong * fastq2koc (char* seqfname, llong *co, int Q)
{

  llong tuple = 0LLU, crvstuple = 0LLU, unituple, drtuple, pfilter;
  memset(co,0LLU,hashsize*sizeof(llong));

  FILE *infp;
  char fq_fname[PATHLEN];
  sprintf(fq_fname,"%s %s",gzpipe_cmd,seqfname);
  if( (infp=popen(fq_fname,"r")) == NULL ) err(errno,"fastq2co():%s",fq_fname);

  char seq[LEN];
  char qual[LEN];

  fgets(seq,LEN,infp); fgets(seq,LEN,infp);
  fgets(qual,LEN,infp); fgets(qual,LEN,infp);

  llong base = 1; char ch ; int basenum,line_num = 0 ;
  unsigned int keycount =0 ;
  for(int pos = 0; pos < strlen(seq); pos++){
    if(seq[pos] == '\n' ){
      fgets(seq,LEN,infp); fgets(seq,LEN,infp);
      fgets(qual,LEN,infp); fgets(qual,LEN,infp);
      line_num+=4;
      if( !feof(infp) ) {
        base = 1;
        pos = -1;
        continue ;
      }
      else break;
    }
    else{
      ch = seq[pos];
      basenum = Basemap[(int)ch];
      if( (basenum != DEFAULT) && ( qual[pos] >= Q ) ){
        tuple = ( ( tuple<< 2 ) | (llong)basenum ) & tupmask ;
        crvstuple = ( crvstuple >> 2 ) + (((llong)basenum^3LLU) << crvsaddmove);
        base++;
      }

      else {
        base = 1;
        continue;
      };
    };

    if( base > TL ){
      unituple = tuple < crvstuple ? tuple:crvstuple;
      unsigned int dim_tup = ((unituple & domask) >> ( (half_outctx_len)*2 ) ) ;
      pfilter = dim_shuf_arr[dim_tup];
      if( ( pfilter >= dim_end) || (pfilter < dim_start ) ) continue;
      pfilter = pfilter - dim_start;
      drtuple = ( ( (unituple & undomask)
              + ( ( unituple & ( ( 1LLU<< ( half_outctx_len*2) ) - 1)) << (TL*2 - half_outctx_len*4) ) )
              >> ( drlevel*4 ) )
              + pfilter ;

      unsigned int i,n ;
      for(i=0;i<hashsize;i++) {
        n = HASH(drtuple, i, hashsize);
        if (co[n] == 0LLU){
          co[n] = (drtuple << OCCRC_BIT) + 1LLU ;
          keycount++;
          if( keycount > hashlimit )
            err(errno,"the context space is too crowd, try rerun the program using -k%d", half_ctx_len + 1);
          break;
        } else if ( ( co[n] >> OCCRC_BIT ) == drtuple ) {

            if( (co[n] & OCCRC_MAX) < OCCRC_MAX )
              co[n]+=1LLU;
            else
              co[n]|= OCCRC_MAX ;
          break ;
        };
      };
    };
  }
  pclose(infp);
  return co;
};

llong write_fqkoc2file(char* cofilename, llong *co)
{
  int comp_code_bits = half_ctx_len - drlevel > COMPONENT_SZ ? 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 0 ;
  FILE **outf;
  outf = malloc(component_num * sizeof(FILE *));
  char cofilename_with_component[PATHLEN];
  for(int i=0;i<component_num ;i++)
  {
    sprintf(cofilename_with_component,"%s.%d",cofilename,i);
    if ( (outf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"write_fqkoc2file()") ;
  };

  unsigned int count, wr = 0;
 llong newid;
  for(count=0;count < hashsize; count++)
  {
    if( co[count] > 0 ) {
   newid = ((co[count] >> (comp_code_bits + OCCRC_BIT)) << OCCRC_BIT ) | (co[count] & OCCRC_MAX) ;
      fwrite( &newid, sizeof(llong),1,outf[ (co[count] >> OCCRC_BIT ) % component_num ] );
      wr++;
    }
  }

  for(int i=0;i<component_num ;i++)
    fclose(outf[i]);

 free(outf);
  return wr;
};


llong write_fqco2file(char* cofilename, llong *co)
{
 int comp_code_bits = half_ctx_len - drlevel > COMPONENT_SZ ? 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 0 ;
 FILE **outf;
 outf = malloc(component_num * sizeof(FILE *));
 char cofilename_with_component[PATHLEN];
 for(int i=0;i<component_num ;i++)
  {
    sprintf(cofilename_with_component,"%s.%d",cofilename,i);
    if ( (outf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"write_fqco2file()") ;
  };
 unsigned int count, wr = 0, newid;
 for(count=0;count < hashsize; count++)
 {
  if( (co[count] & CT_MAX) == CT_MAX ) {
   newid = (unsigned int)( co[count] >> (comp_code_bits + CT_BIT) ) ;
   fwrite( &newid, sizeof(unsigned int),1,outf[(int)( (co[count] >> CT_BIT ) % component_num )] );
      wr++;
  }
 }
 for(int i=0;i<component_num ;i++)
   fclose(outf[i]);
 free(outf);
 return wr;
};



llong wrt_co2cmpn_use_inn_subctx(char* cofilename, llong *co)
{
 int comp_code_bits = half_ctx_len - drlevel > COMPONENT_SZ ? 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 0 ;
 FILE **outf;
 outf = malloc(component_num * sizeof(FILE *));
 char cofilename_with_component[PATHLEN];
  for(int i=0;i<component_num ;i++)
  {
    sprintf(cofilename_with_component,"%s.%d",cofilename,i);
    if ( (outf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"wrt_co2cmpn_use_inn_subctx()") ;
  };
  unsigned int count, wr = 0, newid;
 for(count=0;count < hashsize; count++)
 {
  if( co[count] != 0 )
  {
   newid = (unsigned int)( co[count] >> comp_code_bits ) ;
   fwrite( &newid, sizeof(unsigned int),1,outf[(int)( co[count] % component_num )] );
   wr++;
  }
 }
 for(int i=0;i<component_num ;i++)
   fclose(outf[i]);

 free(outf);
  return wr;
};
