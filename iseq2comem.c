//	Copyright 2019 Huiguang Yi. All Rights Reservered.
//
//	Licensed under the Apache License, Version 2.0 (the "License");
//	you may not use this file except in compliance with the License.
//	You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
//	Unless required by applicable law or agreed to in writing, software
//	distributed under the License is distributed on an "AS IS" BASIS,
//	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//	See the License for the specific language governing permissions and
//	limitations under the License.
   
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
int reads2mco(char* seqfname,const char *co_dir, char * pipecmd){
#define unit_incrs 1000000
int comp_code_bits = half_ctx_len - drlevel > COMPONENT_SZ ? 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 0 ;
 size_t **cof_count = malloc( component_num * sizeof (size_t*) );
 FILE **outindf = malloc(component_num * sizeof(FILE *));
 FILE **outf = malloc(component_num * sizeof(FILE *));
 char indexfname[PATHLEN]; char combined_cof[PATHLEN];
 size_t cof_count_sz = unit_incrs;
 for(int i=0;i<component_num ;i++)
  {
  cof_count[i] = (size_t *)calloc( cof_count_sz , sizeof(size_t) );
  sprintf(combined_cof,"%s/combco.%d",co_dir,i);
  sprintf(indexfname,"%s/combco.index.%d",co_dir,i);
   if( (outf[i] = fopen(combined_cof,"wb")) == NULL) err(errno,"%s",combined_cof);
    if( (outindf[i] = fopen(indexfname,"wb")) == NULL) err(errno,"%s",indexfname);
  };
 FILE *infp;
 char fas_fname[PATHLEN];
 if(pipecmd[0] != '\0'){
  sprintf(fas_fname,"%s %s",pipecmd,seqfname);
  if( (infp=popen(fas_fname,"r")) == NULL ) err(errno,"reads2mco():%s",fas_fname);
 }
 else
  if( (infp=fopen(seqfname,"r")) == NULL ) err(errno,"reads2mco():%s",fas_fname);;
 char seqin_buff[ READSEQ_BUFFSZ + 1 ];
 int newLen = fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);
 if(! (newLen >0) ) err(errno,"reads2mco():eof or fread error file=%s",seqfname);
 llong base = 1; char ch; int basenum;
 llong tuple = 0LLU, crvstuple = 0LLU,
  unituple, drtuple, pfilter;
 llong readn = 0;
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
   if( readn > cof_count_sz ){
    for(int i=0;i<component_num ;i++){
     size_t *newtmp = (size_t *)realloc(cof_count[i], (cof_count_sz + unit_incrs) * sizeof(size_t));
     if (newtmp != NULL) {
       cof_count[i] = newtmp;
      memset(cof_count[i] + cof_count_sz ,0, unit_incrs * sizeof(size_t) );
     }
     else err(errno,"cof_count[%d] realloc failed",i);
    }
    cof_count_sz += unit_incrs;
   }
   readn++;
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
   cof_count[drtuple % component_num][readn]++;
   unsigned int newid = (unsigned int)( drtuple >> comp_code_bits ) ;
   fwrite( &newid, sizeof(unsigned int),1,outf[(int)( drtuple % component_num )] );
   };
  };
 for(int i=0;i<component_num ;i++){
  llong cumval = 0 ;
  for(llong n=0; n<=readn;n++ ){
   cumval += cof_count[i][n];
   fwrite(&cumval, sizeof(llong),1,outindf[i]);
  }
  fclose(outf[i]);
  fclose(outindf[i]);
 }
 printf("decomposing %s by reads is complete!\n",seqfname);
 return 1;
}
const char gzpipe_cmd[]= "zcat -fc";
llong * fasta2co(char* seqfname, llong *co, char * pipecmd)
{
 llong tuple = 0LLU, crvstuple = 0LLU,
    unituple, drtuple, pfilter;
 memset(co,0LLU,hashsize*sizeof(llong));
 char seqin_buff[ READSEQ_BUFFSZ + 1 ];
 FILE *infp;
 char fas_fname[PATHLEN];
 if(pipecmd[0] != '\0')
  sprintf(fas_fname,"%s %s",pipecmd,seqfname);
 else
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
#define LEN 20000
#define CT_BIT 4
#define CT_MAX 0xfLLU
llong * fastq2co(char* seqfname, llong *co, char *pipecmd, int Q, int M )
{
 if(M >= CT_MAX) err(errno,"fastq2co(): Occurence num should smaller than %d", (int)CT_MAX);
 llong tuple = 0LLU, crvstuple = 0LLU, unituple, drtuple, pfilter;
 memset(co,0LLU,hashsize*sizeof(llong));
 FILE *infp;
 char fq_fname[PATHLEN];
 if(pipecmd[0] != '\0')
    sprintf(fq_fname,"%s %s",pipecmd,seqfname);
  else
    sprintf(fq_fname,"%s %s",gzpipe_cmd,seqfname);
 if( (infp=popen(fq_fname,"r")) == NULL ) err(errno,"fastq2co():%s",fq_fname);
 char *seq = malloc(LEN+10);
 char *qual = malloc(LEN+10);
 fgets(seq,LEN,infp); fgets(seq,LEN,infp);
 fgets(qual,LEN,infp); fgets(qual,LEN,infp);
 llong base = 1; char ch ; int basenum,line_num = 0 ;
 unsigned int keycount =0 ;
 int sl = strlen(seq);
 for(int pos = 0; pos < sl; pos++){
  if(seq[pos] == '\n' ){
   fgets(seq,LEN,infp); fgets(seq,LEN,infp);
   fgets(qual,LEN,infp); fgets(qual,LEN,infp);
   sl = strlen(seq);
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
     if( M == 1) co[n] = (drtuple << CT_BIT) | CT_MAX;
     else co[n] = (drtuple << CT_BIT) + 1LLU;
          if( keycount > hashlimit)
            err(errno,"the context space is too crowd, try rerun the program using -k%d", half_ctx_len + 1);
     break;
    }
    else if ( ( co[n] >> CT_BIT ) == drtuple ) {
           if( (co[n] & CT_MAX) == CT_MAX ) break;
      co[n] += 1LLU;
      if( !((co[n] & CT_MAX) < M) ) co[n]|= CT_MAX ;
      break ;
        };
   };
  };
 }
 printf("%d reads detected\n",line_num);
 pclose(infp);
 free(seq);
  free(qual);
 return co;
};
#define OCCRC_BIT 16
#define OCCRC_MAX 0xffffLLU
llong * fastq2koc (char* seqfname, llong *co, char *pipecmd, int Q)
{
  llong tuple = 0LLU, crvstuple = 0LLU, unituple, drtuple, pfilter;
  memset(co,0LLU,hashsize*sizeof(llong));
  FILE *infp;
  char fq_fname[PATHLEN];
 if(pipecmd[0] != '\0')
    sprintf(fq_fname,"%s %s",pipecmd,seqfname);
  else
    sprintf(fq_fname,"%s %s",gzpipe_cmd,seqfname);
  if( (infp=popen(fq_fname,"r")) == NULL ) err(errno,"fastq2koc():%s",fq_fname);
 char *seq = malloc(LEN+10);
  char *qual = malloc(LEN+10);
  fgets(seq,LEN,infp); fgets(seq,LEN,infp);
  fgets(qual,LEN,infp); fgets(qual,LEN,infp);
  llong base = 1; char ch ; int basenum,line_num = 0 ;
  unsigned int keycount =0 ;
 int sl = strlen(seq);
  for(int pos = 0; pos < sl; pos++){
    if(seq[pos] == '\n' ){
      fgets(seq,LEN,infp); fgets(seq,LEN,infp);
      fgets(qual,LEN,infp); fgets(qual,LEN,infp);
   sl = strlen(seq);
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
          break ;
        };
      };
    };
  }
  pclose(infp);
 free(seq);
 free(qual);
  return co;
};
unsigned int write_fqkoc2files(char* cofilename, llong *co)
{
  int comp_code_bits = half_ctx_len - drlevel > COMPONENT_SZ ? 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 0 ;
  FILE **outf,**abdf;
  outf = malloc(component_num * sizeof(FILE *));
 abdf = malloc(component_num * sizeof(FILE *));
  char cofilename_with_component[PATHLEN];
  for(int i=0;i<component_num ;i++)
  {
    sprintf(cofilename_with_component,"%s.%d",cofilename,i);
    if ( (outf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"write_fqkoc2files()") ;
  sprintf(cofilename_with_component,"%s.%d.a",cofilename,i);
  if ( (abdf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"write_fqkoc2files()") ;
  };
  unsigned int count, wr = 0, newid,compi;
 unsigned short abdc;
  for(count=0;count < hashsize; count++)
  {
    if( co[count] > 0 ) {
   compi = (co[count] >> OCCRC_BIT ) % component_num;
      newid = (unsigned int)(co[count] >> (comp_code_bits + OCCRC_BIT));
      fwrite( &newid, sizeof(newid),1,outf[compi] );
   abdc = co[count] & OCCRC_MAX;
   fwrite( &abdc, sizeof(abdc),1,abdf[compi] );
      wr++;
    }
  }
  for(int i=0;i<component_num ;i++){
    fclose(outf[i]);
  fclose(abdf[i]);
 }
  free(outf);
 free(abdf);
  return wr;
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
#define THREAD_MAX 65536
#define FQ_LEN 4096
llong * mt_shortreads2koc (char* seqfname, llong *co, char *pipecmd,int p){
printf("running mt_shortreads2koc()\n");
char (*fq_buff)[FQ_LEN] = malloc( THREAD_MAX * FQ_LEN );
char tmp[FQ_LEN];
int l;
unsigned int keycount =0 ;
  memset(co,0LLU,hashsize*sizeof(llong));
 FILE *infp;
  char fq_fname[PATHLEN];
  if(pipecmd[0] != '\0') sprintf(fq_fname,"%s %s",pipecmd,seqfname);
 else sprintf(fq_fname,"%s %s",gzpipe_cmd,seqfname);
  if( (infp=popen(fq_fname,"r")) == NULL ) err(errno,"mtfastq2koc():%s",fq_fname);
 while (!feof(infp)){
  for (l = 0; l < THREAD_MAX && fgets(tmp,FQ_LEN,infp) && fgets(fq_buff[l],FQ_LEN,infp) && fgets(tmp,FQ_LEN,infp) && fgets(tmp,FQ_LEN,infp); l++) ;
#pragma omp parallel for num_threads(p) schedule(guided)
  for (int t = 0 ; t < l ; t++ ){
   int base = 1; char ch;
      llong tuple = 0LLU;
      llong crvstuple = 0LLU;
      llong unituple = 0LLU;
   for(int pos = 0; (ch = fq_buff[t][pos]) != '\n'; pos++){
    int basenum = Basemap[(int)ch];
    if(basenum != DEFAULT ){
     tuple = ( ( tuple<< 2 ) | (llong)basenum ) & tupmask ;
     crvstuple = ( crvstuple >> 2 ) + (((llong)basenum^3LLU) << crvsaddmove);
     base++;
    }else{ base = 1;continue;}
    if( base > TL ){
         unituple = tuple < crvstuple ? tuple:crvstuple;
        unsigned int dim_tup = ((unituple & domask) >> ( (half_outctx_len)*2 ) ) ;
     llong pfilter = dim_shuf_arr[dim_tup];
     if( ( pfilter >= dim_end) || (pfilter < dim_start ) ) continue;
     pfilter = pfilter - dim_start;
     llong drtuple = ( ( (unituple & undomask)
      + ( ( unituple & ( ( 1LLU<< ( half_outctx_len*2) ) - 1)) << (TL*2 - half_outctx_len*4) ) )
      >> ( drlevel*4 ) )
      + pfilter ;
     for(unsigned int i=0;i<hashsize;i++){
      unsigned int n = HASH(drtuple, i, hashsize);
      if (co[n] == 0LLU){
#pragma omp atomic write
      co[n] = (drtuple << OCCRC_BIT) + 1LLU ;
#pragma omp atomic
       keycount++;
       if( keycount > hashlimit )
         err(errno,"the context space is too crowd, try rerun the program using -k%d", half_ctx_len + 1);
       break;
      }else if ( ( co[n] >> OCCRC_BIT ) == drtuple ) {
       if( (co[n] & OCCRC_MAX) < OCCRC_MAX )
#pragma omp atomic
         co[n]+=1LLU;
        break ;
      }
     }
    }
   }
  }
 }
 pclose(infp);
 free(fq_buff);
 return co;
}
