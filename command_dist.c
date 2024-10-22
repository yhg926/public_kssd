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
   
#include "command_dist.h"
#include "command_shuffle.h"
#include "iseq2comem.h"
#include "co2mco.h"
#include "mytime.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <math.h>
#include <tgmath.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>
#include <stdbool.h>
#include <malloc.h>
#include <libgen.h>
#ifdef _OPENMP
 #include <omp.h>
#endif
mem_dispatch_t mem_dispatch = {0,0,0,0,0,0,0, OBJ_TRADITION , STD , MAKE_REF ,0};
mem_usage_stat_t mem_usage_stat = { 0, 0, 1e7 };
static llint real_time_mem;
const char co_dstat[] = "cofiles.stat";
const char mco_dstat[] = "mcofiles.stat";
const char distoutdir[] = "dist";
typedef int mco_co_dist_t[BIN_SZ];
static inline unsigned int * mco_co_mmpdist_core(gidobj_t** unit_arrmco, char *co_fcode_in, unsigned int *ctx_obj_ct_in );
static inline void mco_co_dist_core(gidobj_t** unit_arrmco, char *co_fcode_in, int bin_sz, mco_co_dist_t shared_ctx_num_in);
#define LINE_LEN 1024
typedef struct prt_line { int len; char line[LINE_LEN]; } prt_line_t ;
static inline void output_ctrl (unsigned int X_size, unsigned int XnY_size, print_ctrl_t* outfield, char *rname, prt_line_t* linebuf ) ;
int dist_dispatch(dist_opt_val_t *opt_val)
{
 if( opt_val->mmry == 0 )
  opt_val->mmry = get_sys_mmry();
 llint base_mem = (llint)( (opt_val->mmry - get_sys_mmry()) * BBILLION ) - mem_usage_stat.others ;
 real_time_mem = base_mem + (llint)(get_sys_mmry()*BBILLION) ;
 int p_fit_mem = 1;
 if( ( opt_val->refpath[0] != '\0' ) )
  {
  const char *refco_dstat_fpath = test_get_fullpath(opt_val->refpath,co_dstat);
  const char *refmco_dstat_fpath = test_get_fullpath(opt_val->refpath,mco_dstat);
  if( (refco_dstat_fpath == NULL ) && ( refmco_dstat_fpath == NULL) ){
   infile_tab_t *ref_stat = dist_organize_refpath(opt_val);
   if(ref_stat->infile_num == 0)
    err(errno,"no valid input .fas/.fq file or absent %s | %s in %s \n", co_dstat, mco_dstat, opt_val->refpath ) ;
   infile_fmt_count_t* ref_fmt_count;
      ref_fmt_count = infile_fmt_count(ref_stat);
   if( (ref_stat->infile_num == 0) || (ref_fmt_count->fasta + ref_fmt_count->fastq != ref_stat->infile_num ))
    err(errno,"not a valid input files: make sure input files are .fas/.fq format"
     "or .co.num |.mco.num files with %s | %s in %s \n", co_dstat, mco_dstat, opt_val->refpath );
   const char *dist_refco_dir = opt_val->outdir;
   mkdir(dist_refco_dir,0777);
   int *shuffled_refname_ind = shuffleN( ref_stat->infile_num, 0 );
   mem_usage_stat.input_file_name_sz = ref_stat->infile_num * ( sizeof(llong) + PATHLEN*sizeof(char) );
real_time_mem -= mem_usage_stat.input_file_name_sz;
   dim_shuffle = get_dim_shuffle(opt_val);
    hashsize = get_hashsz(dim_shuffle);
    seq2co_global_var_initial();
    mem_usage_stat.shuffled_subctx_arr_sz = (1LLU << 4*( dim_shuffle->dim_shuffle_stat.subk) )*sizeof(int);
real_time_mem -= mem_usage_stat.shuffled_subctx_arr_sz;
   p_fit_mem = real_time_mem / ( (hashsize + 1) * sizeof(llong));
   if( opt_val->p < p_fit_mem )
    p_fit_mem = opt_val->p;
   else if(p_fit_mem < 1)
    err(errno,"dist_dispatch():\n"
    " Kmer hashing need mem.(%lf G) exceed the mem. system or user provide (%lf G)\n"
    " user can either consider specify more mem.(-m ) or use smaller k value ( -k)\n"
    " or increase dimension reduction level ( -L)\n",
    (double)hashsize * sizeof(llong)/1e9, opt_val->mmry);
   size_t phash_mem = (llong)p_fit_mem * (hashsize + 1) * sizeof(llong);
real_time_mem -= phash_mem;
   opt_val->abundance = false;
   const char *refcostat = run_stageI(opt_val, ref_stat, shuffled_refname_ind, dist_refco_dir, p_fit_mem);
   run_stageII(refcostat, opt_val->outdir , opt_val->p);
   free(dim_shuffle->shuffled_dim);
real_time_mem += mem_usage_stat.shuffled_subctx_arr_sz;
  }
  else if ( refco_dstat_fpath != NULL) {
   if (refmco_dstat_fpath == NULL)
    run_stageII(refco_dstat_fpath, opt_val->refpath ,opt_val->p);
  }
  else if(refmco_dstat_fpath != NULL) {;}
  free((char*)refco_dstat_fpath);
  free((char*)refmco_dstat_fpath);
 }
 if ( (opt_val->num_remaining_args >0) || (opt_val->fpath[0] != '\0' ) )
 {
  const char *qryco_dstat_fpath = NULL;
  const char *qrymco_dstat_fpath = NULL;
  if((opt_val->pipecmd[0]=='\0') && (opt_val->num_remaining_args >0)){
   qryco_dstat_fpath = test_get_fullpath(opt_val->remaining_args[0],co_dstat);
   qrymco_dstat_fpath = test_get_fullpath(opt_val->remaining_args[0],mco_dstat);
  }
  if(opt_val->refpath[0] != '\0'){
   const char *ref_db = test_get_fullpath(opt_val->refpath, mco_dstat);
   if(ref_db==NULL) err(errno,"need speficy the ref-sketch path for -r to run the query-ref search model");
   FILE * ref_mco_stat_fp;
    if ((ref_mco_stat_fp = fopen(ref_db,"rb")) == NULL) err(errno,"mco stat file:%s",ref_db );
   mco_dstat_t mco_ref_dstat;
   fread( &mco_ref_dstat, sizeof(mco_dstat_t),1,ref_mco_stat_fp );
   fclose(ref_mco_stat_fp);
   if( qryco_dstat_fpath != NULL ){
    FILE *qry_co_stat_fp;
    if (( qry_co_stat_fp = fopen(qryco_dstat_fpath,"rb")) == NULL) err(errno,"qry co stat file:%s",qryco_dstat_fpath);
    co_dstat_t co_qry_dstat;
    fread( &co_qry_dstat, sizeof(co_dstat_t), 1, qry_co_stat_fp);
    if(co_qry_dstat.shuf_id != mco_ref_dstat.shuf_id )
     err(errno, "qry shuf_id: %d not match ref shuf_id: %d\ntry regenerate .co dir and feed -s the .shuf"
           "file used to generated ref database",co_qry_dstat.shuf_id,mco_ref_dstat.shuf_id);
    else if(co_qry_dstat.comp_num != mco_ref_dstat.comp_num)
     err(errno, "qry comp_num: %d not match ref comp_num: %d",co_qry_dstat.comp_num, mco_ref_dstat.comp_num);
    mco_cbdco_nobin_dist(opt_val);
    }
   else if (qrymco_dstat_fpath != NULL)
    err(errno,"when -r specified, the query sould not be .mco format, the valid query format shoulde be .fas/.fq file or .co");
   else{
    infile_tab_t *infile_stat = dist_organize_infiles(opt_val);
    infile_fmt_count_t* qry_fmt_count;
      qry_fmt_count = infile_fmt_count(infile_stat);
    bool is_valid_fas_fq_in = (infile_stat->infile_num != 0) &&
       (qry_fmt_count->fasta + qry_fmt_count->fastq == infile_stat->infile_num );
    if(is_valid_fas_fq_in){
    }
    else err(errno,"please specify valid query genomes seq or .co file for database search");
   }
  }
  else if( qryco_dstat_fpath != NULL) {
   if( opt_val->num_remaining_args == 1 ){
       run_stageII(qryco_dstat_fpath,opt_val->outdir, opt_val->p);
   }
   else if (opt_val->num_remaining_args > 1){
    combine_queries(opt_val);
   }
    }
  else if ( (qrymco_dstat_fpath != NULL) && (opt_val->num_remaining_args > 1 )){
  }
  else {
    infile_tab_t *infile_stat = dist_organize_infiles(opt_val);
        infile_fmt_count_t* qry_fmt_count;
        qry_fmt_count = infile_fmt_count(infile_stat);
        bool is_valid_fas_fq_in = (infile_stat->infile_num != 0) &&
          (qry_fmt_count->fasta + qry_fmt_count->fastq == infile_stat->infile_num );
        if(is_valid_fas_fq_in || (opt_val->pipecmd[0] != '\0') ){
          const char *dist_co_dir = opt_val->outdir;
     mkdir(dist_co_dir,0777);
     int *shuffled_refname_ind = shuffleN( infile_stat->infile_num, 0 );
     mem_usage_stat.input_file_name_sz = infile_stat->infile_num * ( sizeof(llong) + PATHLEN*sizeof(char) );
real_time_mem -= mem_usage_stat.input_file_name_sz;
     dim_shuffle = get_dim_shuffle(opt_val);
        hashsize = get_hashsz(dim_shuffle);
        seq2co_global_var_initial();
        mem_usage_stat.shuffled_subctx_arr_sz = (1LLU << 4*( dim_shuffle->dim_shuffle_stat.subk) )*sizeof(int);
real_time_mem -= mem_usage_stat.shuffled_subctx_arr_sz;
     p_fit_mem = real_time_mem / ( (hashsize + 1) * sizeof(llong));
     if( opt_val->p < p_fit_mem ) p_fit_mem = opt_val->p;
        else if(p_fit_mem < 1)
          err(errno,"dist_dispatch():\n"
          " Kmer hashing need mem.(%lf G) exceed the mem. system or user provide (%lf G)\n"
          " user can either consider specify more mem.(-m ) or use smaller k value ( -k)\n"
          " or increase dimension reduction level ( -L)\n",
          (double)hashsize * sizeof(llong)/1e9, opt_val->mmry);
     size_t phash_mem = (llong)p_fit_mem * (hashsize + 1) * sizeof(llong);
real_time_mem -= phash_mem ;
        run_stageI(opt_val,infile_stat,shuffled_refname_ind,dist_co_dir,p_fit_mem);
        }
    else err(errno,"not valid raw seq format");
  }
 }
 return 0;
};
dim_shuffle_t *get_dim_shuffle( dist_opt_val_t *opt_val_in )
{
 char shuf_infile_name_prefix[PATHLEN+9];
  char shuf_infile_name[PATHLEN];
  strcpy(shuf_infile_name, opt_val_in->dr_file);
  if( strcmp( shuf_infile_name, "" ) == 0 )
  {
    srand ( time(NULL) );
    dim_shuffle_stat_t dim_shuffle_stat =
    {
      rand(),
      opt_val_in->k,
      opt_val_in->dr_level + add_len_drlevel2subk(),
      opt_val_in->dr_level,
    };
  struct stat outd;
  if( (stat(opt_val_in->outdir, &outd) != 0 ) || (! S_ISDIR(outd.st_mode) ) )
   mkdir(opt_val_in->outdir,0777);
  sprintf(shuf_infile_name_prefix, "%s/default",opt_val_in->outdir);
    write_dim_shuffle_file( &dim_shuffle_stat,shuf_infile_name_prefix);
    sprintf(shuf_infile_name,"%s.shuf",shuf_infile_name_prefix);
  };
  return read_dim_shuffle_file(shuf_infile_name);
};
int get_hashsz(dim_shuffle_t *dim_shuffle_in )
{
 llong ctx_space_sz = 1LLU << 4*( dim_shuffle_in->dim_shuffle_stat.k - dim_shuffle_in->dim_shuffle_stat.drlevel );
  int primer_ind = 4*( dim_shuffle_in->dim_shuffle_stat.k - dim_shuffle_in->dim_shuffle_stat.drlevel ) - CTX_SPC_USE_L - 7;
  if(primer_ind < 0 || primer_ind > 24 ){
    int k_add = primer_ind < 0 ? (1 + (0-primer_ind)/4) : - (1 + ( primer_ind - 24 )/4) ;
    err(errno,"get_hashsz(): primer_ind: %d out of range(0 ~ 24), by formula:\n"
    "int primer_ind = 4*(opt_val->k - dim_shuffle->dim_shuffle_stat.drlevel) - CTX_SPC_USE_L - 7\n"
    "this might caused by too small or too large k\n"
    "kmer length = %d\n"
    "dim reduction level = %d\n"
    "ctx_space size = %llu\n"
    "CTX space usage limit = %lf\n\n"
    "try rerun the program with option -k = %d",
    primer_ind, dim_shuffle_in->dim_shuffle_stat.k, dim_shuffle_in->dim_shuffle_stat.drlevel,ctx_space_sz,
       (double)1/(1 << CTX_SPC_USE_L),dim_shuffle_in->dim_shuffle_stat.k + k_add );
  };
  int hashsize_get = primer[primer_ind];
 return hashsize_get ;
}
const char* test_get_fullpath(const char *parent_path, const char *dstat_f)
{
 struct stat path_stat;
  if( stat(parent_path, &path_stat) < 0 )
  err(errno,"test_get_fullpath()::%s",parent_path);
  if( S_ISDIR(path_stat.st_mode) ){
  char* fullpath = malloc(PATHLEN+1);
  sprintf((char *)fullpath,"%s/%s", parent_path, dstat_f);
  FILE *fp;
  if ( (fp = fopen(fullpath,"rb")) != NULL ){
   fclose(fp);
   return fullpath;
  }
  else{
   free((char*)fullpath);
   return NULL;
  }
 }
 else
  return NULL;
};
const char * run_stageI (dist_opt_val_t *opt_val, infile_tab_t *seqfile_stat,
      int* shuffled_seqfname_ind, const char *co_dir, int p_fit_mem)
{
 llong **CO = malloc( p_fit_mem * sizeof(llong *) );
  for(int i = 0; i< p_fit_mem; i++ ){
   CO[i] = (llong *)malloc(hashsize * sizeof(llong) );
 }
 llong all_ctx_ct = 0 ;
 ctx_obj_ct_t *ctx_ct_list = malloc(sizeof(ctx_obj_ct_t) * seqfile_stat->infile_num);
 if(opt_val->byread){
  for(int i = 0; i< seqfile_stat->infile_num; i++){
   char* seqfname = seqfile_stat->organized_infile_tab[ shuffled_seqfname_ind[i] ].fpath;
   printf("decomposing %s by reads\n",seqfname) ;
   reads2mco(seqfname, co_dir, opt_val->pipecmd);
  }
 }
 else {
  int num_threads = seqfile_stat->infile_num > p_fit_mem ? p_fit_mem : 1;
  int fcount = 0;
#pragma omp parallel for num_threads(num_threads) reduction(+:all_ctx_ct) schedule(guided)
   for(int i = 0; i< seqfile_stat->infile_num; i++){
     int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      char* seqfname = seqfile_stat->organized_infile_tab[ shuffled_seqfname_ind[i] ].fpath;
      char cofname[PATHLEN];
      sprintf(cofname,"%s/%d.co",co_dir,i);
      llong *co;
      if(isOK_fmt_infile(seqfname,fastq_fmt,FQ_FMT_SZ) || opt_val->pipecmd[0]!='\0'){
    if(opt_val->abundance){
     co = mt_shortreads2koc(seqfname,CO[tid],opt_val->pipecmd, p_fit_mem);
      ctx_ct_list[i] = write_fqkoc2files(cofname,co);
    }
    else{
     co = fastq2co(seqfname,CO[tid],opt_val->pipecmd,opt_val->kmerqlty,opt_val->kmerocrs);
         ctx_ct_list[i] = write_fqco2file(cofname,co);
    }
       }
       else{
    if(opt_val->abundance) {
     opt_val->abundance = 0;
     printf("Warning: close abundance mode (-A) since non-fastq file input.\n");
    }
    if(opt_val->u)
     co = uniq_fasta2co(seqfname,CO[tid],opt_val->pipecmd);
    else
      co = fasta2co(seqfname,CO[tid],opt_val->pipecmd);
        ctx_ct_list[i] = wrt_co2cmpn_use_inn_subctx(cofname,co);
       }
   all_ctx_ct += ctx_ct_list[i] ;
#pragma omp atomic
   fcount++;
   printf("%d/%d decomposing %s\r",fcount,seqfile_stat->infile_num,seqfname) ;
  }
  printf("\n");
#pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
    for(int c = 0; c < component_num; c++){
   size_t *cof_index_in_cbdco = malloc( (seqfile_stat->infile_num + 1) * sizeof(size_t) );
      char tmpfname[PATHLEN];
      FILE *com_cofp,*indexfp, *com_abund_fp;
     sprintf(tmpfname,"%s/%s.%d",co_dir,skch_prefix,c);
   if( (com_cofp = fopen(tmpfname,"wb")) == NULL) err(errno,"%s",tmpfname);
      sprintf(tmpfname,"%s/%s.%d",co_dir,idx_prefix,c);
      if( (indexfp = fopen(tmpfname,"wb")) == NULL) err(errno,"%s",tmpfname);
   if(opt_val->abundance) {
    sprintf(tmpfname,"%s/%s.%d.a",co_dir,skch_prefix,c);
    com_abund_fp = fopen(tmpfname,"wb");
    if( com_abund_fp == NULL) err(errno,"%s",tmpfname);
   }
      void *tmp_mem = malloc( LD_FCTR * 2 * (1LLU << (4*COMPONENT_SZ - CTX_SPC_USE_L)) *sizeof(unsigned int) );
      struct stat tmpstat;
      FILE *tmpfp;
      cof_index_in_cbdco[0] = 0;
      for(int i = 0; i< seqfile_stat->infile_num; i++){
       sprintf(tmpfname,"%s/%d.co.%d", co_dir, i, c);
        if( (tmpfp = fopen(tmpfname,"rb")) == NULL) err(errno,"%s",tmpfname);
        stat(tmpfname,&tmpstat);
        int tmpkmerct = tmpstat.st_size/sizeof(unsigned int);
        cof_index_in_cbdco[i+1] = (size_t)cof_index_in_cbdco[i] + tmpkmerct;
        fread(tmp_mem,tmpstat.st_size,1,tmpfp);
        fwrite(tmp_mem,tmpstat.st_size,1,com_cofp);
        fclose(tmpfp);
        remove(tmpfname);
    if(opt_val->abundance) {
     sprintf(tmpfname,"%s/%d.co.%d.a", co_dir, i, c);
     if( (tmpfp = fopen(tmpfname,"rb")) == NULL) err(errno,"%s",tmpfname);
     fread(tmp_mem, sizeof(unsigned short)*tmpkmerct, 1 , tmpfp);
     fwrite(tmp_mem, sizeof(unsigned short)*tmpkmerct, 1 , com_abund_fp );
     fclose(tmpfp);
     remove(tmpfname);
    }
      }
      fclose(com_cofp);
   if(opt_val->abundance) fclose(com_abund_fp) ;
   free(tmp_mem);
      fwrite(cof_index_in_cbdco,sizeof(size_t), seqfile_stat->infile_num + 1, indexfp);
      fclose(indexfp);
      free(cof_index_in_cbdco);
    }
 }
  for(int i = 0; i< p_fit_mem; i++ ) free(CO[i]);
  free(CO);
  co_dstat_t co_dstat_wrout;
 co_dstat_wrout.shuf_id = dim_shuffle->dim_shuffle_stat.id ;
 co_dstat_wrout.koc = opt_val->abundance ;
 co_dstat_wrout.kmerlen = dim_shuffle->dim_shuffle_stat.k * 2;
 co_dstat_wrout.dim_rd_len = dim_shuffle->dim_shuffle_stat.drlevel * 2 ;
  co_dstat_wrout.comp_num = component_num ;
  co_dstat_wrout.infile_num = seqfile_stat->infile_num;
 co_dstat_wrout.all_ctx_ct = all_ctx_ct;
  char *co_dstat_fullname = malloc(PATHLEN*sizeof(char) );
  sprintf(co_dstat_fullname, "%s/%s",co_dir,co_dstat);
  FILE *coutfp;
  if ( ( coutfp = fopen(co_dstat_fullname,"wb")) == NULL ) err(errno,"%s",co_dstat_fullname);
   fwrite(&co_dstat_wrout,sizeof(co_dstat_wrout),1,coutfp);
 fwrite(ctx_ct_list,sizeof(ctx_obj_ct_t),co_dstat_wrout.infile_num,coutfp);
 free(ctx_ct_list);
 for(int i = 0; i< co_dstat_wrout.infile_num; i++)
  fwrite(seqfile_stat->organized_infile_tab[ shuffled_seqfname_ind[i] ].fpath,PATHLEN,1,coutfp);
  fclose(coutfp);
 return (const char *)co_dstat_fullname;
}
void run_stageII(const char * co_dstat_fpath, const char* dist_mco_dir, int p_fit_mem)
{
    const char* dist_co_dir = get_pathname(co_dstat_fpath,co_dstat);
   if ( mkdir(dist_mco_dir, 0700) ) {
    if(errno == EEXIST)
     printf("Warning: write mco file to an exists outdir:%s\n",dist_mco_dir);
    else
     err(errno,"run_stageII(): mkdir %s error", dist_mco_dir);
   }
      const char* mco_dstat_fpath = malloc(PATHLEN*sizeof(char));
   sprintf((char*)mco_dstat_fpath,"%s/%s",dist_mco_dir,mco_dstat);
   FILE *co_stat_fp,*mco_stat_fp;
   if( ( co_stat_fp = fopen(co_dstat_fpath,"rb")) == NULL ) err(errno,"run_stageII(():%s",co_dstat_fpath);
   co_dstat_t co_dstat_readin;
   fread( &co_dstat_readin, sizeof(co_dstat_t),1,co_stat_fp );
   if( ( mco_stat_fp = fopen(mco_dstat_fpath,"wb")) == NULL ) err(errno,"run_stageII(():%s",mco_dstat_fpath);
   mco_dstat_t mco_dstat_writeout;
   mco_dstat_writeout.shuf_id = co_dstat_readin.shuf_id;
   mco_dstat_writeout.kmerlen = co_dstat_readin.kmerlen ;
      mco_dstat_writeout.dim_rd_len = co_dstat_readin.dim_rd_len ;
   mco_dstat_writeout.infile_num = co_dstat_readin.infile_num;
   mco_dstat_writeout.comp_num = co_dstat_readin.comp_num ;
   fwrite(&mco_dstat_writeout,sizeof(mco_dstat_writeout),1, mco_stat_fp );
   ctx_obj_ct_t *tmp_ctx_ct = malloc(sizeof(ctx_obj_ct_t)*co_dstat_readin.infile_num);
   fread(tmp_ctx_ct,sizeof(ctx_obj_ct_t),co_dstat_readin.infile_num,co_stat_fp);
   fwrite(tmp_ctx_ct,sizeof(ctx_obj_ct_t),co_dstat_readin.infile_num,mco_stat_fp);
   char (*tmpname)[PATHLEN] = malloc( PATHLEN * co_dstat_readin.infile_num );
   fread(tmpname,PATHLEN,co_dstat_readin.infile_num,co_stat_fp);
   fwrite(tmpname,PATHLEN,co_dstat_readin.infile_num,mco_stat_fp);
   free(tmp_ctx_ct);
   free(tmpname);
   fclose(co_stat_fp);
   fclose(mco_stat_fp);
   combco2mco(dist_mco_dir,dist_co_dir,co_dstat_readin.infile_num,co_dstat_readin.comp_num,p_fit_mem);
   free((char*)dist_co_dir );
   free((char*)mco_dstat_fpath);
}
int alp_size = 4 ;
static ctx_obj_ct_t initial_dist[BIN_SZ];
static int ref_seq_num,qry_seq_num,kmerlen,dim_reduct_len ;
void mco_co_dist( char *refmco_dname, char *qryco_dname, const char *distout_dir, int p_fit_mem)
{
 FILE *refmco_dstat_fp, *qryco_dstat_fp;
 char *refmco_dstat_fpath = malloc(PATHLEN*sizeof(char));
 char *qryco_dstat_fpath = malloc(PATHLEN*sizeof(char));
 sprintf(refmco_dstat_fpath,"%s/%s",refmco_dname,mco_dstat);
 sprintf(qryco_dstat_fpath,"%s/%s",qryco_dname,co_dstat);
 if( (refmco_dstat_fp = fopen(refmco_dstat_fpath,"rb")) == NULL )
  err(errno,"need provied mco dir path for mco_co_dist() arg 1. refmco_dstat_fpath");
 if( (qryco_dstat_fp = fopen(qryco_dstat_fpath,"rb")) == NULL )
  err(errno,"need provied co dir path for mco_co_dist() arg 2.  qryco_dstat_fpath");
 mco_dstat_t mco_dstat_readin ;
  co_dstat_t co_dstat_readin ;
  fread(&mco_dstat_readin,sizeof(mco_dstat_readin),1,refmco_dstat_fp);
 fread(&co_dstat_readin,sizeof(co_dstat_readin),1,qryco_dstat_fp);
 unsigned int * qry_ctx_ct_list = malloc(co_dstat_readin.infile_num * sizeof(unsigned int));
 unsigned int * ref_ctx_ct_list = malloc(mco_dstat_readin.infile_num * sizeof(unsigned int));
 fread(qry_ctx_ct_list,sizeof(unsigned int),co_dstat_readin.infile_num,qryco_dstat_fp);
 fread(ref_ctx_ct_list,sizeof(unsigned int),mco_dstat_readin.infile_num,refmco_dstat_fp);
 char (*cofname)[PATHLEN] = malloc(co_dstat_readin.infile_num * PATHLEN);
 char (*mcofname)[PATHLEN] = malloc(mco_dstat_readin.infile_num * PATHLEN);
 fread(cofname,PATHLEN,co_dstat_readin.infile_num,qryco_dstat_fp);
 fread(mcofname,PATHLEN,mco_dstat_readin.infile_num,refmco_dstat_fp);
 fclose(refmco_dstat_fp);
  fclose(qryco_dstat_fp);
 if( !(mco_dstat_readin.comp_num == co_dstat_readin.comp_num) )
  err(errno,"query args not match ref args: ref.comp_num = %d vs. %d = qry.comp_num",
    mco_dstat_readin.comp_num, co_dstat_readin.comp_num);
 if(!(mco_dstat_readin.shuf_id == co_dstat_readin.shuf_id))
  err(errno,"query args not match ref args: ref.shuf_id = %d vs. %d = qry.shuf_id",
    mco_dstat_readin.shuf_id, co_dstat_readin.shuf_id);
 int ref_bin_num = mco_dstat_readin.infile_num / BIN_SZ;
 int binsz;
 for(int i=0; i<=ref_bin_num;i++ ){
  if( i == ref_bin_num ){
     binsz = mco_dstat_readin.infile_num % BIN_SZ;
          if(binsz == 0) continue;
    }else binsz = BIN_SZ;
#pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
  for( int k=0; k< co_dstat_readin.infile_num; k++){
   if(qry_ctx_ct_list[k]==0){
    warnx("%dth co file is empty",k);
    continue;
   }
   char dist_fcode[PATHLEN];
   sprintf(dist_fcode,"%s/%d.%d.dist", distout_dir, i, k );
   FILE *distfp;
   if ( (distfp = fopen(dist_fcode,"wb")) == NULL) err(errno,"mco_co_dist()::%s",dist_fcode);
   fwrite(initial_dist, sizeof(ctx_obj_ct_t), binsz, distfp);
   fclose(distfp);
  }
  for ( int j = 0; j < mco_dstat_readin.comp_num; j++ ){
   char mco_fcode[PATHLEN];
   sprintf(mco_fcode,"%s/%d.mco.%d",refmco_dname,i,j);
   gidobj_t** unit_arrmco_readin = read_unit_arrmco_file(mco_fcode);
#pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
   for(int k=0; k< co_dstat_readin.infile_num; k++){
    if(qry_ctx_ct_list[k]==0) continue;
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
    char co_fcode[PATHLEN]; char dist_fcode[PATHLEN];
    sprintf(co_fcode,"%s/%d.%d.co.%d",qryco_dname,k/BIN_SZ, k % BIN_SZ, j);
    sprintf(dist_fcode,"%s/%d.%d.dist",distout_dir, i , k );
    int fd;
    if( ( (fd = open(dist_fcode,O_RDWR, 0600) ) == -1) )
         err(errno,"mco_co_dist()::distfile = %s[tid = %d]",dist_fcode,tid);
    ctx_obj_ct_t * ctx_obj_ct = mmap(NULL, binsz*(sizeof(ctx_obj_ct_t)), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if(ctx_obj_ct == MAP_FAILED) err(errno,"ctx_obj_ct mmap error");
    ctx_obj_ct = mco_co_mmpdist_core(unit_arrmco_readin,co_fcode,ctx_obj_ct);
    if ( msync( ctx_obj_ct, binsz*(sizeof(ctx_obj_ct_t)), MS_SYNC ) < 0 )
     err(errno,"mco_co_dist()::ctx_obj_ct msync failed");
    munmap(ctx_obj_ct,binsz*(sizeof(ctx_obj_ct_t)));
    close(fd);
   }
   free_unit_arrmco(unit_arrmco_readin);
  }
 }
 ref_seq_num = mco_dstat_readin.infile_num ;
  qry_seq_num = co_dstat_readin.infile_num ;
 kmerlen = co_dstat_readin.kmerlen;
  dim_reduct_len = co_dstat_readin.dim_rd_len;
 char distf[PATHLEN];
 sprintf(distf, "%s/distance.out", distout_dir);
 FILE *distfp;
 if( (distfp = fopen(distf,"a")) == NULL ) err(errno,"mco_co_dist():%s",distf);
 for(int i=0; i<=ref_bin_num;i++ ){
  for ( int k = 0; k < co_dstat_readin.infile_num; k++ ){
   if(qry_ctx_ct_list[k]>0)
    fname_dist_print(i,k,distout_dir,ref_ctx_ct_list,qry_ctx_ct_list,mcofname,cofname,distfp);
  }
 }
 fclose(distfp);
 free(ref_ctx_ct_list);
 free(qry_ctx_ct_list);
 free(mcofname);
 free(cofname);
}
void mco_cbd_co_dist(dist_opt_val_t *opt_val_in)
{
 int p_fit_mem = opt_val_in->p;
 llong mem_limit = (llong)opt_val_in->mmry*BBILLION;
 char *refmco_dname = opt_val_in->refpath;
 char *qryco_dname = opt_val_in->remaining_args[0];
 const char *distout_dir = opt_val_in->outdir;
  printf("run mco_cbd_co_dist(), %fG memory used\t%d threads used\n",opt_val_in->mmry,p_fit_mem);
  FILE *refmco_dstat_fp, *qryco_dstat_fp;
  char *refmco_dstat_fpath = malloc(PATHLEN*sizeof(char));
  char *qryco_dstat_fpath = malloc(PATHLEN*sizeof(char));
  sprintf(refmco_dstat_fpath,"%s/%s",refmco_dname,mco_dstat);
  sprintf(qryco_dstat_fpath,"%s/%s",qryco_dname,co_dstat);
  if( (refmco_dstat_fp = fopen(refmco_dstat_fpath,"rb")) == NULL )
    err(errno,"need provied mco dir path for mco_co_dist() arg 1. refmco_dstat_fpath");
  if( (qryco_dstat_fp = fopen(qryco_dstat_fpath,"rb")) == NULL )
    err(errno,"need provied co dir path for mco_co_dist() arg 2.  qryco_dstat_fpath");
  mco_dstat_t mco_dstat_readin ;
  co_dstat_t co_dstat_readin ;
  fread(&mco_dstat_readin,sizeof(mco_dstat_readin),1,refmco_dstat_fp);
  fread(&co_dstat_readin,sizeof(co_dstat_readin),1,qryco_dstat_fp);
  ctx_obj_ct_t * qry_ctx_ct_list = malloc(co_dstat_readin.infile_num * sizeof(ctx_obj_ct_t));
  ctx_obj_ct_t * ref_ctx_ct_list = malloc(mco_dstat_readin.infile_num * sizeof(ctx_obj_ct_t));
  fread(qry_ctx_ct_list,sizeof(ctx_obj_ct_t),co_dstat_readin.infile_num,qryco_dstat_fp);
  fread(ref_ctx_ct_list,sizeof(ctx_obj_ct_t),mco_dstat_readin.infile_num,refmco_dstat_fp);
  char (*cofname)[PATHLEN] = malloc(co_dstat_readin.infile_num * PATHLEN);
  char (*mcofname)[PATHLEN] = malloc(mco_dstat_readin.infile_num * PATHLEN);
  fread(cofname,PATHLEN,co_dstat_readin.infile_num,qryco_dstat_fp);
  fread(mcofname,PATHLEN,mco_dstat_readin.infile_num,refmco_dstat_fp);
  fclose(refmco_dstat_fp);
  fclose(qryco_dstat_fp);
  if( !(mco_dstat_readin.comp_num == co_dstat_readin.comp_num) )
    err(errno,"query args not match ref args: ref.comp_num = %d vs. %d = qry.comp_num",
    mco_dstat_readin.comp_num, co_dstat_readin.comp_num);
  if(!(mco_dstat_readin.shuf_id == co_dstat_readin.shuf_id))
    err(errno,"query args not match ref args: ref.shuf_id = %d vs. %d = qry.shuf_id",
    mco_dstat_readin.shuf_id, co_dstat_readin.shuf_id);
  int ref_bin_num = mco_dstat_readin.infile_num / BIN_SZ;
 if( mco_dstat_readin.infile_num % BIN_SZ > 0 ) ref_bin_num +=1;
 char onedist[PATHLEN];
 sprintf(onedist,"%s/sharedk_ct.dat",distout_dir);
  int dist_bfp = open(onedist,O_RDWR,0600) ;
  if (dist_bfp == -1) {
  close(dist_bfp);
  dist_bfp = open(onedist,O_RDWR|O_CREAT, 0600) ;
  if (dist_bfp == -1) err(errno,"mco_cbd_co_dist()::%s",onedist);
 }
 else{
  errno = EEXIST;
  err(errno,"mco_cbd_co_dist()::%s",onedist);
 }
 size_t disf_sz = (size_t)mco_dstat_readin.infile_num*co_dstat_readin.infile_num*sizeof(ctx_obj_ct_t) ;
 if(ftruncate(dist_bfp, disf_sz) == -1) err(errno,"mco_cbd_co_dist()::ftruncate");
 close(dist_bfp);
 dist_bfp = open(onedist,O_RDWR, 0600);
 if (dist_bfp == -1) err(errno,"mco_cbd_co_dist()::%s",onedist);
 ctx_obj_ct_t *ctx_obj_ct = mmap(NULL,disf_sz,PROT_READ | PROT_WRITE,MAP_SHARED ,dist_bfp,0);
 if(ctx_obj_ct == MAP_FAILED) err(errno,"ctx_obj_ct mmap error");
 close(dist_bfp);
 int page_sz = sysconf(_SC_PAGESIZE);
 size_t comp_sz = (1LLU << 4*COMPONENT_SZ);
 if( comp_sz % page_sz != 0 ) err(errno,"comp_sz %d is not multiple of page_sz %d ",comp_sz,page_sz );
 int num_unit_mem = mem_limit / (mco_dstat_readin.infile_num*sizeof(ctx_obj_ct_t) * page_sz);
 if(num_unit_mem < 1) err(errno,"at least %fG memory needed to map ./onedist, specify more memory use -m",
  (float)mco_dstat_readin.infile_num*sizeof(ctx_obj_ct_t) * page_sz/1073741824 );
 int num_cof_batch = num_unit_mem*page_sz;
 size_t unitsz_distf_mapped = (size_t)num_cof_batch * mco_dstat_readin.infile_num * sizeof(ctx_obj_ct_t) ;
 int num_mapping_distf = co_dstat_readin.infile_num / num_cof_batch ;
 size_t maplength;
 int bnum_infile;
 FILE *cbd_fcode_comp_fp,*cbd_fcode_comp_index_fp;
  struct stat cbd_fcode_stat;
 size_t *fco_pos = malloc(sizeof(size_t) * (co_dstat_readin.infile_num + 1) );
 size_t *mco_offset_index = malloc(sizeof(size_t) * comp_sz);
 unsigned int *mco_bin_index = malloc( sizeof(unsigned int) * comp_sz * ref_bin_num );
 gidobj_t* mco_mem = malloc( sizeof(gidobj_t) * 442317172 );
 char mco_fcode[PATHLEN]; char mco_index_fcode[PATHLEN];
 char co_cbd_fcode[PATHLEN];char co_cbd_index_fcode[PATHLEN];
 for(int b=0;b<=num_mapping_distf;b++){
  if(b==num_mapping_distf){
   bnum_infile = co_dstat_readin.infile_num % num_cof_batch ;
   if( bnum_infile == 0 ) continue;
  }else bnum_infile = num_cof_batch;
  maplength = (size_t)bnum_infile * mco_dstat_readin.infile_num * sizeof(ctx_obj_ct_t);
  printf("disf_sz=%lu\trefnum=%d\tqrynum=%d\tnum_mapping_distf=%dbnum_infile=%d\t\t%lu\t%lu\tflag1:Ok\n",
   disf_sz,mco_dstat_readin.infile_num,co_dstat_readin.infile_num,num_mapping_distf,bnum_infile,maplength,(size_t)b*unitsz_distf_mapped);
    for ( int j = 0; j < mco_dstat_readin.comp_num; j++ ) {
   sprintf(mco_index_fcode,"%s/mco.index.%d",refmco_dname,j);
      sprintf(mco_fcode,"%s/mco.%d",refmco_dname,j);
   FILE *indexfp, *mcofp;
   if( (indexfp = fopen(mco_index_fcode,"rb"))==NULL) err(errno,"mco_cbd_co_dist()::%s",mco_index_fcode);
   fread(mco_offset_index,sizeof(size_t),comp_sz,indexfp);
   fread(mco_bin_index,sizeof(unsigned int),comp_sz * ref_bin_num,indexfp);
   fclose(indexfp);
   struct stat s;
   if( (mcofp = fopen(mco_fcode,"rb"))==NULL) err(errno,"mco_cbd_co_dist()::%s",mco_fcode);
   stat(mco_fcode, &s);
   fread(mco_mem,sizeof(gidobj_t),s.st_size/sizeof(gidobj_t),mcofp);
   fclose(mcofp);
      sprintf(co_cbd_fcode,"%s/combco.%d",qryco_dname,j);
      if( (cbd_fcode_comp_fp = fopen(co_cbd_fcode,"rb"))==NULL) err(errno,"mco_cbd_co_dist()::%s",co_cbd_fcode);
      stat(co_cbd_fcode, &cbd_fcode_stat);
   if(co_dstat_readin.koc){
   }
      unsigned int *cbd_fcode_mem = malloc(cbd_fcode_stat.st_size);
      fread(cbd_fcode_mem,sizeof(unsigned int),cbd_fcode_stat.st_size/sizeof(unsigned int),cbd_fcode_comp_fp);
      fclose(cbd_fcode_comp_fp);
   sprintf(co_cbd_index_fcode,"%s/combco.index.%d",qryco_dname,j);
   if( (cbd_fcode_comp_index_fp = fopen(co_cbd_index_fcode,"rb"))==NULL)
    err(errno,"mco_cbd_co_dist()::%s",co_cbd_index_fcode);
   fread(fco_pos,sizeof(size_t),co_dstat_readin.infile_num + 1 ,cbd_fcode_comp_index_fp);
   fclose(cbd_fcode_comp_index_fp);
#pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
      for(int kind = 0; kind < bnum_infile; kind++){
    int k = b*num_cof_batch + kind;
        if(qry_ctx_ct_list[k]==0) continue;
    unsigned int ind, mcogid, pos;
    llong distf_offset = (size_t)k * mco_dstat_readin.infile_num;
     for(int n = 0; n < fco_pos[k+1] - fco_pos[k]; n++){
       ind = cbd_fcode_mem[ fco_pos[k] + n ];
     pos = 0;
     for(int bin=0; bin < ref_bin_num; bin++){
      int bin_gnum = mco_bin_index[ ind *ref_bin_num + bin];
      for(int g = 0; g < bin_gnum ; g++ ){
       mcogid = bin*BIN_SZ + mco_mem[ mco_offset_index[ind] + pos ];
       ctx_obj_ct[distf_offset + mcogid ]++ ;
       pos++;
      }
     }
     }
     }
    free(cbd_fcode_mem);
    munmap(mco_offset_index, comp_sz * sizeof(size_t));
    munmap(mco_bin_index,comp_sz*ref_bin_num*sizeof(unsigned int) );
   }
    if ( msync( ctx_obj_ct + (size_t)b*num_cof_batch*mco_dstat_readin.infile_num, maplength, MS_SYNC ) < 0 )
     err(errno,"mco_cbd_co_dist()::ctx_obj_ct msync failed");
    munmap(ctx_obj_ct + (size_t)b*num_cof_batch*mco_dstat_readin.infile_num, maplength);
  }
 free(fco_pos);
  ref_seq_num = mco_dstat_readin.infile_num ;
  qry_seq_num = co_dstat_readin.infile_num ;
  kmerlen = co_dstat_readin.kmerlen;
  dim_reduct_len = co_dstat_readin.dim_rd_len;
 dist_print_nobin(distout_dir,ref_seq_num, qry_seq_num, ref_ctx_ct_list, qry_ctx_ct_list,num_cof_batch,mcofname, cofname, opt_val_in);
  free(ref_ctx_ct_list);
  free(qry_ctx_ct_list);
  free(mcofname);
  free(cofname);
}
void mco_cbdco_nobin_dist(dist_opt_val_t *opt_val_in)
{
  int p_fit_mem = opt_val_in->p;
  llong mem_limit = (llong)opt_val_in->mmry*BBILLION;
  char *refmco_dname = opt_val_in->refpath;
  char *qryco_dname = opt_val_in->remaining_args[0];
  const char *distout_dir = opt_val_in->outdir;
  mkdir(distout_dir,0700);
  FILE *refmco_dstat_fp, *qryco_dstat_fp;
  char *refmco_dstat_fpath = malloc(PATHLEN*sizeof(char));
  char *qryco_dstat_fpath = malloc(PATHLEN*sizeof(char));
  sprintf(refmco_dstat_fpath,"%s/%s",refmco_dname,mco_dstat);
  sprintf(qryco_dstat_fpath,"%s/%s",qryco_dname,co_dstat);
  if( (refmco_dstat_fp = fopen(refmco_dstat_fpath,"rb")) == NULL )
    err(errno,"need provied mco dir path for mco_co_dist() arg 1. refmco_dstat_fpath");
  if( (qryco_dstat_fp = fopen(qryco_dstat_fpath,"rb")) == NULL )
    err(errno,"need provied co dir path for mco_co_dist() arg 2.  qryco_dstat_fpath");
  mco_dstat_t mco_dstat_readin ;
  co_dstat_t co_dstat_readin ;
  fread(&mco_dstat_readin,sizeof(mco_dstat_readin),1,refmco_dstat_fp);
  fread(&co_dstat_readin,sizeof(co_dstat_readin),1,qryco_dstat_fp);
  ctx_obj_ct_t * qry_ctx_ct_list = malloc(co_dstat_readin.infile_num * sizeof(ctx_obj_ct_t));
  ctx_obj_ct_t * ref_ctx_ct_list = malloc(mco_dstat_readin.infile_num * sizeof(ctx_obj_ct_t));
  fread(qry_ctx_ct_list,sizeof(ctx_obj_ct_t),co_dstat_readin.infile_num,qryco_dstat_fp);
  fread(ref_ctx_ct_list,sizeof(ctx_obj_ct_t),mco_dstat_readin.infile_num,refmco_dstat_fp);
  char (*cofname)[PATHLEN] = malloc(co_dstat_readin.infile_num * PATHLEN);
  char (*mcofname)[PATHLEN] = malloc(mco_dstat_readin.infile_num * PATHLEN);
  fread(cofname,PATHLEN,co_dstat_readin.infile_num,qryco_dstat_fp);
  fread(mcofname,PATHLEN,mco_dstat_readin.infile_num,refmco_dstat_fp);
  fclose(refmco_dstat_fp);
  fclose(qryco_dstat_fp);
  if( !(mco_dstat_readin.comp_num == co_dstat_readin.comp_num) )
    err(errno,"query args not match ref args: ref.comp_num = %d vs. %d = qry.comp_num",
    mco_dstat_readin.comp_num, co_dstat_readin.comp_num);
  if(!(mco_dstat_readin.shuf_id == co_dstat_readin.shuf_id))
    err(errno,"query args not match ref args: ref.shuf_id = %d vs. %d = qry.shuf_id",
    mco_dstat_readin.shuf_id, co_dstat_readin.shuf_id);
  char onedist[PATHLEN];
  sprintf(onedist,"%s/sharedk_ct.dat",distout_dir);
  int dist_bfp = open(onedist,O_RDWR,0600) ;
  if (dist_bfp == -1) {
    close(dist_bfp);
    dist_bfp = open(onedist,O_RDWR|O_CREAT, 0600) ;
    if (dist_bfp == -1) err(errno," mco_cbdco_nobin_dist()::%s",onedist);
  }else err(EEXIST," mco_cbdco_nobin_dist():%s",onedist);
  int page_sz = sysconf(_SC_PAGESIZE);
  size_t comp_sz = (1LLU << 4*COMPONENT_SZ);
  if( comp_sz % page_sz != 0 ) err(errno,"comp_sz %d is not multiple of page_sz %d ",comp_sz,page_sz );
  size_t maplength;
  int bnum_infile;
  size_t *mco_offset_index = malloc(sizeof(size_t) * comp_sz);
#define MCOMM_SZ 442317172
  gidobj_t* mco_mem = malloc( sizeof(gidobj_t) * MCOMM_SZ );
  char mco_fcode[PATHLEN]; char mco_index_fcode[PATHLEN];
  char co_cbd_fcode[PATHLEN];char co_cbd_index_fcode[PATHLEN];
  ref_seq_num = mco_dstat_readin.infile_num ;
  qry_seq_num = co_dstat_readin.infile_num ;
  kmerlen = co_dstat_readin.kmerlen;
  dim_reduct_len = co_dstat_readin.dim_rd_len;
  char distf[PATHLEN];
  sprintf(distf, "%s/distance.out", distout_dir);
  int num_unit_mem = mem_limit / (mco_dstat_readin.infile_num*sizeof(ctx_obj_ct_t) * page_sz);
  if(num_unit_mem < 1) err(errno,"at least %fG memory needed to map ./onedist, specify more memory use -m",
      (float)mco_dstat_readin.infile_num*sizeof(ctx_obj_ct_t) * page_sz/1073741824 );
  int num_cof_batch = num_unit_mem*page_sz;
  if( opt_val_in->shared_kmerpath[0] != '\0'){
    dist_print_nobin(distout_dir,ref_seq_num, qry_seq_num, ref_ctx_ct_list, qry_ctx_ct_list,num_cof_batch,mcofname, cofname,opt_val_in);
    return;
  }
  size_t unitsz_distf_mapped = (size_t)num_cof_batch * mco_dstat_readin.infile_num * sizeof(ctx_obj_ct_t) ;
  int num_mapping_distf = co_dstat_readin.infile_num / num_cof_batch ;
  size_t disf_sz = (size_t)mco_dstat_readin.infile_num*co_dstat_readin.infile_num*sizeof(ctx_obj_ct_t) ;
  if(ftruncate(dist_bfp, disf_sz) == -1) err(errno,"mco_cbdco_nobin_dist()::ftruncate");
  close(dist_bfp);
  dist_bfp = open(onedist,O_RDWR, 0600);
  if (dist_bfp == -1) err(errno,"mco_cbdco_nobin_dist()::%s",onedist);
  ctx_obj_ct_t *ctx_obj_ct = mmap(NULL,disf_sz,PROT_READ | PROT_WRITE,MAP_SHARED,dist_bfp,0);
  if(ctx_obj_ct == MAP_FAILED) err(errno,"ctx_obj_ct mmap error");
  close(dist_bfp);
 mmp_any_t* mco_offset_index_comps = malloc(sizeof( mmp_any_t) * mco_dstat_readin.comp_num);
 mmp_any_t* mco_mem_comps = malloc(sizeof( mmp_any_t) * mco_dstat_readin.comp_num);
 mmp_any_t* fco_pos_comps = malloc(sizeof( mmp_any_t) * mco_dstat_readin.comp_num);
 mmp_any_t* cbd_fcode_mem_comps = malloc(sizeof( mmp_any_t) * mco_dstat_readin.comp_num);
 for ( int j = 0; j < mco_dstat_readin.comp_num; j++ ) {
  sprintf(mco_index_fcode,"%s/%s.%d",refmco_dname,mco_idx_prefix,j);
  mco_offset_index_comps[j] = mmp_any(mco_index_fcode);
  sprintf(mco_fcode,"%s/%s.%d",refmco_dname,mco_gids_prefix,j);
  mco_mem_comps[j] = mmp_any(mco_fcode);
  sprintf(co_cbd_index_fcode,"%s/%s.%d",qryco_dname,idx_prefix,j);
  fco_pos_comps[j] = mmp_any(co_cbd_index_fcode);
  sprintf(co_cbd_fcode,"%s/%s.%d",qryco_dname,skch_prefix,j);
  cbd_fcode_mem_comps[j] = mmp_any(co_cbd_fcode);
 }
  for(int b=0;b<=num_mapping_distf;b++){
    if(b==num_mapping_distf) bnum_infile = co_dstat_readin.infile_num % num_cof_batch ;
  else bnum_infile = num_cof_batch;
    if( bnum_infile == 0 ) continue;
    maplength = (size_t)bnum_infile * mco_dstat_readin.infile_num * sizeof(ctx_obj_ct_t);
printf("disf_sz=%lu\trefnum=%d\tqrynum=%d\tnum_mapping_distf=%d\tbatch_qrynum=%d\t%lu\t%lu\n",disf_sz,mco_dstat_readin.infile_num,co_dstat_readin.infile_num,num_mapping_distf,bnum_infile,maplength,(size_t)b*unitsz_distf_mapped);
    for ( int j = 0; j < mco_dstat_readin.comp_num; j++ ) {
   size_t *mco_offset_index= (size_t *) mco_offset_index_comps[j].mmp;
   gidobj_t* mco_mem = (gidobj_t*) mco_mem_comps[j].mmp;
   size_t *fco_pos = (size_t *)fco_pos_comps[j].mmp;
   unsigned int *cbd_fcode_mem = (unsigned int *)cbd_fcode_mem_comps[j].mmp;
#pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
      for(int kind = 0; kind < bnum_infile; kind++){
        int k = b*num_cof_batch + kind;
        if(qry_ctx_ct_list[k]==0) continue;
        size_t distf_offset = (size_t)k * mco_dstat_readin.infile_num;
        for(size_t n = fco_pos[k]; n < fco_pos[k+1]; n++){
          unsigned int ind = cbd_fcode_mem[n];
     size_t s = ind > 0 ? mco_offset_index[ind-1] : 0 ;
     for(size_t g = s; g < mco_offset_index[ind] ; g++)
      ctx_obj_ct[distf_offset + mco_mem[g]]++ ;
        }
      }
    }
    if ( msync( ctx_obj_ct + (size_t)b*num_cof_batch*mco_dstat_readin.infile_num, maplength, MS_SYNC ) < 0 )
      err(errno,"mco_cbdco_nobin_dist()::ctx_obj_ct msync failed");
    munmap(ctx_obj_ct + (size_t)b*num_cof_batch*mco_dstat_readin.infile_num, maplength);
  }
 for ( int j = 0; j < mco_dstat_readin.comp_num; j++ ){
  munmap(mco_offset_index_comps[j].mmp,mco_offset_index_comps[j].fsize);
  munmap(mco_mem_comps[j].mmp,mco_mem_comps[j].fsize);
  munmap(fco_pos_comps[j].mmp,fco_pos_comps[j].fsize);
  munmap(cbd_fcode_mem_comps[j].mmp,cbd_fcode_mem_comps[j].fsize);
 }
  free(mco_offset_index_comps);
  free( mco_mem_comps);
  free( fco_pos_comps);
  free (cbd_fcode_mem_comps);
  dist_print_nobin(distout_dir,ref_seq_num, qry_seq_num, ref_ctx_ct_list, qry_ctx_ct_list,num_cof_batch,mcofname,cofname,opt_val_in);
  free(mco_offset_index);
  free(mco_mem);
  free(ref_ctx_ct_list);
  free(qry_ctx_ct_list);
  free(mcofname);
  free(cofname);
}
typedef struct koc_dist {
  llong shared_koc_ct;
  ctx_obj_ct_t shared_k_ct;
} koc_dist_t;
void mco_cbd_koc_compatible_dist(dist_opt_val_t *opt_val_in)
{
  int p_fit_mem = opt_val_in->p;
  llong mem_limit = (llong)opt_val_in->mmry*BBILLION;
  char *refmco_dname = opt_val_in->refpath;
  char *qryco_dname = opt_val_in->remaining_args[0];
  const char *distout_dir = opt_val_in->outdir;
 mkdir(distout_dir,0700);
  FILE *refmco_dstat_fp, *qryco_dstat_fp;
  char *refmco_dstat_fpath = malloc(PATHLEN*sizeof(char));
  char *qryco_dstat_fpath = malloc(PATHLEN*sizeof(char));
  sprintf(refmco_dstat_fpath,"%s/%s",refmco_dname,mco_dstat);
  sprintf(qryco_dstat_fpath,"%s/%s",qryco_dname,co_dstat);
  if( (refmco_dstat_fp = fopen(refmco_dstat_fpath,"rb")) == NULL )
    err(errno,"need provied mco dir path for mco_co_dist() arg 1. refmco_dstat_fpath");
  if( (qryco_dstat_fp = fopen(qryco_dstat_fpath,"rb")) == NULL )
    err(errno,"need provied co dir path for mco_co_dist() arg 2.  qryco_dstat_fpath");
  mco_dstat_t mco_dstat_readin ;
  co_dstat_t co_dstat_readin ;
  fread(&mco_dstat_readin,sizeof(mco_dstat_readin),1,refmco_dstat_fp);
  fread(&co_dstat_readin,sizeof(co_dstat_readin),1,qryco_dstat_fp);
  ctx_obj_ct_t * qry_ctx_ct_list = malloc(co_dstat_readin.infile_num * sizeof(ctx_obj_ct_t));
  ctx_obj_ct_t * ref_ctx_ct_list = malloc(mco_dstat_readin.infile_num * sizeof(ctx_obj_ct_t));
  fread(qry_ctx_ct_list,sizeof(ctx_obj_ct_t),co_dstat_readin.infile_num,qryco_dstat_fp);
  fread(ref_ctx_ct_list,sizeof(ctx_obj_ct_t),mco_dstat_readin.infile_num,refmco_dstat_fp);
  char (*cofname)[PATHLEN] = malloc(co_dstat_readin.infile_num * PATHLEN);
  char (*mcofname)[PATHLEN] = malloc(mco_dstat_readin.infile_num * PATHLEN);
  fread(cofname,PATHLEN,co_dstat_readin.infile_num,qryco_dstat_fp);
  fread(mcofname,PATHLEN,mco_dstat_readin.infile_num,refmco_dstat_fp);
  fclose(refmco_dstat_fp);
  fclose(qryco_dstat_fp);
  if( !(mco_dstat_readin.comp_num == co_dstat_readin.comp_num) )
    err(errno,"query args not match ref args: ref.comp_num = %d vs. %d = qry.comp_num",
    mco_dstat_readin.comp_num, co_dstat_readin.comp_num);
  if(!(mco_dstat_readin.shuf_id == co_dstat_readin.shuf_id))
    err(errno,"query args not match ref args: ref.shuf_id = %d vs. %d = qry.shuf_id",
    mco_dstat_readin.shuf_id, co_dstat_readin.shuf_id);
  int ref_bin_num = mco_dstat_readin.infile_num / BIN_SZ;
  if( mco_dstat_readin.infile_num % BIN_SZ > 0 ) ref_bin_num +=1;
  char onedist[PATHLEN];
  sprintf(onedist,"%s/sharedk_ct.dat",distout_dir);
  int dist_bfp = open(onedist,O_RDWR,0600) ;
  if (dist_bfp == -1) {
    close(dist_bfp);
    dist_bfp = open(onedist,O_RDWR|O_CREAT, 0600) ;
    if (dist_bfp == -1) err(errno,"mco_cbd_co_dist()::%s",onedist);
  }else err(EEXIST,"mco_cbd_koc_compatible_dist:%s",onedist);
 int page_sz = sysconf(_SC_PAGESIZE);
  size_t comp_sz = (1LLU << 4*COMPONENT_SZ);
  if( comp_sz % page_sz != 0 ) err(errno,"comp_sz %d is not multiple of page_sz %d ",comp_sz,page_sz );
 size_t maplength;
  int bnum_infile;
 FILE *cbd_fcode_comp_fp,*cbd_fcode_comp_index_fp;
  struct stat cbd_fcode_stat;
  size_t *fco_pos = malloc(sizeof(size_t) * (co_dstat_readin.infile_num + 1) );
  size_t *mco_offset_index = malloc(sizeof(size_t) * comp_sz);
  unsigned int *mco_bin_index = malloc( (size_t)comp_sz * sizeof(unsigned int) * ref_bin_num );
#define MCOMM_SZ 442317172
  gidobj_t* mco_mem = malloc( sizeof(gidobj_t) * MCOMM_SZ );
  char mco_fcode[PATHLEN]; char mco_index_fcode[PATHLEN];
  char co_cbd_fcode[PATHLEN];char co_cbd_index_fcode[PATHLEN];
 ref_seq_num = mco_dstat_readin.infile_num ;
  qry_seq_num = co_dstat_readin.infile_num ;
  kmerlen = co_dstat_readin.kmerlen;
  dim_reduct_len = co_dstat_readin.dim_rd_len;
  char distf[PATHLEN];
  sprintf(distf, "%s/distance.out", distout_dir);
 co_dstat_readin.koc = 0;
  int num_unit_mem = mem_limit / (mco_dstat_readin.infile_num*sizeof(ctx_obj_ct_t) * page_sz);
  if(num_unit_mem < 1) err(errno,"at least %fG memory needed to map ./onedist, specify more memory use -m",
      (float)mco_dstat_readin.infile_num*sizeof(ctx_obj_ct_t) * page_sz/1073741824 );
  int num_cof_batch = num_unit_mem*page_sz;
 if( opt_val_in->shared_kmerpath[0] != '\0'){
  dist_print_nobin(distout_dir,ref_seq_num, qry_seq_num, ref_ctx_ct_list, qry_ctx_ct_list,num_cof_batch,mcofname, cofname,opt_val_in);
  return;
 }
  size_t unitsz_distf_mapped = (size_t)num_cof_batch * mco_dstat_readin.infile_num * sizeof(ctx_obj_ct_t) ;
  int num_mapping_distf = co_dstat_readin.infile_num / num_cof_batch ;
 size_t disf_sz = (size_t)mco_dstat_readin.infile_num*co_dstat_readin.infile_num*sizeof(ctx_obj_ct_t) ;
  if(ftruncate(dist_bfp, disf_sz) == -1) err(errno,"mco_cbd_koc_compatible_dist()::ftruncate");
  close(dist_bfp);
 dist_bfp = open(onedist,O_RDWR, 0600);
  if (dist_bfp == -1) err(errno,"mco_cbd_koc_compatible_dist()::%s",onedist);
  ctx_obj_ct_t *ctx_obj_ct = mmap(NULL,disf_sz,PROT_READ | PROT_WRITE,MAP_SHARED,dist_bfp,0);
  if(ctx_obj_ct == MAP_FAILED) err(errno,"ctx_obj_ct mmap error");
  close(dist_bfp);
  for(int b=0;b<=num_mapping_distf;b++){
    if(b==num_mapping_distf){
      bnum_infile = co_dstat_readin.infile_num % num_cof_batch ;
      if( bnum_infile == 0 ) continue;
    }else bnum_infile = num_cof_batch;
    maplength = (size_t)bnum_infile * mco_dstat_readin.infile_num * sizeof(ctx_obj_ct_t);
    printf("disf_sz=%lu\trefnum=%d\tqrynum=%d\tnum_mapping_distf=%dbnum_infile=%d\t\t%lu\t%lu\tflag1:Ok\n",
      disf_sz,mco_dstat_readin.infile_num,co_dstat_readin.infile_num,num_mapping_distf,bnum_infile,maplength,(size_t)b*unitsz_distf_mapped);
    for ( int j = 0; j < mco_dstat_readin.comp_num; j++ ) {
      sprintf(mco_index_fcode,"%s/mco.index.%d",refmco_dname,j);
      sprintf(mco_fcode,"%s/mco.%d",refmco_dname,j);
      FILE *indexfp, *mcofp;
      if( (indexfp = fopen(mco_index_fcode,"rb"))==NULL) err(errno,"mco_cbd_koc_compatible_dist()::%s",mco_index_fcode);
      fread(mco_offset_index,sizeof(size_t),comp_sz,indexfp);
      fread(mco_bin_index,sizeof(unsigned int),(size_t)comp_sz * ref_bin_num,indexfp);
      fclose(indexfp);
      struct stat s;
      if( (mcofp = fopen(mco_fcode,"rb"))==NULL) err(errno,"mco_cbd_koc_compatible_dist()::%s",mco_fcode);
      stat(mco_fcode, &s);
   if(s.st_size > sizeof(gidobj_t) * MCOMM_SZ) mco_mem = (gidobj_t *)realloc(mco_mem, s.st_size );
      fread(mco_mem,sizeof(gidobj_t),s.st_size/sizeof(gidobj_t),mcofp);
      fclose(mcofp);
      sprintf(co_cbd_fcode,"%s/combco.%d",qryco_dname,j);
      if( (cbd_fcode_comp_fp = fopen(co_cbd_fcode,"rb"))==NULL) err(errno,"mco_cbd_koc_compatible_dis()::%s",co_cbd_fcode);
      stat(co_cbd_fcode, &cbd_fcode_stat);
      unsigned int *cbd_fcode_mem = malloc(cbd_fcode_stat.st_size);
      fread(cbd_fcode_mem,sizeof(unsigned int),cbd_fcode_stat.st_size/sizeof(unsigned int),cbd_fcode_comp_fp);
      fclose(cbd_fcode_comp_fp);
      sprintf(co_cbd_index_fcode,"%s/combco.index.%d",qryco_dname,j);
      if( (cbd_fcode_comp_index_fp = fopen(co_cbd_index_fcode,"rb"))==NULL)
        err(errno,"mco_cbd_koc_compatible_dis()::%s",co_cbd_index_fcode);
      fread(fco_pos,sizeof(size_t),co_dstat_readin.infile_num + 1, cbd_fcode_comp_index_fp);
      fclose(cbd_fcode_comp_index_fp);
#pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
      for(int kind = 0; kind < bnum_infile; kind++){
        int k = b*num_cof_batch + kind;
        if(qry_ctx_ct_list[k]==0) continue;
        size_t ind, mcogid,pos;
        size_t distf_offset = (size_t)k * mco_dstat_readin.infile_num;
        for(size_t n = 0; n < fco_pos[k+1] - fco_pos[k]; n++){
          ind = cbd_fcode_mem[ fco_pos[k] + n ];
          pos = 0;
          for(int bin=0; bin < ref_bin_num; bin++){
            int bin_gnum = mco_bin_index[ (size_t)ind *ref_bin_num + bin ];
            for(int g = 0; g < bin_gnum ; g++ ){
              mcogid = bin*BIN_SZ + mco_mem[ (size_t)mco_offset_index[ind] + pos ];
              ctx_obj_ct[distf_offset + mcogid]++ ;
              pos++;
            }
          }
        }
      }
      free(cbd_fcode_mem);
    }
    if ( msync( ctx_obj_ct + (size_t)b*num_cof_batch*mco_dstat_readin.infile_num, maplength, MS_SYNC ) < 0 )
      err(errno,"mco_cbd_co_dist()::ctx_obj_ct msync failed");
    munmap(ctx_obj_ct + (size_t)b*num_cof_batch*mco_dstat_readin.infile_num, maplength);
  }
    free(fco_pos);
  dist_print_nobin(distout_dir,ref_seq_num, qry_seq_num, ref_ctx_ct_list, qry_ctx_ct_list,num_cof_batch,mcofname, cofname,opt_val_in);
 free(mco_offset_index);
  free(mco_bin_index);
 free(mco_mem);
  free(ref_ctx_ct_list);
  free(qry_ctx_ct_list);
  free(mcofname);
  free(cofname);
}
static inline ctx_obj_ct_t * mco_co_mmpdist_core(gidobj_t** unit_arrmco, char *co_fcode_in, ctx_obj_ct_t * ctx_obj_ct_in )
{
  mmp_uint_t mmpcofile;
  mmpcofile = mmp_uint_arr(co_fcode_in);
  unsigned int ind,mcogid;
  int ctx_num = mmpcofile.fsize/sizeof(unsigned int);
 for(int n = 0; n < ctx_num; n++){
    ind = mmpcofile.mmpco[n];
    if(unit_arrmco[ind] != NULL){
      for(unsigned int k = 1; k< unit_arrmco[ind][0] + 1; k++ ){
        mcogid = unit_arrmco[ind][k] ;
        ctx_obj_ct_in[mcogid]++ ;
     }
    }
  }
  munmap(mmpcofile.mmpco, mmpcofile.fsize);
 return ctx_obj_ct_in;
}
static inline void mco_co_dist_core(gidobj_t** unit_arrmco, char *co_fcode_in, int bin_sz,
      mco_co_dist_t shared_ctx_num_in )
{
 mmp_uint_t mmpcofile;
 mmpcofile = mmp_uint_arr(co_fcode_in);
 memset(shared_ctx_num_in,0, bin_sz);
 unsigned int ind,mcogid;
 int ctx_num = mmpcofile.fsize/sizeof(unsigned int);
 for(int n = 0; n < ctx_num; n++){
  ind = mmpcofile.mmpco[n] ;
  if(unit_arrmco[ind] != NULL){
   for(int k = 1; k< unit_arrmco[ind][0] + 1; k++ ){
    mcogid = unit_arrmco[ind][k];
     shared_ctx_num_in[mcogid]++ ;
   }
  }
 }
 munmap(mmpcofile.mmpco, mmpcofile.fsize);
 printf("%s:",co_fcode_in);
 for(int e=0; e < bin_sz; e++ ){
 }
 printf("\n");
}
void dist_print( const char *distf, FILE *dist_fp )
{
 ctx_obj_ct_t *ctx_obj_ct;
 int fd;
 struct stat s;
 fd = open (distf, O_RDONLY);
 check (fd < 0, "open %s failed: %s", distf, strerror (errno));
 fstat (fd, & s);
 ctx_obj_ct = mmap(0, s.st_size , PROT_READ, MAP_PRIVATE, fd, 0);
 check ( ctx_obj_ct == MAP_FAILED, "mmap %s failed: %s", distf, strerror (errno));
  fprintf(dist_fp,"output %s\n",distf);
 close(fd);
 munmap(ctx_obj_ct, s.st_size);
}
char full_distfcode[PATHLEN];
void fname_dist_print(int ref_bin_code, int qry_fcode, const char *distout_dir, unsigned int*ref_ctx_ct_list,
   unsigned int*qry_ctx_ct_list, char (*refname)[PATHLEN], char (*qryfname)[PATHLEN], FILE *dout_fp)
{
 sprintf(full_distfcode,"%s/%d.%d.dist",distout_dir,ref_bin_code,qry_fcode);
  ctx_obj_ct_t *ctx_obj_ct;
 int fd;
 struct stat s;
 fd = open (full_distfcode, O_RDONLY);
 check (fd < 0, "open %s failed: %s", full_distfcode, strerror (errno));
 fstat (fd, & s);
 ctx_obj_ct = mmap(0, s.st_size , PROT_READ, MAP_PRIVATE, fd, 0);
 check ( ctx_obj_ct == MAP_FAILED, "mmap %s failed: %s", full_distfcode, strerror (errno));
 double jac_ind,contain_ind,Dm,Da,P_K_in_X_XnY, P_K_in_Y_XnY,
 j_prim, c_prim,
 sd_j_prim, sd_c_prim,
 CI95_j_prim1, CI95_j_prim2, CI95_c_prim1, CI95_c_prim2,
 CI95_Dm_prim1,CI95_Dm_prim2,CI95_Da_prim1,CI95_Da_prim2;
 int Min_XY_size, X_size, Y_size, XnY_size, XuY_size, X_XnY_size, Y_XnY_size ;
 for(llong i = 0;i < s.st_size/sizeof(ctx_obj_ct_t); i++) {
  X_size = ref_ctx_ct_list[ref_bin_code*BIN_SZ + i];
  Y_size = qry_ctx_ct_list[qry_fcode];
  Min_XY_size = X_size < Y_size ? X_size : Y_size ;
  XnY_size = ctx_obj_ct[i];
  XuY_size = X_size + Y_size - XnY_size ;
  X_XnY_size = X_size - XnY_size;
  Y_XnY_size = Y_size- XnY_size;
  jac_ind = (double)XnY_size / XuY_size;
  contain_ind = (double)XnY_size / Min_XY_size ;
  Dm = jac_ind ==1? 0: -log(2*jac_ind/(1+jac_ind)) / kmerlen ;
  Da = contain_ind==1? 0: -log(contain_ind) / kmerlen ;
  P_K_in_X_XnY = 1 - pow( (1- 1/pow(alp_size,(kmerlen - dim_reduct_len) )), X_XnY_size );
  P_K_in_Y_XnY = 1 - pow( (1- 1/pow(alp_size,(kmerlen - dim_reduct_len) )), Y_XnY_size );
  double rs = P_K_in_X_XnY * P_K_in_Y_XnY * ( X_XnY_size + Y_XnY_size )
       /(P_K_in_X_XnY + P_K_in_Y_XnY - 2*P_K_in_X_XnY * P_K_in_Y_XnY);
  j_prim = ((double)XnY_size - rs) / XuY_size ;
  c_prim = ((double)XnY_size - rs) / Min_XY_size ;
  sd_j_prim = pow(j_prim*(1 - j_prim) / XuY_size, 0.5) ;
  sd_c_prim = pow(c_prim*(1 - c_prim) / Min_XY_size,0.5) ;
  CI95_j_prim1 = j_prim - 1.96*sd_j_prim;
  CI95_j_prim2 = j_prim + 1.96*sd_j_prim;
  CI95_c_prim1 = c_prim - 1.96*sd_c_prim;
  CI95_c_prim2 = c_prim + 1.96*sd_c_prim;
  CI95_Dm_prim1 = CI95_j_prim2 == 1? 0:-log(2*CI95_j_prim2/(1+CI95_j_prim2)) / kmerlen ;
  CI95_Dm_prim2 = CI95_j_prim1 == 1? 0:-log(2*CI95_j_prim1/(1+CI95_j_prim1)) / kmerlen ;
  CI95_Da_prim1 = CI95_c_prim2 == 1? 0:-log(CI95_c_prim2) / kmerlen ;
  CI95_Da_prim2 = CI95_c_prim1 == 1? 0:-log(CI95_c_prim1) / kmerlen ;
  double pv_j_prim = 0.5 * erfc( j_prim / sd_j_prim * pow(0.5,0.5) );
  double pv_c_prim = 0.5 * erfc( c_prim/sd_c_prim * pow(0.5,0.5) );
  double qv_j_prim = pv_j_prim * ref_seq_num*qry_seq_num ;
  double qv_c_prim = pv_c_prim * ref_seq_num*qry_seq_num ;
  fprintf(dout_fp,"%s\t%s\t%u-%u|%u|%u\t%lf\t%lf\t%lf\t%lf\t[%lf,%lf]\t[%lf,%lf]\t[%lf,%lf]\t[%lf,%lf]\t%E\t%E\t%E\t%E\n", qryfname[qry_fcode],refname[ref_bin_code*BIN_SZ + i],XnY_size,(unsigned int)rs,X_size,Y_size,jac_ind,Dm,
        contain_ind,Da,CI95_j_prim1,CI95_j_prim2,CI95_Dm_prim1,CI95_Dm_prim2,CI95_c_prim1,
    CI95_c_prim2,CI95_Da_prim1,CI95_Da_prim2,pv_j_prim,pv_c_prim,qv_j_prim,qv_c_prim);
 }
 close(fd);
  munmap(ctx_obj_ct, s.st_size);
}
void koc_dist_print_nobin ( const char *distout_dir,unsigned int ref_num, unsigned int qry_num, unsigned int*ref_ctx_ct_list,
      unsigned int*qry_ctx_ct_list, int num_cof_batch, char (*refname)[PATHLEN], char (*qryfname)[PATHLEN])
{
  sprintf(full_distfcode,"%s/sharedk_ct.dat",distout_dir);
  int fd;
  struct stat s;
  fd = open (full_distfcode, O_RDONLY);
  check (fd < 0, "open %s failed: %s", full_distfcode, strerror (errno));
  fstat (fd, & s);
  koc_dist_t * ctx_obj_ct = mmap(0, s.st_size , PROT_READ, MAP_PRIVATE, fd, 0);
  check ( ctx_obj_ct == MAP_FAILED, "mmap %s failed: %s", full_distfcode, strerror (errno));
  close(fd);
  char distf[PATHLEN];
  sprintf(distf, "%s/distance.out", distout_dir);
  FILE *distfp = fopen(distf,"a") ;
  if( distfp == NULL ) err(errno,"dist_print_nobin():%s",distf);
  double jac_ind,contain_ind,Dm,Da,P_K_in_X_XnY, P_K_in_Y_XnY,
  j_prim, c_prim, Dm_prim, Da_prim, sd_j_prim, sd_c_prim,
  CI95_j_prim1, CI95_j_prim2, CI95_c_prim1, CI95_c_prim2,
  CI95_Dm_prim1,CI95_Dm_prim2,CI95_Da_prim1,CI95_Da_prim2;
  int Min_XY_size, X_size, Y_size, XnY_size, XuY_size, X_XnY_size, Y_XnY_size ;
  int num_mapping_distf = qry_num / num_cof_batch; int bnum_infile;
 double abundence_rowsum;
  for(int b=0;b<=num_mapping_distf;b++){
    if(b==num_mapping_distf){
      bnum_infile = qry_num % num_cof_batch ;
      if( bnum_infile == 0 ) continue;
    }else bnum_infile = num_cof_batch;
    size_t maplength = (size_t)bnum_infile * ref_num * sizeof(koc_dist_t);
    for(int qid = 0; qid < bnum_infile; qid++) {
   abundence_rowsum = 0;
   for(int rid = 0; rid < ref_num; rid++)
    if(ctx_obj_ct[ ((llong)b*num_cof_batch + qid)*ref_num + rid ].shared_k_ct > 0)
     abundence_rowsum += (double)ctx_obj_ct[ ((llong)b*num_cof_batch + qid) * ref_num + rid ].shared_koc_ct
              / (double)ctx_obj_ct[ ((llong)b*num_cof_batch + qid)*ref_num + rid ].shared_k_ct;
   Y_size = qry_ctx_ct_list[b*num_cof_batch + qid];
      for(int rid = 0; rid < ref_num; rid++) {
        X_size = ref_ctx_ct_list[rid];
        Min_XY_size = X_size < Y_size ? X_size : Y_size ;
        XnY_size = ctx_obj_ct[ ((llong)b*num_cof_batch + qid)*ref_num + rid ].shared_k_ct;
    double abundence_pct = (double)ctx_obj_ct[ ((llong)b*num_cof_batch + qid)*ref_num + rid ].shared_koc_ct
               / XnY_size;
        XuY_size = X_size + Y_size - XnY_size ;
        X_XnY_size = X_size - XnY_size;
        Y_XnY_size = Y_size- XnY_size;
        jac_ind = (double)XnY_size / XuY_size;
        contain_ind = (double)XnY_size / Min_XY_size ;
        Dm = jac_ind == 1? 0: -log(2*jac_ind/(1+jac_ind)) / kmerlen ;
        Da = contain_ind== 1? 0: -log(contain_ind) / kmerlen ;
        P_K_in_X_XnY = 1 - pow( (1- 1/pow(alp_size,(kmerlen - dim_reduct_len) )), X_XnY_size );
        P_K_in_Y_XnY = 1 - pow( (1- 1/pow(alp_size,(kmerlen - dim_reduct_len) )), Y_XnY_size );
        double rs = P_K_in_X_XnY * P_K_in_Y_XnY * ( X_XnY_size + Y_XnY_size )
              /(P_K_in_X_XnY + P_K_in_Y_XnY - 2*P_K_in_X_XnY * P_K_in_Y_XnY);
        j_prim = ((double)XnY_size - rs) / XuY_size ;
        c_prim = ((double)XnY_size - rs) / Min_XY_size ;
        Dm_prim = j_prim == 1? 0:-log(2*j_prim/(1+j_prim)) / kmerlen ;
        Da_prim = c_prim==1? 0:-log(c_prim) / kmerlen ;
        sd_j_prim = pow(j_prim*(1 - j_prim) / XuY_size, 0.5) ;
        sd_c_prim = pow(c_prim*(1 - c_prim) / Min_XY_size,0.5) ;
        CI95_j_prim1 = j_prim - 1.96*sd_j_prim;
        CI95_j_prim2 = j_prim + 1.96*sd_j_prim;
        CI95_c_prim1 = c_prim - 1.96*sd_c_prim;
        CI95_c_prim2 = c_prim + 1.96*sd_c_prim;
        CI95_Dm_prim1 = CI95_j_prim2 == 1? 0:-log(2*CI95_j_prim2/(1+CI95_j_prim2)) / kmerlen ;
        CI95_Dm_prim2 = CI95_j_prim1 == 1? 0:-log(2*CI95_j_prim1/(1+CI95_j_prim1)) / kmerlen ;
        CI95_Da_prim1 = CI95_c_prim2 == 1? 0:-log(CI95_c_prim2) / kmerlen ;
        CI95_Da_prim2 = CI95_c_prim1 == 1? 0:-log(CI95_c_prim1) / kmerlen ;
        double pv_j_prim = 0.5 * erfc( j_prim / sd_j_prim * pow(0.5,0.5) );
        double pv_c_prim = 0.5 * erfc( c_prim/sd_c_prim * pow(0.5,0.5) );
        double qv_j_prim = pv_j_prim * ref_seq_num*qry_seq_num ;
        double qv_c_prim = pv_c_prim * ref_seq_num*qry_seq_num ;
        fprintf(distfp,"%s\t%s\t%lf\t%u-%u|%u|%u\t%lf\t%lf\t%lf\t%lf\t%lf[%lf,%lf]\t%lf[%lf,%lf]\t%lf[%lf,%lf]\t%lf[%lf,%lf]\t%E\t%E\t%E\t%E\n",
          qryfname[b*num_cof_batch + qid],refname[rid],abundence_pct,XnY_size,(unsigned int)rs,X_size,Y_size,jac_ind,Dm,
          contain_ind,Da,j_prim,CI95_j_prim1,CI95_j_prim2,Dm_prim,CI95_Dm_prim1,CI95_Dm_prim2,c_prim,CI95_c_prim1,
          CI95_c_prim2,Da_prim,CI95_Da_prim1,CI95_Da_prim2,pv_j_prim,pv_c_prim,qv_j_prim,qv_c_prim);
      }
    }
    munmap(ctx_obj_ct + (size_t)b*num_cof_batch*ref_num, maplength);
  }
  fclose(distfp);
}
void dist_print_nobin (const char *distout_dir,unsigned int ref_num, unsigned int qry_num, unsigned int*ref_ctx_ct_list,
      unsigned int*qry_ctx_ct_list, int num_cof_batch, char (*refname)[PATHLEN], char (*qryfname)[PATHLEN],dist_opt_val_t *opt_val)
{
 if( opt_val->shared_kmerpath[0] != '\0') strcpy(full_distfcode, opt_val->shared_kmerpath );
  else sprintf(full_distfcode,"%s/sharedk_ct.dat",distout_dir);
  int fd;
  struct stat s;
  fd = open (full_distfcode, O_RDONLY);
  check (fd < 0, "open %s failed: %s", full_distfcode, strerror (errno));
  fstat (fd, & s);
  ctx_obj_ct_t * ctx_obj_ct = mmap(0, s.st_size , PROT_READ, MAP_PRIVATE, fd, 0);
  check ( ctx_obj_ct == MAP_FAILED, "mmap %s failed: %s", full_distfcode, strerror (errno));
  close(fd);
 int p_fit_mem= opt_val->p;
 prt_line_t * prt_buf = malloc(ref_num * sizeof(prt_line_t));
 char distf[PATHLEN];
  sprintf(distf, "%s/distance.out", distout_dir);
 FILE *distfp = fopen(distf,"w") ;
  if( distfp == NULL ) err(errno,"dist_print_nobin():%s",distf);
 int num_mapping_distf = qry_num / num_cof_batch; int bnum_infile;
 print_ctrl_t outfield;
 outfield.metric = opt_val->metric ;
 outfield.pfield = opt_val->outfields;
 outfield.correction = opt_val->correction;
 outfield.dthreshold = opt_val->mut_dist_max;
 outfield.cmprsn_num = ref_num*qry_num;
#define MAX_CNAME_SIZE 30
 char header[2][3][MAX_CNAME_SIZE] = {
  {"Jaccard\tMashD","P-value(J)\tFDR(J)","Jaccard_CI\tMashD_CI" },
  {"ContainmentM\tAafD","P-value(C)\tFDR(C)","ContainmentM_CI\tAafD_CI" }
 };
 fprintf(distfp,"Qry\tRef\tShared_k|Ref_s|Qry_s");
 for( int i = 0 ; i<= (int)outfield.pfield ; i++)
  fprintf(distfp,"\t%s", header[outfield.metric][i]);
 fprintf(distfp,"\n");
#define NREF 1024
 int N_max = opt_val->num_neigb;
 if( (N_max > NREF) || (N_max > ref_num) ) { err(errno,"neighborN_max %d should smaller than NREF %d and ref_num %d",N_max,NREF,ref_num); } ;
 typedef struct { double metric; int rid; } Nref_stuct;
 Nref_stuct bestNref[NREF];
 for(int b=0;b<=num_mapping_distf;b++){
    if(b==num_mapping_distf){
      bnum_infile = qry_num % num_cof_batch ;
      if( bnum_infile == 0 ) continue;
    }else bnum_infile = num_cof_batch;
    size_t maplength = (size_t)bnum_infile * ref_num * sizeof(ctx_obj_ct_t);
  for(int qid = 0; qid < bnum_infile; qid++) {
   outfield.Y_size = qry_ctx_ct_list[b*num_cof_batch + qid];
   outfield.qname = qryfname[b*num_cof_batch + qid];
   outfield.qry_len = strlen(outfield.qname);
   llong offset = ((llong)b*num_cof_batch + qid)*ref_num;
   if( N_max ) {
    for (int i=0; i< N_max ;i++ ) bestNref[i] = (Nref_stuct){0, -1};
    for(int rid = 0; rid < ref_num; rid++) {
     unsigned int X_size = ref_ctx_ct_list[rid];
         unsigned int XnY_size = ctx_obj_ct[offset + rid];
         double metric = outfield.metric == Ctm ?
        (double) XnY_size / (X_size < outfield.Y_size ? X_size : outfield.Y_size) :
        (double) XnY_size / (X_size + outfield.Y_size - XnY_size) ;
     for(int i = N_max - 1 ; i>=0; i-- ){
      if(metric > bestNref[i].metric){
        bestNref[i+1] = bestNref[i];
        bestNref[i] = (Nref_stuct){ metric, rid };
      }
      else break;
     }
    }
#pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
    for( int i = 0 ; i< N_max;i++ ){
     if (bestNref[i].rid != -1)
       output_ctrl( ref_ctx_ct_list[ bestNref[i].rid ], ctx_obj_ct[offset + bestNref[i].rid], &outfield, refname[bestNref[i].rid], &prt_buf[i]);
     else prt_buf[i].len = 0 ;
    }
    for( int i = 0 ; i< N_max;i++ )
          if(prt_buf[i].len >1) fwrite(prt_buf[i].line, prt_buf[i].len, 1, distfp);;
   }
   else {
#pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
    for(int rid = 0; rid < ref_num; rid++)
     output_ctrl( ref_ctx_ct_list[ rid ], ctx_obj_ct[ offset + rid ], &outfield, refname[rid], &prt_buf[rid]);
    for(int rid = 0; rid < ref_num; rid++)
     if(prt_buf[rid].len >1) fwrite(prt_buf[rid].line, prt_buf[rid].len, 1, distfp);;
   }
  }
  munmap(ctx_obj_ct + (size_t)b*num_cof_batch*ref_num, maplength);
 }
 fclose(distfp);
  free(prt_buf);
 if(!opt_val->keep_shared_kmer) remove(full_distfcode);
}
#define GET_MATRIC(X,Y) ( (X) == Jcd ? 1 / (2*(Y)) + 0.5 : 1 / (Y) )
static inline void output_ctrl (unsigned int X_size, unsigned int XnY_size, print_ctrl_t* outfield, char *rname, prt_line_t* linebuf ){
 double rs = 0;
 if(outfield->correction){
  unsigned int X_XnY_size = X_size - XnY_size;
  unsigned int Y_XnY_size = outfield->Y_size - XnY_size;
  double P_K_in_X_XnY = 1 - pow( (1- 1/pow(alp_size,(kmerlen - dim_reduct_len) )), X_XnY_size );
  double P_K_in_Y_XnY = 1 - pow( (1- 1/pow(alp_size,(kmerlen - dim_reduct_len) )), Y_XnY_size );
  rs = P_K_in_X_XnY * P_K_in_Y_XnY * ( X_XnY_size + Y_XnY_size )
       /(P_K_in_X_XnY + P_K_in_Y_XnY - 2*P_K_in_X_XnY * P_K_in_Y_XnY);
 }
 unsigned int tmp = outfield->metric == Jcd ? X_size + outfield->Y_size - XnY_size
                   :(X_size < outfield->Y_size ? X_size : outfield->Y_size);
 double metric = ((double)XnY_size - rs) / tmp;
 double dist = log( GET_MATRIC(outfield->metric,metric) ) / kmerlen;
 if (dist > 1) dist = 1 ;
  if( dist > outfield->dthreshold ) { linebuf->len = 0 ; return; };
#define LST_PRT_LEN 27
 snprintf(linebuf->line,LINE_LEN,"%s\t%s\t%u-%u|%u|%u\t%.6lf\t%.6lf",outfield->qname,rname,XnY_size,(unsigned int)rs,X_size, outfield->Y_size, metric, dist);
 linebuf->len = LST_PRT_LEN + outfield->qry_len + strlen(linebuf->line + LST_PRT_LEN + outfield->qry_len);
 if(outfield->pfield > Dst) {
  double sd = pow(metric*(1 - metric) / tmp, 0.5) ;
  double pv = 0.5 * erfc( metric / sd * pow(0.5,0.5) ) ;
  snprintf(linebuf->line + linebuf->len, LINE_LEN - linebuf->len,"\t%E\t%E", pv, pv * outfield->cmprsn_num );
  linebuf->len += 4 + strlen(linebuf->line + linebuf->len + 4);
  if(outfield->pfield >Qv){
   double CI95_mtrc_1 = metric - 1.96*sd;
   double CI95_mtrc_2 = metric + 1.96*sd;
   double CI95_dist_1 = log( GET_MATRIC(outfield->metric,CI95_mtrc_2) ) / kmerlen;
   double CI95_dist_2 = log( GET_MATRIC(outfield->metric,CI95_mtrc_1) ) / kmerlen;
   snprintf(linebuf->line + linebuf->len, LINE_LEN - linebuf->len,"\t[%.6lf,%.6lf]\t[%.6lf,%.6lf]", CI95_mtrc_1, CI95_mtrc_2, CI95_dist_1, CI95_dist_2);
   linebuf->len += 20 + strlen(linebuf->line + linebuf->len + 20);
  }
 }
 snprintf(linebuf->line + linebuf->len, LINE_LEN - linebuf->len,"\n");
 linebuf->len += 1;
}
infile_tab_t* dist_organize_infiles (dist_opt_val_t *opt_val)
{
 int fmt_ck;
 if(opt_val->pipecmd[0]=='\0')
  fmt_ck = 1;
 else
  fmt_ck = 0;
 if( strcmp( opt_val->fpath, "" ) !=0 )
 {
  return organize_infile_list(opt_val->fpath,fmt_ck);
 }
 else if( opt_val->num_remaining_args > 0 )
 {
  return organize_infile_frm_arg(opt_val->num_remaining_args, opt_val->remaining_args,fmt_ck);
 }
 else
 {
  perror("please specify the input/query files");
    return NULL;
 }
};
extern const char *acpt_infile_fmt[ACPT_FMT_SZ] ;
infile_tab_t* dist_organize_refpath( dist_opt_val_t *opt_val){
 struct stat path_stat;
 if( stat(opt_val->refpath, &path_stat) < 0)
  err(errno,"dist_organize_refpath():%s",opt_val->refpath);
 if( S_ISDIR(path_stat.st_mode) || isOK_fmt_infile(opt_val->refpath, acpt_infile_fmt,ACPT_FMT_SZ) ){
  char * tmp_arg[] = {opt_val->refpath};
  return organize_infile_frm_arg(1, tmp_arg,1);
 }
 else if(S_ISREG(path_stat.st_mode))
  return organize_infile_list(opt_val->refpath,1);
 else
  return NULL;
}
const char * combine_queries(dist_opt_val_t *opt_val)
{
 int p_fit_mem = opt_val->p;
 const char* co_dir = opt_val->outdir;
 mkdir(co_dir,0700);
 const char kmerct_list_fname[] = "tmp_kmerct_list";
 const char fname_fname[] = "tmp_infiles_fname";
 if(opt_val->abundance){err(errno,"combine_queries(): abundance model not supported yet");}
 FILE *co_stat_fp;
 const char *qryco_dstat_fpath = NULL;
 qryco_dstat_fpath = test_get_fullpath(opt_val->remaining_args[0], co_dstat);
 if( ( co_stat_fp = fopen(qryco_dstat_fpath,"rb")) == NULL ){
  err(errno,"combine_queries():%s", qryco_dstat_fpath);
 }
 co_dstat_t co_dstat_one, co_dstat_it;
 fread(&co_dstat_one, sizeof(co_dstat_t), 1, co_stat_fp);
 if (co_dstat_one.koc) { err(errno,"combine_queries(): abundance model not supported yet"); }
 char one_kmerct_list_name[PATHLEN];
 sprintf(one_kmerct_list_name,"%s/%s",co_dir,kmerct_list_fname);
 FILE *kmerct_list_fp;
 if( (kmerct_list_fp = fopen(one_kmerct_list_name,"wb") ) == NULL ) {err(errno,"%s",one_kmerct_list_name);}
 ctx_obj_ct_t * tmp_ct_list = malloc(sizeof(ctx_obj_ct_t) * co_dstat_one.infile_num);
 fread(tmp_ct_list,sizeof(ctx_obj_ct_t),co_dstat_one.infile_num,co_stat_fp);
 fwrite(tmp_ct_list,sizeof(ctx_obj_ct_t),co_dstat_one.infile_num,kmerct_list_fp);
 char one_infilename_name[PATHLEN];
 sprintf(one_infilename_name,"%s/%s",co_dir,fname_fname);
 FILE *infilename_name_fp;
 if( ( infilename_name_fp = fopen(one_infilename_name,"wb") ) == NULL ){err(errno,"%s", one_infilename_name);}
 char (*tmpname)[PATHLEN] = malloc(PATHLEN * co_dstat_one.infile_num);
 fread(tmpname,PATHLEN,co_dstat_one.infile_num,co_stat_fp);
 fwrite(tmpname,PATHLEN,co_dstat_one.infile_num,infilename_name_fp);
 fclose(co_stat_fp);
 FILE** com_cofp = malloc( sizeof(FILE*) * co_dstat_one.comp_num);
  FILE** indexfp = malloc( sizeof(FILE*) * co_dstat_one.comp_num);
 size_t *index_offset = malloc(sizeof(size_t) * co_dstat_one.comp_num);
#pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
 for(int c = 0; c < co_dstat_one.comp_num; c++) {
  struct stat file_stat;
  char combined_cof[PATHLEN];
  sprintf(combined_cof,"%s/combco.%d",co_dir,c);
  if( (com_cofp[c] = fopen(combined_cof,"wb")) == NULL) err(errno,"%s",combined_cof);
  sprintf(combined_cof,"%s/combco.%d",opt_val->remaining_args[0],c);
  stat(combined_cof, &file_stat);
  unsigned int *tmp_combco = malloc(file_stat.st_size);
  FILE *com_cofp_it;
  if( (com_cofp_it = fopen(combined_cof,"rb")) == NULL) err(errno,"%s",combined_cof);
  fread(tmp_combco, file_stat.st_size, 1, com_cofp_it);
    fwrite(tmp_combco, file_stat.st_size, 1, com_cofp[c]);
    fclose(com_cofp_it);
    free(tmp_combco);
  char indexfname[PATHLEN];
  sprintf(indexfname,"%s/combco.index.%d",co_dir,c);
    if( (indexfp[c] = fopen(indexfname,"wb")) == NULL) err(errno,"%s",indexfname);
  sprintf(indexfname, "%s/combco.index.%d", opt_val->remaining_args[0], c);
  stat(indexfname, &file_stat);
  size_t *tmp_index_combco = malloc(file_stat.st_size);
  FILE *com_indexfp_it;
  if((com_indexfp_it = fopen(indexfname,"rb") ) == NULL) err(errno,"%s",indexfname);
  fread(tmp_index_combco,file_stat.st_size, 1,com_indexfp_it);
  fwrite(tmp_index_combco,file_stat.st_size, 1,indexfp[c]);
  index_offset[c] = tmp_index_combco[file_stat.st_size/sizeof(size_t) - 1] ;
  fclose(com_indexfp_it);
  free(tmp_index_combco);
 }
 for (int i = 1; i< opt_val->num_remaining_args; i++){
  qryco_dstat_fpath = test_get_fullpath(opt_val->remaining_args[i],co_dstat);
  if(qryco_dstat_fpath == NULL){
     printf("%dth query %s is not a valid query: no %s file\n", i, opt_val->remaining_args[i], co_dstat);
      continue;
    }
  else if( ( co_stat_fp = fopen(qryco_dstat_fpath,"rb")) == NULL ){
     printf("combine_queries(): %dth query can not open %s\n", i, qryco_dstat_fpath);
   continue;
   }
  fread(&co_dstat_it, sizeof(co_dstat_t), 1, co_stat_fp);
  if(co_dstat_one.shuf_id != co_dstat_it.shuf_id){
      printf("combine_queries(): %dth shuf_id: %u not match 0th shuf_id: %u\n",i, co_dstat_it.shuf_id, co_dstat_one.shuf_id);
      fclose(co_stat_fp);
      continue;
    }
  else if (co_dstat_it.koc){
   printf("combine_queries(): %dth query abundance model not supported yet \n", i);
      fclose(co_stat_fp);
      continue;
  }
    co_dstat_one.all_ctx_ct += co_dstat_it.all_ctx_ct;
    co_dstat_one.infile_num += co_dstat_it.infile_num;
  tmp_ct_list = realloc(tmp_ct_list, sizeof(ctx_obj_ct_t) * co_dstat_it.infile_num);
  fread(tmp_ct_list,sizeof(ctx_obj_ct_t),co_dstat_it.infile_num,co_stat_fp);
  fwrite(tmp_ct_list,sizeof(ctx_obj_ct_t),co_dstat_it.infile_num,kmerct_list_fp);
  tmpname = realloc(tmpname, PATHLEN * co_dstat_it.infile_num ) ;
  fread(tmpname,PATHLEN,co_dstat_it.infile_num,co_stat_fp);
  fwrite(tmpname,PATHLEN,co_dstat_it.infile_num,infilename_name_fp);
    fclose(co_stat_fp);
#pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
  for(int c = 0; c < co_dstat_one.comp_num; c++) {
   char combined_cof_it[PATHLEN];
   FILE *com_cofp_it;
   struct stat file_stat;
   sprintf(combined_cof_it,"%s/combco.%d",opt_val->remaining_args[i],c);
   stat(combined_cof_it,&file_stat);
   unsigned int *tmpco = malloc(file_stat.st_size);
      if( (com_cofp_it = fopen(combined_cof_it,"rb")) == NULL) err(errno,"%s",combined_cof_it);
   fread(tmpco,file_stat.st_size, 1, com_cofp_it);
   fwrite(tmpco,file_stat.st_size, 1, com_cofp[c]);
   fclose(com_cofp_it);
   free(tmpco);
   char indexfname_it[PATHLEN];
   sprintf(indexfname_it,"%s/combco.index.%d",opt_val->remaining_args[i],c);
   stat(indexfname_it, &file_stat);
   size_t * tmpindex = malloc(file_stat.st_size);
   FILE *indexfp_it;
   if( (indexfp_it = fopen(indexfname_it,"rb")) == NULL) err(errno,"%s",indexfname_it);
   fread(tmpindex,file_stat.st_size, 1,indexfp_it);
   fclose(indexfp_it);
   int tmp_infile_num = file_stat.st_size/sizeof(size_t);
   for(int i=1; i< tmp_infile_num; i++ ) { tmpindex[i] += index_offset[c];}
   fwrite( tmpindex + 1, sizeof(size_t), tmp_infile_num - 1, indexfp[c]);
   index_offset[c] = tmpindex[tmp_infile_num-1];
  }
 }
 fclose(kmerct_list_fp);
 fclose(infilename_name_fp);
#pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
  for(int c = 0; c < co_dstat_one.comp_num; c++) {
  fclose(com_cofp[c]);
  fclose(indexfp[c]);
 }
 FILE *one_co_stat_fp;
 char one_stat_name[PATHLEN];
 sprintf(one_stat_name,"%s/%s",co_dir,co_dstat);
 if( (one_co_stat_fp = fopen(one_stat_name,"wb")) == NULL) { err(errno,"%s",one_stat_name) ;};
 fwrite(&co_dstat_one,sizeof(co_dstat_t),1,one_co_stat_fp);
 struct stat file_stat;
 stat(one_kmerct_list_name, &file_stat);
 tmp_ct_list = realloc(tmp_ct_list, file_stat.st_size);
 if( (kmerct_list_fp = fopen(one_kmerct_list_name,"rb") ) == NULL ) {err(errno,"%s",one_kmerct_list_name);}
 fread(tmp_ct_list,file_stat.st_size,1,kmerct_list_fp);
 fwrite(tmp_ct_list,file_stat.st_size,1,one_co_stat_fp);
 stat(one_infilename_name,&file_stat);
 tmpname = realloc(tmpname,file_stat.st_size);
 if( (infilename_name_fp = fopen(one_infilename_name,"rb") ) == NULL ){err(errno,"%s", one_infilename_name);}
 fread(tmpname,file_stat.st_size,1,infilename_name_fp);
 fwrite(tmpname,file_stat.st_size,1,one_co_stat_fp);
 fclose(one_co_stat_fp);
 remove(one_kmerct_list_name);
 remove(one_infilename_name);
 free(tmp_ct_list);
 free(tmpname);
 free(com_cofp);
  free(indexfp);
 return co_dir;
}
