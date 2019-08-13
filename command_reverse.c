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
   
#include "command_reverse.h"
#include "global_basic.h"
#include "command_shuffle.h"
#include "command_dist.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <argp.h>
#include <argz.h>
#include <err.h>
#include <errno.h>
#include <math.h>
struct arg_reverse
{
  struct arg_global* global;
  char* name;
};
static struct argp_option opt_reverse[] =
{
 {"shufFile",'L',"<path>",0,"provide .shuf file.\v"},
 {"outdir",'o',"<path>",0,"path for recovered k-mer files.\v"},
 {"threads",'p',"INT",0,"threads num.\v"},
  { 0 }
};
static char doc_reverse[] =
  "\n"
  "The reverse doc prefix."
  "\v"
  "The reverse doc suffix."
  ;
reverse_opt_val_t reverse_opt_val =
{
"",
".",
};
static error_t parse_reverse(int key, char* arg, struct argp_state* state) {
  struct arg_reverse* reverse = state->input;
  assert( reverse );
  assert( reverse->global );
  switch(key)
  {
    case 'L':
  {
   struct stat path_stat;
   if( stat(arg,&path_stat) >=0 && S_ISREG(path_stat.st_mode)){
        if(strlen(arg) < PATHLEN )
          strcpy(reverse_opt_val.shufile,arg);
        else
           err(errno,"-L argument path should not longer than %d",PATHLEN);
      }
   break;
  }
  case 'o':
    {
      if(strlen(arg) > PATHLEN) {
        err(errno,"the outdir path should not longer than %d", PATHLEN);
        exit(EXIT_FAILURE);
      };
      strcpy(reverse_opt_val.outdir,arg);
      break;
    }
  case 'p':
    {
#ifdef _OPENMP
        reverse_opt_val.p = atoi(arg) ;
#else
      warnx("This version of kssd was built without OpenMP and "
          "thus does not support multi threading. Ignoring -p %d",atoi(arg));
      break;
#endif
    }
    break;
  case ARGP_KEY_ARGS:
    {
    reverse_opt_val.num_remaining_args = state->argc - state->next;
     reverse_opt_val.remaining_args = state->argv + state->next;
     break;
  }
    case ARGP_KEY_NO_ARGS:
    {
      if(state->argc < 2)
        {
        printf("\v");
        argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
        printf("\v");
        argp_state_help(state,stdout,ARGP_HELP_LONG);
        printf("\v");
        exit(0);
        };
        return EINVAL;
    }
   break;
   default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
};
static struct argp argp_reverse =
{
  opt_reverse,
  parse_reverse,
 "<co dir>",
  doc_reverse
};
int cmd_reverse(struct argp_state* state)
{
  int argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  char* argv0 = argv[0];
 argv[0] = malloc(strlen(state->name) + strlen(" reverse") + 1);
  if(!argv[0])
    argp_failure(state, 1, ENOMEM, 0);
 sprintf(argv[0], "%s reverse", state->name);
  argp_parse(&argp_reverse, argc, argv, ARGP_IN_ORDER, &argc, &shuffle);
  free(argv[0]);
  argv[0] = argv0;
  state->next += argc - 1;
  return co_reverse2kmer(&reverse_opt_val);
}
typedef unsigned int ctx_obj_ct_t;
int co_reverse2kmer(reverse_opt_val_t *opt_val)
{
 dim_shuffle_t* shuf_arr = read_dim_shuffle_file(opt_val->shufile);
 int shuf_arr_len = 1LLU << (4 * shuf_arr->dim_shuffle_stat.subk) ;
 unsigned int rev_shuf_arr[MIN_SUBCTX_DIM_SMP_SZ];
 int count = 0;
 for(unsigned int i=0; i< shuf_arr_len; i++) {
  if( shuf_arr->shuffled_dim[i] < MIN_SUBCTX_DIM_SMP_SZ ){
   rev_shuf_arr[shuf_arr->shuffled_dim[i]] = i ;
   count++;
  }
 }
 if(count != MIN_SUBCTX_DIM_SMP_SZ) err(errno,"count %d not match MIN_SUBCTX_DIM_SMP_SZ %d",count,MIN_SUBCTX_DIM_SMP_SZ);
 int comp_code_bits = shuf_arr->dim_shuffle_stat.k - shuf_arr->dim_shuffle_stat.drlevel > COMPONENT_SZ
           ? 4*(shuf_arr->dim_shuffle_stat.k - shuf_arr->dim_shuffle_stat.drlevel - COMPONENT_SZ ) : 0 ;
 int inner_ctx_bits = shuf_arr->dim_shuffle_stat.subk * 4;
 int half_outer_ctx_bits = (shuf_arr->dim_shuffle_stat.k - shuf_arr->dim_shuffle_stat.subk) *2 ;
 int pf_bits = ( shuf_arr->dim_shuffle_stat.subk - shuf_arr->dim_shuffle_stat.drlevel ) * 4;
 int TL = shuf_arr->dim_shuffle_stat.k * 2;
  if (!(opt_val->num_remaining_args >0 ))
  err(errno,"need speficy the query path");
 const char *qryco_dstat_fpath = NULL;
  qryco_dstat_fpath = test_get_fullpath(opt_val->remaining_args[0],co_dstat);
 if( qryco_dstat_fpath == NULL )
  err(errno,"%s is not a valid query folder",opt_val->remaining_args[0]);
 FILE *qry_co_stat_fp;
 if (( qry_co_stat_fp = fopen(qryco_dstat_fpath,"rb")) == NULL) err(errno,"qry co stat file:%s",qryco_dstat_fpath);
 char *qryco_dname = opt_val->remaining_args[0];
 co_dstat_t co_qry_dstat;
 fread( &co_qry_dstat, sizeof(co_dstat_t), 1, qry_co_stat_fp);
 ctx_obj_ct_t * qry_ctx_ct_list = malloc(co_qry_dstat.infile_num * sizeof(ctx_obj_ct_t));
 fread(qry_ctx_ct_list,sizeof(ctx_obj_ct_t),co_qry_dstat.infile_num,qry_co_stat_fp);
 char (*cofname)[PATHLEN] = malloc(co_qry_dstat.infile_num * PATHLEN);
 fread(cofname,PATHLEN,co_qry_dstat.infile_num,qry_co_stat_fp);
 fclose(qry_co_stat_fp);
 FILE *cbd_fcode_comp_fp,*cbd_fcode_comp_index_fp;
  struct stat cbd_fcode_stat;
 size_t *fco_pos = malloc(sizeof(size_t) * (co_qry_dstat.infile_num + 1) );
 char co_cbd_fcode[PATHLEN];char co_cbd_index_fcode[PATHLEN];
 llong **kmer = malloc( co_qry_dstat.infile_num * sizeof(llong*) );
 int *filled_len = calloc( co_qry_dstat.infile_num, sizeof(int));
 for(int k = 0; k < co_qry_dstat.infile_num; k++){
  kmer[k] = malloc( sizeof(llong) * qry_ctx_ct_list[k] );
 }
 int p_fit_mem = opt_val->p ;
 for ( int j = 0; j < co_qry_dstat.comp_num; j++ ) {
  sprintf(co_cbd_fcode,"%s/combco.%d",qryco_dname,j);
  if( (cbd_fcode_comp_fp = fopen(co_cbd_fcode,"rb"))==NULL) err(errno,"co_reverse2kmer()::%s",co_cbd_fcode);
    stat(co_cbd_fcode, &cbd_fcode_stat);
  unsigned int *cbd_fcode_mem = malloc(cbd_fcode_stat.st_size);
  fread(cbd_fcode_mem,sizeof(unsigned int),cbd_fcode_stat.st_size/sizeof(unsigned int),cbd_fcode_comp_fp);
  fclose(cbd_fcode_comp_fp);
  sprintf(co_cbd_index_fcode,"%s/combco.index.%d",qryco_dname,j);
  if( (cbd_fcode_comp_index_fp = fopen(co_cbd_index_fcode,"rb"))==NULL)
        err(errno,"co_reverse2kmer()::%s",co_cbd_index_fcode);
  fread(fco_pos,sizeof(size_t),co_qry_dstat.infile_num + 1 ,cbd_fcode_comp_index_fp);
  fclose(cbd_fcode_comp_index_fp);
  #pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
  for(int k = 0; k < co_qry_dstat.infile_num; k++){
   if(qry_ctx_ct_list[k]==0) continue;
   char *kstring = malloc(TL + 1);
   kstring[TL] = '\0';
   for(int n = 0; n < fco_pos[k+1] - fco_pos[k]; n++){
    int ind = cbd_fcode_mem[ fco_pos[k] + n ];
    kmer[k][filled_len[k] + n] = core_reverse2unituple(ind,j,comp_code_bits,pf_bits,inner_ctx_bits,half_outer_ctx_bits,rev_shuf_arr);
   }
   filled_len[k] += (fco_pos[k+1] - fco_pos[k]) ;
  }
 }
 #pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
 for(int k = 0; k < co_qry_dstat.infile_num; k++) {
  if(qry_ctx_ct_list[k]==0) continue;
    char *kstring = malloc(TL + 1);
    kstring[TL] = '\0';
  char *filename;
  (filename = strrchr(cofname[k],'/') ) ? ++filename : (filename = cofname[k]);
  char fullname[PATHLEN];
  sprintf(fullname,"%s/%s",opt_val->outdir,filename);
  FILE *kmerf;
  if ( ( kmerf = fopen(fullname,"w")) == NULL ) err(errno,"%s",filename);
  for (int n =0 ; n< qry_ctx_ct_list[k]; n++){
   llong unituple = kmer[k][n];
   for(int i=0; i<TL;i++){
    kstring[TL-i-1] = Mapbase[unituple % 4] ;
    unituple >>= 2 ;
   }
   fprintf(kmerf,"%s\n",kstring);
  }
  fclose(kmerf);
 }
 return co_qry_dstat.infile_num;
}
llong core_reverse2unituple(unsigned int kid, int compid, int compbit, int pf_bits, int inner_ctx_bits, int half_outer_ctx_bits, unsigned int *rev_shuf_arr)
{
 llong drtuple = ( ((llong)kid) << compbit ) + compid ;
 unsigned int ind = rev_shuf_arr[drtuple % MIN_SUBCTX_DIM_SMP_SZ];
 llong tuple = ((drtuple >> pf_bits) << inner_ctx_bits) + (llong)ind;
 llong half_outer_ctx_mask = ( (1LLU << half_outer_ctx_bits) - 1 ) << inner_ctx_bits ;
  llong unituple = (tuple & (half_outer_ctx_mask << half_outer_ctx_bits))
   + ((tuple & half_outer_ctx_mask) >> inner_ctx_bits )
  + ( (tuple & ( (1LLU << inner_ctx_bits) - 1 )) << half_outer_ctx_bits );
 return unituple;
}
