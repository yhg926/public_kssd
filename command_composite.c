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
   
#include "command_composite.h"
#include "global_basic.h"
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
struct arg_composite
{
  struct arg_global* global;
  char* name;
};
static struct argp_option opt_composite[] =
{
 {"ref",'r',"<DIR>", 0, "Path of species specific pan uniq kmer database (reference) \v",1},
 {"query",'q',"<DIR>", 0, "Path of query sketches with abundances \v",2},
 {"outfile",'o',"<DIR>",0,"Output path \v",3},
 {"threads",'p',"<INT>", 0, "Threads number to use \v",4},
  { 0 }
};
static char doc_composite[] =
  "\n"
  "The composite doc prefix."
  "\v"
  "The composite doc suffix."
  ;
composite_opt_t composite_opt ={
 .p = 1,
 .refdir[0] = '\0',
 .qrydir[0] = '\0',
 .outdir = "./"
};
static error_t parse_composite(int key, char* arg, struct argp_state* state) {
  struct arg_composite* composite = state->input;
  assert( composite );
  assert( composite->global );
  switch(key)
  {
  case 'p':
  {
   composite_opt.p = atoi(arg);
   break;
  }
  case 'r':
  {
   strcpy(composite_opt.refdir, arg);
   break;
  }
  case 'q':
  {
   strcpy(composite_opt.qrydir, arg);
   break;
  }
  case 'o':
  {
   strcpy(composite_opt.outdir, arg);
   break;
  }
    case ARGP_KEY_NO_ARGS:
    {
   if(state->argc<2)
   {
       printf("\v");
    argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
    printf("\v");
       argp_state_help(state,stdout,ARGP_HELP_LONG);
       printf("\v");
       return EINVAL;
   }
    }
  break;
    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}
static struct argp argp_composite =
{
  opt_composite,
  parse_composite,
 0,
  doc_composite
};
int cmd_composite(struct argp_state* state)
{
  struct arg_composite composite = { 0, };
  int argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  char* argv0 = argv[0];
  composite.global = state->input;
 argv[0] = malloc(strlen(state->name) + strlen(" composite") + 1);
  if(!argv[0])
    argp_failure(state, 1, ENOMEM, 0);
 sprintf(argv[0], "%s composite", state->name);
  argp_parse(&argp_composite, argc, argv, ARGP_IN_ORDER, &argc, &composite);
  free(argv[0]);
  argv[0] = argv0;
  state->next += argc - 1;
 if(argc >1)
  return get_species_abundance (&composite_opt);
 else
  return -1;
}
int **ref_abund;
int get_species_abundance (composite_opt_t * composite_opt) {
 if ( composite_opt->refdir[0] == '\0' || composite_opt->qrydir[0] == '\0' || strcmp (composite_opt->refdir, composite_opt->qrydir) == 0 )
  err(errno, "get_species_abundance(): refdir or qrydir is not initialized\n" );
 const char *ref_dstat_fpath = test_get_fullpath(composite_opt->refdir,co_dstat);
 if (ref_dstat_fpath == NULL) err(errno,"cannot find %s under %s ",co_dstat, ref_dstat_fpath);
  const char *qry_dstat_fpath = test_get_fullpath(composite_opt->qrydir,co_dstat);
 if (qry_dstat_fpath == NULL) err(errno,"cannot find %s under %s ",co_dstat, qry_dstat_fpath);
 FILE *ref_dstat_fp, *qry_dstat_fp ;
 if( (ref_dstat_fp = fopen(ref_dstat_fpath,"rb")) == NULL ) err(errno, "get_species_abundance():%s",ref_dstat_fpath);
 if( (qry_dstat_fp = fopen(qry_dstat_fpath,"rb")) == NULL ) err(errno, "get_species_abundance():%s",qry_dstat_fpath);
 co_dstat_t ref_dstat, qry_dstat ;
 fread( &ref_dstat, sizeof(co_dstat_t),1,ref_dstat_fp);
 fread( &qry_dstat, sizeof(co_dstat_t),1,qry_dstat_fp);
 if(!qry_dstat.koc) err(errno, "get_species_abundance(): query has not abundance");
 if(qry_dstat.shuf_id != ref_dstat.shuf_id)
  printf("get_species_abundance(): qry shuf_id %u not match ref shuf_id: %u\n",qry_dstat.shuf_id, ref_dstat.shuf_id);
 ctx_obj_ct_t* tmp_ct_list = malloc(sizeof(ctx_obj_ct_t) * ref_dstat.infile_num);
  fread(tmp_ct_list,sizeof(ctx_obj_ct_t),ref_dstat.infile_num,ref_dstat_fp);
 char (*refname)[PATHLEN] = malloc(PATHLEN * ref_dstat.infile_num);
  fread(refname,PATHLEN,ref_dstat.infile_num,ref_dstat_fp);
  fread(tmp_ct_list,sizeof(ctx_obj_ct_t),qry_dstat.infile_num,qry_dstat_fp);
  free(tmp_ct_list);
  char (*qryname)[PATHLEN] = malloc(PATHLEN * qry_dstat.infile_num);
  fread(qryname,PATHLEN,qry_dstat.infile_num, qry_dstat_fp);
 char tmpfname[PATHLEN];
 struct stat tmpstat;
 FILE *tmpfp;
#define REALLOC_US 8
 ref_abund = (int **)malloc(ref_dstat.infile_num* sizeof(int *));
 for (int i = 0; i < ref_dstat.infile_num; i++) ref_abund[i] = (int *)malloc(REALLOC_US *sizeof(int));
 for(int qn = 0; qn < qry_dstat.infile_num; qn++) {
  for (int i = 0; i < ref_dstat.infile_num; i++) ref_abund[i][0] = 0 ;
  for(int c = 0; c < ref_dstat.comp_num; c++){
   sprintf(tmpfname,"%s/%s.%d",composite_opt->refdir,skch_prefix,c);
   if( (tmpfp = fopen(tmpfname,"rb"))==NULL) err(errno,"get_species_abundance():%s",tmpfname);
   stat(tmpfname, &tmpstat);
   unsigned int *ref_combco = malloc(tmpstat.st_size);
   fread(ref_combco, tmpstat.st_size, 1, tmpfp);
   fclose(tmpfp);
   sprintf(tmpfname,"%s/%s.%d",composite_opt->refdir,idx_prefix,c);
   if( (tmpfp = fopen(tmpfname,"rb"))==NULL) err(errno,"get_species_abundance():%s",tmpfname);
   stat(tmpfname, &tmpstat);
   size_t *ref_index_combco = malloc(tmpstat.st_size);
   fread(ref_index_combco,tmpstat.st_size, 1, tmpfp);
   fclose(tmpfp);
   sprintf(tmpfname,"%s/%s.%d",composite_opt->qrydir,skch_prefix,c);
   if( (tmpfp = fopen(tmpfname,"rb"))==NULL) err(errno,"get_species_abundance():%s",tmpfname);
   stat(tmpfname, &tmpstat);
   unsigned int *qry_combco = malloc(tmpstat.st_size);
   fread(qry_combco, tmpstat.st_size, 1, tmpfp);
   fclose(tmpfp);
     sprintf(tmpfname,"%s/%s.%d",composite_opt->qrydir,idx_prefix,c);
     if( (tmpfp = fopen(tmpfname,"rb"))==NULL) err(errno,"get_species_abundance():%s",tmpfname);
     stat(tmpfname, &tmpstat);
     size_t *qry_index_combco = malloc(tmpstat.st_size);
     fread(qry_index_combco,tmpstat.st_size, 1, tmpfp);
     fclose(tmpfp);
   sprintf(tmpfname,"%s/%s.%d.a",composite_opt->qrydir,skch_prefix,c);
   if( (tmpfp = fopen(tmpfname,"rb"))==NULL) err(errno,"get_species_abundance():%s",tmpfname);
   stat(tmpfname, &tmpstat);
   unsigned short* qry_abund = malloc(tmpstat.st_size);
   fread(qry_abund, tmpstat.st_size, 1, tmpfp);
   fclose(tmpfp);
   int hash_sz = nextPrime( (int)((double)(qry_index_combco[qn+1] - qry_index_combco[qn]) / LD_FCTR) );
   size_t *km2abund_idx = calloc( hash_sz, sizeof(size_t) ) ;
   for (size_t idx = qry_index_combco[qn] ; idx < qry_index_combco[qn+1]; idx++){
    for(int i = 0; i< hash_sz; i++){
     unsigned int tmphv = HASH(qry_combco[idx],i, hash_sz);
     if(km2abund_idx[tmphv] == 0) {
      km2abund_idx[tmphv] = idx + 1;
      break;
     }
    }
   }
#pragma omp parallel for num_threads(composite_opt->p) schedule(guided)
   for(int rn = 0; rn < ref_dstat.infile_num; rn++) {
    for( size_t ri = ref_index_combco[rn]; ri < ref_index_combco[rn+1]; ri++){
     for (int i = 0; i< hash_sz; i++) {
      unsigned int hv = HASH(ref_combco[ri],i, hash_sz);
      if(km2abund_idx[hv] == 0) break;
      else if( qry_combco[km2abund_idx[hv] - 1] == ref_combco[ri] ){
       ref_abund[rn][0]++;
       ref_abund[rn][ref_abund[rn][0]] = qry_abund[km2abund_idx[hv] - 1];
       if (ref_abund[rn][0] % REALLOC_US == REALLOC_US - 1) {
              int newsize = ( (ref_abund[rn][0] / REALLOC_US) + 2 ) * REALLOC_US ;
              ref_abund[rn] = (int*)realloc( ref_abund[rn] , newsize * sizeof(int) );
       }
       break;
      }
     }
    }
   }
   free(km2abund_idx);
   free(ref_combco);
   free(ref_index_combco);
    free(qry_combco);
     free(qry_index_combco);
   free(qry_abund);
  }
#define MIN_KM_S 6
#define ST_PCTL (0.98)
#define ED_PCTL (0.99)
  int *sort_ref = malloc(ref_dstat.infile_num* sizeof(int));
   for(int i = 0; i< ref_dstat.infile_num; i++) sort_ref[i] = i;
  qsort(sort_ref, ref_dstat.infile_num, sizeof(sort_ref[0]), comparator_idx);
  for ( int i = 0; i< ref_dstat.infile_num; i++ ){
   int kmer_num = ref_abund[sort_ref[i]][0];
   if (kmer_num < MIN_KM_S) break;
   qsort(ref_abund[sort_ref[i]] + 1, kmer_num, sizeof(int),comparator);
   int sum = 0;
   for(int n = 1; n <= kmer_num; n++) sum+= ref_abund[sort_ref[i]][n];
   int median_idx = kmer_num /2;
   int pct09_idx = kmer_num * ST_PCTL ;
   int lastsum = 0; int lastn = 0;
      for(int n = pct09_idx ; n <= kmer_num*ED_PCTL; n++) {
    lastsum += ref_abund[sort_ref[i]][n];
    lastn++;
   }
   printf("%s\t%s\t%d\t%f\t%f\t%d\t%d\n",qryname[qn],refname[sort_ref[i]], kmer_num, (float)sum/kmer_num,(float)lastsum/lastn,ref_abund[sort_ref[i]][median_idx], ref_abund[sort_ref[i]][pct09_idx]);
  }
  free(sort_ref);
 }
 for (int i = 0; i < ref_dstat.infile_num; i++) free(ref_abund[i]);
  free(ref_abund) ;
 free(refname);
 free(qryname);
 return 1;
}
int comparator_idx (const void *a, const void *b){
 return ( ref_abund[*(int*)b][0] - ref_abund[*(int*)a][0] );
}
int comparator (const void *a, const void *b){
  return ( *(int*)a - *(int*)b );
}
