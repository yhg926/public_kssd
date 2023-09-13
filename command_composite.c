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
#include <libgen.h>
#include <dirent.h>
const char binVec_suffix[] = "abv";
const char abunMtx_suffix[] = "abm";
const char abunMtx_idx_suffix[] = "abmi";
const char abunMtx_name_suffix[] = "name";
const char binVec_dirname[] = "abundance_Vec";
const char y_l2n_suffix[] = "yl2n";
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
 {"binVec",'b',0,0,"Output species abundances in Binary Vector format (.abv) \v",5},
 {"idxbv",'i',0,0,"build index of abundance Binary Vector \v",6},
 {"search",'s',"<0-2>",0,"search for similar abundance Binary Vectors using L1norm(1)/L2norm(2)/cosine(0)\v",7},
 {"readabv",'d',0,0,"Read .abv file \v",8},
  { 0 }
};
static char doc_composite[] =
  "\n"
  "The composite doc prefix."
  "\v"
  "The composite doc suffix."
  ;
composite_opt_t composite_opt ={
 .b = 0,
 .i = 0,
 .s = -1,
 .d = 0,
 .p = 1,
 .refdir[0] = '\0',
 .qrydir[0] = '\0',
 .outdir = "./",
 .num_remaining_args = 0,
 .remaining_args = NULL
};
static error_t parse_composite(int key, char* arg, struct argp_state* state) {
  struct arg_composite* composite = state->input;
  assert( composite );
  assert( composite->global );
  switch(key)
  {
  case 'b':
  {
   composite_opt.b = 1;
   break;
  }
  case 'i':
  {
    composite_opt.i = 1;
    break;
  }
  case 's':
  {
    composite_opt.s = atoi(arg);
    break;
  }
  case 'd':
  {
   composite_opt.d = 1;
   break;
  }
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
  case ARGP_KEY_ARGS:
  {
   composite_opt.num_remaining_args = state->argc - state->next;
   composite_opt.remaining_args = state->argv + state->next;
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
 if( composite_opt.refdir[0] != '\0' ){
  if (composite_opt.qrydir[0] != '\0')
   return get_species_abundance (&composite_opt);
  else if (composite_opt.i)
   return index_abv (&composite_opt);
  else if (composite_opt.s != -1){
   if( composite_opt.s >=0 && composite_opt.s <3 && composite_opt.remaining_args >0 ) return abv_search(&composite_opt);
   else printf("\vUsage: kssd composite -r <ref> -s <0|1|2> <query.abv>\n\v");
  }
  else printf("\vUsage: kssd composite -r <ref> < mode: -q | -i | -s >\n\v");
  return 1;
 }
 else if(composite_opt.d){
  if(composite_opt.remaining_args < 1) printf("\vUsage: kssd composite -d <query.abv>\n\v");
  else return read_abv (&composite_opt);
 }
 else{
  printf("\vUsage: kssd composite -r <ref> < mode: -q | -i | -s >\n\v");
  return -1;
 }
}
int read_abv (composite_opt_t *composite_opt){
 for(int i =0; i<composite_opt->num_remaining_args;i++){
    const char *ext = strrchr(composite_opt->remaining_args[i],'.');
    if(strcmp(ext + 1, binVec_suffix) != 0){
      printf("%dth argument %s is not a .abv file, skipped\n",i,composite_opt->remaining_args[i]);
      continue;
    }
  struct stat st;
    stat(composite_opt->remaining_args[i], &st);
    binVec_t tmp_abv ;
  FILE *abvfh;;
    if(( abvfh = fopen(composite_opt->remaining_args[i],"rb") ) == NULL) err(errno, "read_abv():%s",composite_opt->remaining_args[i]);
  for(int l=0; l< st.st_size; l+=sizeof(binVec_t)) {
     fread(&tmp_abv,sizeof(binVec_t),1,abvfh);
   printf("%d\t%f\n",tmp_abv.ref_idx,tmp_abv.pct);
  }
  fclose(abvfh);
 }
 return 1;
}
float *tmp_measure;
typedef struct xny_pct { float x; float y;} xny_pct_t ;
int abv_search (composite_opt_t *composite_opt){
 char* abv_fpath = malloc(PATHLEN);
 FILE *abvfh,*namefh,*abunMtx_idxfh,*abunMtx_fh,*y_l2n_fh;
 sprintf(abv_fpath,"%s/%s.%s",composite_opt->refdir,binVec_dirname,abunMtx_name_suffix);
 if(( namefh = fopen(abv_fpath,"r") ) == NULL) err(errno, "abv_search():%s",abv_fpath);
 int abv_fn = 0; int size = 100;
 char (*tmpfname)[PATHLEN] = malloc(size * PATHLEN);
 while( fgets(tmpfname[abv_fn],PATHLEN,namefh) != NULL) {
  tmpfname[abv_fn][strcspn(tmpfname[abv_fn], "\n")] = 0;
  if (abv_fn == size -1){
   size+=100;
    tmpfname = realloc( tmpfname, size * PATHLEN);
  }
  abv_fn++;
 }
 fclose(namefh);
 double *y_l2n = malloc(sizeof(double) * abv_fn);
 sprintf(abv_fpath,"%s/%s.%s",composite_opt->refdir,binVec_dirname,y_l2n_suffix);
 if((y_l2n_fh = fopen(abv_fpath,"rb") ) == NULL) err(errno, "abv_search():%s",abv_fpath);
 fread(y_l2n,sizeof(double),abv_fn,y_l2n_fh);
 fclose(y_l2n_fh);
 struct stat st;
 sprintf(abv_fpath,"%s/%s.%s",composite_opt->refdir,binVec_dirname,abunMtx_idx_suffix);
 stat(abv_fpath, &st);
 int num_otu = st.st_size / sizeof(int);
 int *abunMtx_idx= malloc(st.st_size);
 if((abunMtx_idxfh = fopen(abv_fpath,"rb") ) == NULL) err(errno, "abv_search():%s",abv_fpath);
 fread(abunMtx_idx,sizeof(int),num_otu,abunMtx_idxfh);
 fclose(abunMtx_idxfh);
  sprintf(abv_fpath,"%s/%s.%s",composite_opt->refdir,binVec_dirname,abunMtx_suffix);
  stat(abv_fpath, &st);
 binVec_t *abunMtx = malloc(st.st_size);
 if((abunMtx_fh = fopen(abv_fpath,"rb") ) == NULL) err(errno, "abv_search():%s",abv_fpath);
 fread(abunMtx,sizeof(binVec_t),st.st_size/sizeof(binVec_t),abunMtx_fh);
 fclose(abunMtx_fh);
 xny_pct_t *xny_pct ;
 int *abv_ids = malloc(sizeof(int) * abv_fn);
#define DFLT (-2)
 tmp_measure = malloc(sizeof(float) * abv_fn);
 if(composite_opt->s == 1) xny_pct = malloc(sizeof(xny_pct_t) *abv_fn);
 for(int i =0; i<composite_opt->num_remaining_args;i++){
  const char *ext = strrchr(composite_opt->remaining_args[i],'.');
  if(strcmp(ext + 1, binVec_suffix) != 0){
   printf("%dth argument %s is not a .abv file, skipped\n",i,composite_opt->remaining_args[i]);
   continue;
  }
  int n_abv_match = 0;
  sprintf(abv_fpath,"%s/%s/%s", composite_opt->refdir,binVec_dirname,composite_opt->remaining_args[i]);
  stat(abv_fpath, &st);
  int abv_len = st.st_size / sizeof(binVec_t);
  binVec_t *tmp_abv = malloc(sizeof(binVec_t)*(abv_len+1));
  if(( abvfh = fopen(abv_fpath,"rb") ) == NULL) err(errno, "abv_search():%s",abv_fpath);
  fread(tmp_abv,sizeof(binVec_t),abv_len,abvfh);
  for(int i = 0; i<abv_fn ; i++) tmp_measure[i] = DFLT;
  if(composite_opt->s == 1) memset(xny_pct,0, sizeof(xny_pct_t) *abv_fn );
  float xl2n = 0;
  for(int d = 0; d< abv_len; d++){
   int ref_idx = tmp_abv[d].ref_idx;
   xl2n += tmp_abv[d].pct*tmp_abv[d].pct;
   for(int j = ref_idx==0?0:abunMtx_idx[ref_idx-1]; j< abunMtx_idx[ref_idx]; j++){
    if(tmp_measure[abunMtx[j].ref_idx] == DFLT) {
     tmp_measure[abunMtx[j].ref_idx] = 0;
     abv_ids[n_abv_match] = abunMtx[j].ref_idx;
     n_abv_match++;
    }
    if(composite_opt->s == 1){
     tmp_measure[abunMtx[j].ref_idx] += (float)fabs(abunMtx[j].pct - tmp_abv[d].pct);
     xny_pct[abunMtx[j].ref_idx].x += tmp_abv[d].pct;
     xny_pct[abunMtx[j].ref_idx].y += abunMtx[j].pct;
    }
    else if (composite_opt->s == 2)
     tmp_measure[abunMtx[j].ref_idx] += (abunMtx[j].pct - tmp_abv[d].pct) * (abunMtx[j].pct - tmp_abv[d].pct);
    else
     tmp_measure[abunMtx[j].ref_idx] += abunMtx[j].pct * tmp_abv[d].pct;
   }
  }
  if(composite_opt->s == 0) {
   for(int n = 0; n < n_abv_match; n++)
    tmp_measure[abv_ids[n]] = tmp_measure[abv_ids[n]] / (sqrt(xl2n)*y_l2n[abv_ids[n]]) ;
  }
  if(composite_opt->s == 1) {
   for(int n = 0; n < n_abv_match; n++)
    tmp_measure[abv_ids[n]] += (2*100 - xny_pct[abv_ids[n]].x - xny_pct[abv_ids[n]].y);
   qsort(abv_ids,n_abv_match,sizeof(int),comparator_measure);
   printf("Sample\t#OTU\tL1norm\n");
   for(int n = 0; n < n_abv_match; n++) printf("%s\t%lf\n",tmpfname[abv_ids[n]],tmp_measure[abv_ids[n]]);
  }
  else if(composite_opt->s == 2){
   qsort(abv_ids,n_abv_match,sizeof(int),comparator_measure);
   for(int n = 0; n < n_abv_match; n++) printf("%s\t%lf\n",tmpfname[abv_ids[n]],sqrt(tmp_measure[abv_ids[n]]));
  }
  else{
   printf("Sample\t#OTU\tCosineXY\n");
   qsort(abv_ids,n_abv_match,sizeof(int),comparator_measure);
   for(int n = n_abv_match -1 ; n >= 0; n--) printf("%s\t%lf\n",tmpfname[abv_ids[n]],tmp_measure[abv_ids[n]]);
  }
  fclose(abvfh);
  free(tmp_abv);
 }
 free(tmp_measure);
 free(tmpfname);
 free(abunMtx);
 free(abunMtx_idx);
 free(y_l2n);
 if(composite_opt->s == 1) free(xny_pct);
 return 1;
}
int index_abv (composite_opt_t *composite_opt){
 char* abv_dpath = malloc(PATHLEN);
 sprintf(abv_dpath,"%s/%s", composite_opt->refdir,binVec_dirname);
 DIR *abv_dir = opendir(abv_dpath);
 if(!abv_dir) err(errno, "index_abv(): %s does not exists! \n",abv_dpath);
 const char *ref_dstat_fpath = test_get_fullpath(composite_opt->refdir,co_dstat);
 if (ref_dstat_fpath == NULL) err(errno,"cannot find %s under %s ",co_dstat, composite_opt->refdir);
 FILE *ref_dstat_fp;
 if( (ref_dstat_fp = fopen(ref_dstat_fpath,"rb")) == NULL ) err(errno, "index_abv():%s",ref_dstat_fpath);
 co_dstat_t ref_dstat;
 fread( &ref_dstat, sizeof(co_dstat_t),1,ref_dstat_fp);
 fclose(ref_dstat_fp);
 int *abunMtx_count = calloc(ref_dstat.infile_num,sizeof(int));
 binVec_t **abunMtx_dynamic = malloc(sizeof(binVec_t *) * ref_dstat.infile_num);
 for(int i =0; i< ref_dstat.infile_num;i++) abunMtx_dynamic[i] = malloc(sizeof(binVec_t));
 int abv_fcnt = 0;
 FILE *abvfh,*namefh,*abunMtx_idxfh,*abunMtx_fh,*y_l2n_fh;
 struct dirent *dirent;
 binVec_t binVec_tmp;
 char* abv_fpath = malloc(PATHLEN);
 sprintf(abv_fpath,"%s/%s.%s",composite_opt->refdir,binVec_dirname,abunMtx_name_suffix);
 if(( namefh = fopen(abv_fpath,"w") ) == NULL) err(errno, "index_abv():%s",abv_fpath);
 sprintf(abv_fpath,"%s/%s.%s",composite_opt->refdir,binVec_dirname,y_l2n_suffix);
 if(( y_l2n_fh = fopen(abv_fpath,"wb") ) == NULL) err(errno, "index_abv():%s",abv_fpath);
 while (( dirent = readdir(abv_dir)) != NULL){
  const char *ext = strrchr(dirent->d_name,'.');
  if(strcmp(ext + 1, binVec_suffix) == 0){
   sprintf(abv_fpath,"%s/%s",abv_dpath,dirent->d_name);
   if(( abvfh = fopen(abv_fpath,"rb") ) == NULL) err(errno, "index_abv():%s",abv_fpath);
   double y_l2n = 0;
   struct stat st;
   stat(abv_fpath, &st);
   int abv_len = st.st_size/sizeof(binVec_t) ;
   for(int i =0; i<abv_len ;i++){
    fread(&binVec_tmp , sizeof(binVec_t), 1,abvfh);
     y_l2n+= binVec_tmp.pct*binVec_tmp.pct;
    abunMtx_dynamic[binVec_tmp.ref_idx] = realloc((binVec_t *)abunMtx_dynamic[binVec_tmp.ref_idx], sizeof(binVec_t) * (abunMtx_count[binVec_tmp.ref_idx] + 1) );
        abunMtx_dynamic[binVec_tmp.ref_idx][abunMtx_count[binVec_tmp.ref_idx]].ref_idx = abv_fcnt ;
    abunMtx_dynamic[binVec_tmp.ref_idx][abunMtx_count[binVec_tmp.ref_idx]].pct = binVec_tmp.pct;
    abunMtx_count[binVec_tmp.ref_idx]++;
   }
   fclose(abvfh);
   sprintf(abv_fpath,"%s\n",dirent->d_name);
   fwrite(abv_fpath,strlen(abv_fpath),1,namefh);
   y_l2n = sqrt(y_l2n);
   fwrite(&y_l2n,sizeof(double),1,y_l2n_fh);
   abv_fcnt++;
  }
 };
 closedir(abv_dir);
 fclose(namefh);
 fclose(y_l2n_fh);
  sprintf(abv_fpath,"%s/%s.%s",composite_opt->refdir,binVec_dirname,abunMtx_suffix);
 if(( abunMtx_fh = fopen(abv_fpath,"wb") ) == NULL) err(errno, "index_abv():%s",abv_fpath);
 for(int i = 0 ;i < ref_dstat.infile_num;i++){
  if(abunMtx_count[i])
   fwrite( abunMtx_dynamic[i], sizeof(binVec_t), abunMtx_count[i],abunMtx_fh);
 }
 fclose(abunMtx_fh);
 for(int i = 1 ;i < ref_dstat.infile_num;i++)
  abunMtx_count[i] += abunMtx_count[i-1];
 sprintf(abv_fpath,"%s/%s.%s",composite_opt->refdir,binVec_dirname,abunMtx_idx_suffix);
 if(( abunMtx_idxfh = fopen(abv_fpath,"wb") ) == NULL) err(errno, "index_abv():%s",abv_fpath);
 fwrite(abunMtx_count,sizeof(int),ref_dstat.infile_num,abunMtx_idxfh);
 fclose(abunMtx_idxfh);
 free(abv_dpath);
 free(abv_fpath);
 free(abunMtx_count);
 free(abunMtx_dynamic);
 return 1;
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
 binVec_t* binVec = malloc(ref_dstat.infile_num* sizeof(binVec_t));
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
  if(composite_opt->b){
   sprintf(tmpfname,"%s/%s",composite_opt->refdir,binVec_dirname);
   mkdir(tmpfname,0777);
   sprintf(tmpfname,"%s/%s/%s.%s",composite_opt->refdir,binVec_dirname,basename(qryname[qn]),binVec_suffix);
   if( (tmpfp = fopen(tmpfname,"wb"))==NULL) err(errno,"get_species_abundance():%s",tmpfname);
  }
  int num_pass = 0;
  float binVecsum = 0;
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
   if(composite_opt->b){
    if (ref_abund[sort_ref[i]][median_idx] > 1 && kmer_num > 7){
     binVec[num_pass].ref_idx = sort_ref[i];
     binVec[num_pass].pct = (float)lastsum/lastn;
     binVecsum += binVec[num_pass].pct;
     num_pass++;
    }
   }
   else{
    printf("%s\t%s\t%d\t%f\t%f\t%d\t%d\n",qryname[qn],refname[sort_ref[i]], kmer_num, (float)sum/kmer_num,(float)lastsum/lastn,ref_abund[sort_ref[i]][median_idx], ref_abund[sort_ref[i]][kmer_num-1]);
   }
  }
  if(composite_opt->b){
   for(int i = 0; i< num_pass;i++){
    binVec[i].pct = ( binVec[i].pct - 1)*100/(binVecsum - num_pass) ;
    printf("%d\t%f\n",binVec[i].ref_idx, binVec[i].pct);
   }
   fwrite(binVec,sizeof(binVec_t),num_pass,tmpfp);
   fclose(tmpfp);
  }
  free(sort_ref);
 }
 for (int i = 0; i < ref_dstat.infile_num; i++) free(ref_abund[i]);
  free(ref_abund) ;
 free(binVec);
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
int comparator_measure (const void *a, const void *b){
 float rtv = tmp_measure[*(int*)a] - tmp_measure[*(int*)b] ;
  if(rtv > 0) return 1;
 else if(rtv < 0) return -1;
 else return 0;
}
