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
   
#include "command_set.h"
#include "global_basic.h"
#include "command_dist.h"
#include <sys/stat.h>
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
struct arg_set
{
  struct arg_global* global;
  char* name;
};
static struct argp_option opt_set[] =
{
  {"union",'u', 0, 0, "get union set of the sketches.\v",1 },
 {"subtract",'s',"<pan>", 0,"subtract the pan-sketch from all input sketches.\v",2 },
 {"intsect",'i',"<pan>", 0, "intersect with the pan-sketch for all input sketches.\v",2},
 {"uniq_union",'q',0, 0, "get uniq union set of the sketches.\v",3 },
 {"combin_pan",'c',0, 0, "combine pan files to combco file.\v",4 },
 {"threads",'p',"<INT>", 0, "number of threads.\v",4 },
 {"print",'P',0, 0, "print genome names.\v",4 },
 {"grouping",'g',"<file.tsv>",0,"grouping genomes by input category file.\v",4},
 {"outdir",'o',"<path>",0,"specify the output directory.\v",5},
  { 0 }
};
static char doc_set[] =
  "\n"
  "The set doc prefix."
  "\v"
  "The set doc suffix."
  ;
typedef struct set_opt
{
 int operation;
 int p;
 int P;
 int num_remaining_args;
 char ** remaining_args;
 char insketchpath[PATHLEN];
 char pansketchpath[PATHLEN];
 char subsetf[PATHLEN];
 char outdir[PATHLEN];
} set_opt_t ;
set_opt_t set_opt = {
.operation = -1,
.p = 1,
.P = 0,
.num_remaining_args = 0,
.remaining_args = NULL,
.insketchpath[0] = '\0',
.pansketchpath[0]='\0',
.subsetf[0] = '\0',
.outdir = "./"
};
typedef struct subset
{
 int taxid;
 char *taxname;
 int *gids;
} subset_t ;
typedef struct compan
{
 int taxn;
 int gn;
 subset_t * tax;
} compan_t;
int sketch_union();
int sketch_operate();
int uniq_sketch_union();
int combin_pans();
void print_gnames();
compan_t *organize_taxf(char* taxfile);
int grouping_genomes(char* taxfile) ;
static error_t parse_set(int key, char* arg, struct argp_state* state) {
  struct arg_set* set = state->input;
  assert( set );
  assert( set->global );
  switch(key)
  {
    case 'u':
  {
   if (set_opt.operation != -1 ) printf("set operation is already set, -u is ignored.\n");
   else set_opt.operation = 2 ;
   break;
  }
  case 'q':
    {
      if (set_opt.operation != -1 ) printf("set operation is already set, -q is ignored.\n");
      else set_opt.operation = 3 ;
      break;
    }
  case 's':
  {
   if (set_opt.operation != -1) printf("set operation is already set, -s is ignored.\n");
   else {
    set_opt.operation = 0;
    strcpy(set_opt.pansketchpath, arg);
   }
   break;
  }
  case 'i':
  {
   if (set_opt.operation != -1 ) printf("set operation is already set, -i is ignored.\n");
   else {
    set_opt.operation = 1 ;
    strcpy(set_opt.pansketchpath, arg);
   }
   break;
  }
  case 'c':
  {
   if (set_opt.operation != -1 ) printf("set operation is already set, -c is ignored.\n");
   else set_opt.operation = 4 ;
   break;
  }
  case 'o':
  {
   strcpy(set_opt.outdir, arg);
   break;
  }
  case 'p':
  {
   set_opt.p = atoi(arg);
   break;
  }
  case 'P':
  {
   set_opt.P = 1;
   break;
  }
  case 'g':
  {
   strcpy(set_opt.subsetf, arg);
   break;
  }
  case ARGP_KEY_ARGS:
   strcpy(set_opt.insketchpath, state->argv[state->next]);
   set_opt.num_remaining_args = state->argc - state->next;
   set_opt.remaining_args = state->argv + state->next;
   break;
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
static struct argp argp_set =
{
  opt_set,
  parse_set,
 "<combined sketch>",
  doc_set
};
int cmd_set(struct argp_state* state)
{
  struct arg_set set = { 0, };
  int argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  char* argv0 = argv[0];
  set.global = state->input;
 argv[0] = malloc(strlen(state->name) + strlen(" set") + 1);
  if(!argv[0])
    argp_failure(state, 1, ENOMEM, 0);
 sprintf(argv[0], "%s set", state->name);
  argp_parse(&argp_set, argc, argv, ARGP_IN_ORDER, &argc, &set);
 free(argv[0]);
  argv[0] = argv0;
  state->next += argc - 1;
 if(argc >1){
  if(set_opt.operation == 2)
   return sketch_union() ;
   if(set_opt.operation == 3)
   return uniq_sketch_union() ;
  if(set_opt.operation == 4)
   return combin_pans();
  if(set_opt.operation == 0 || set_opt.operation == 1 )
   return sketch_operate() ;
  else {
   if(set_opt.P) print_gnames();
   else if (set_opt.subsetf[0]!='\0') return grouping_genomes(set_opt.subsetf);
   else printf("set operation use : -u, -q, -i or -s\n");
   return -1 ;
  }
 }
 else
  return -1;
}
const char skch_prefix[]="combco";
const char idx_prefix[]="combco.index";
const char pan_prefix[]="pan";
const char uniq_pan_prefix[]="uniq_pan";
int sketch_union()
{
 const char* co_dstat_fpath = NULL;
 char combco[20];
 char unionco[20];
 co_dstat_fpath = test_get_fullpath(set_opt.insketchpath,co_dstat);
 if(co_dstat_fpath == NULL ) err(errno,"cannot find %s under %s ",co_dstat,set_opt.insketchpath);
 FILE *co_stat_fp;
 if( ( co_stat_fp = fopen(co_dstat_fpath,"rb")) == NULL ) err(errno,"sketch_union():%s",co_dstat_fpath);
 co_dstat_t co_dstat_readin;
 fread( &co_dstat_readin, sizeof(co_dstat_t),1,co_stat_fp );
 if(co_dstat_readin.infile_num == 1){
  char inpbuff;
  printf("only 1 sketch, use %s as pan-sketch?(Y/N)\n",set_opt.insketchpath);
  scanf(" %c", &inpbuff);
  if ( (inpbuff == 'Y') || (inpbuff == 'y') ) {
   chdir(set_opt.insketchpath);
   for(int i=0 ; i < co_dstat_readin.comp_num;i++){
    sprintf(combco,"%s.%d",skch_prefix,i);
    sprintf(unionco,"%s.%d",pan_prefix,i);
    if(rename(combco,unionco) !=0) err(errno,"sketch_union()");
   }
   printf("the union directory: %s created successfully\n", set_opt.insketchpath) ;
   return 1;
  }
 }
 mkdir(set_opt.outdir,0777);
 char outpath[PATHLEN];
 sprintf(outpath,"%s/%s",set_opt.outdir,co_dstat);
 FILE *co_stat_fp2;
 if( ( co_stat_fp2 = fopen(outpath,"wb")) == NULL ) err(errno,"sketch_union():%s",outpath);
 fwrite( &co_dstat_readin,sizeof(co_dstat_t),1, co_stat_fp2 );
  fclose(co_stat_fp);
  fclose(co_stat_fp2);
 size_t comp_sz = (1LLU << 4*COMPONENT_SZ);
 llong* dict = (llong*)malloc(comp_sz/8);
 unsigned int tmpcbdco;
 for(int i=0; i < co_dstat_readin.comp_num; i++){
  memset(dict,0,comp_sz/8);
  sprintf(outpath,"%s/%s.%d",set_opt.insketchpath,skch_prefix,i);
  struct stat s;
   if(stat(outpath, &s) != 0) err(errno,"sketch_union():%s",outpath);
  size_t size = s.st_size / sizeof(unsigned int);
  if( ( co_stat_fp = fopen(outpath,"rb")) == NULL ) err(errno,"sketch_union():%s",outpath);
  for(size_t n=0; n< size ; n++){
   fread(&tmpcbdco,sizeof(unsigned int),1,co_stat_fp);
   dict[tmpcbdco/64] |= ( 0x8000000000000000LLU >> (tmpcbdco % 64) ) ;
  }
  fclose(co_stat_fp);
  sprintf(outpath,"%s/%s.%d",set_opt.outdir,pan_prefix,i);
  if( ( co_stat_fp = fopen(outpath,"wb")) == NULL ) err(errno,"sketch_union():%s",outpath);
  for(unsigned int n=0;n< comp_sz/64; n++){
   if(dict[n]){
    for(int b=0; b< 64; b++){
     if ((0x8000000000000000LLU >> b) & dict[n]){
      unsigned int var = 64*n + b ;
      fwrite(&var,sizeof(unsigned int),1,co_stat_fp);
     }
    }
   }
  }
  fclose(co_stat_fp);
 }
 free(dict);
 return 1;
}
int sketch_operate()
{
 int ret = 1;
 co_dstat_t co_dstat_pan, co_dstat_origin ;
 const char* co_dstat_fpath = NULL;
  co_dstat_fpath = test_get_fullpath(set_opt.pansketchpath,co_dstat);
 if(co_dstat_fpath == NULL ) err(errno,"cannot find %s under %s ",co_dstat,set_opt.pansketchpath);
 FILE *co_stat_fp;
  if( ( co_stat_fp = fopen(co_dstat_fpath,"rb")) == NULL ) err(errno,"sketch_operate():%s",co_dstat_fpath);
  fread( &co_dstat_pan, sizeof(co_dstat_t),1,co_stat_fp );
 fclose(co_stat_fp);
 co_dstat_fpath = test_get_fullpath(set_opt.insketchpath,co_dstat);
 if(co_dstat_fpath == NULL ) err(errno,"cannot find %s under %s ",co_dstat,set_opt.insketchpath);
 struct stat s;
 if(stat(co_dstat_fpath, &s) != 0) err(errno,"sketch_operate():%s",co_dstat_fpath);
 size_t codstat_fz = s.st_size;
 co_dstat_t *tmpmem = malloc(codstat_fz);
 if( ( co_stat_fp = fopen(co_dstat_fpath,"rb")) == NULL ) err(errno,"sketch_operate():%s",co_dstat_fpath);
 fread(tmpmem,codstat_fz, 1, co_stat_fp);
 co_dstat_origin = *tmpmem;
 if (co_dstat_pan.shuf_id != co_dstat_origin.shuf_id) err(errno,"sketcing id not match(%d Vs. %d)",co_dstat_origin.shuf_id,co_dstat_pan.shuf_id);
 fclose(co_stat_fp);
 unsigned int *tmp_ctx_ct = (void *) tmpmem + sizeof(co_dstat_t);
 memset(tmp_ctx_ct,0,co_dstat_origin.infile_num*sizeof(unsigned int)) ;
 mkdir(set_opt.outdir,0777);
 char tmppath[PATHLEN];
 size_t *fco_pos = malloc(sizeof(size_t) * (co_dstat_origin.infile_num + 1) );
 size_t *post_fco_pos = malloc(sizeof(size_t) * (co_dstat_origin.infile_num + 1) );
 post_fco_pos[0] = 0;
 size_t comp_sz = (1LLU << 4*COMPONENT_SZ);
  llong* dict = (llong*)malloc(comp_sz/8);
  unsigned int tmppanco;
 for(int c=0; c< co_dstat_pan.comp_num; c++ ){
  memset(dict,0,comp_sz/8);
  sprintf(tmppath,"%s/%s.%d",set_opt.pansketchpath, pan_prefix, c);
    if(stat(tmppath, &s) != 0){
   sprintf(tmppath,"%s/%s.%d",set_opt.pansketchpath, uniq_pan_prefix, c);
   if(stat(tmppath, &s) != 0) err(errno,"sketch_operate():%s",tmppath);
  }
    size_t size = s.st_size / sizeof(unsigned int);
  if( ( co_stat_fp = fopen(tmppath,"rb")) == NULL ) err(errno,"sketch_operate():%s",tmppath);
  for(size_t n=0; n< size ; n++){
      fread(&tmppanco,sizeof(unsigned int),1,co_stat_fp);
      dict[tmppanco/64] |= ( 0x8000000000000000LLU >> (tmppanco % 64) ) ;
    }
  fclose(co_stat_fp);
  sprintf(tmppath,"%s/%s.index.%d",set_opt.insketchpath,skch_prefix, c);
  if( ( co_stat_fp = fopen(tmppath,"rb")) == NULL ) err(errno,"sketch_operate():%s",tmppath);
  fread(fco_pos,sizeof(size_t),co_dstat_origin.infile_num + 1 ,co_stat_fp);
  fclose(co_stat_fp);
  sprintf(tmppath,"%s/%s.%d",set_opt.insketchpath,skch_prefix, c);
    if( ( co_stat_fp = fopen(tmppath,"rb")) == NULL ) err(errno,"sketch_operate():%s",tmppath);
    unsigned int *cbd_fcode_mem = malloc(fco_pos[co_dstat_origin.infile_num] * sizeof(unsigned int));
    fread(cbd_fcode_mem,sizeof(unsigned int),fco_pos[co_dstat_origin.infile_num],co_stat_fp);
    fclose(co_stat_fp);
  sprintf(tmppath,"%s/%s.%d",set_opt.outdir,skch_prefix, c);
  if( ( co_stat_fp = fopen(tmppath,"wb")) == NULL) err(errno,"sketch_operate():%s",tmppath);
  for(int i = 0 ; i < co_dstat_origin.infile_num; i++){
   post_fco_pos[i+1] = post_fco_pos[i];
   for(int n = 0; n < fco_pos[i+1] - fco_pos[i]; n++){
    if( set_opt.operation == ( (dict[ cbd_fcode_mem[ fco_pos[i] + n ]/64 ] & (0x8000000000000000LLU >> (cbd_fcode_mem[ fco_pos[i] + n ] % 64)) ) > 0 ) ){
     fwrite(cbd_fcode_mem + fco_pos[i] + n, sizeof(unsigned int), 1, co_stat_fp);
     post_fco_pos[i+1]++;
     tmp_ctx_ct[i]++;
    }
   }
  }
  fclose(co_stat_fp);
  sprintf(tmppath,"%s/%s.index.%d",set_opt.outdir,skch_prefix, c);
    if( ( co_stat_fp = fopen(tmppath,"wb")) == NULL) err(errno,"sketch_operate():%s",tmppath);
  fwrite(post_fco_pos,sizeof(size_t),co_dstat_origin.infile_num + 1,co_stat_fp);
  fclose(co_stat_fp);
 }
 sprintf(tmppath,"%s/%s",set_opt.outdir,co_dstat);
  if( ( co_stat_fp = fopen(tmppath,"wb")) == NULL) err(errno,"sketch_operate():%s",tmppath);
  fwrite(tmpmem,codstat_fz,1,co_stat_fp);
  free(tmpmem);
  fclose(co_stat_fp);
 ret = 0;
 return ret ;
}
int uniq_sketch_union()
{
  const char* co_dstat_fpath = NULL;
  char combco[20];
  char unionco[20];
  co_dstat_fpath = test_get_fullpath(set_opt.insketchpath,co_dstat);
  if(co_dstat_fpath == NULL ) err(errno,"cannot find %s under %s ",co_dstat,set_opt.insketchpath);
  FILE *co_stat_fp;
  if( ( co_stat_fp = fopen(co_dstat_fpath,"rb")) == NULL ) err(errno,"uniq_sketch_union():%s",co_dstat_fpath);
  co_dstat_t co_dstat_readin;
  fread( &co_dstat_readin, sizeof(co_dstat_t),1,co_stat_fp );
  if(co_dstat_readin.infile_num == 1){
    char inpbuff;
    printf("only 1 sketch, use %s as pan-sketch?(Y/N)\n",set_opt.insketchpath);
    scanf(" %c", &inpbuff);
    if ( (inpbuff == 'Y') || (inpbuff == 'y') ) {
      chdir(set_opt.insketchpath);
      for(int i=0 ; i < co_dstat_readin.comp_num;i++){
        sprintf(combco,"%s.%d",skch_prefix,i);
        sprintf(unionco,"%s.%d",uniq_pan_prefix,i);
        if(rename(combco,unionco) !=0) err(errno,"uniq_sketch_union()");
      }
      printf("the union directory: %s created successfully\n", set_opt.insketchpath) ;
      return 1;
    }
  }
  mkdir(set_opt.outdir,0777);
  char outpath[PATHLEN];
  sprintf(outpath,"%s/%s",set_opt.outdir,co_dstat);
  FILE *co_stat_fp2;
  if( ( co_stat_fp2 = fopen(outpath,"wb")) == NULL ) err(errno,"uniq_sketch_union():%s",outpath);
  fwrite( &co_dstat_readin,sizeof(co_dstat_t),1, co_stat_fp2 );
  fclose(co_stat_fp);
  fclose(co_stat_fp2);
  size_t comp_sz = (1LLU << 4*COMPONENT_SZ);
  llong* dict = (llong*)malloc(comp_sz/8);
 llong* dict2 = (llong*)malloc(comp_sz/8);
  unsigned int tmpcbdco;
  for(int i=0; i < co_dstat_readin.comp_num; i++){
    memset(dict,0,comp_sz/8);
  memset(dict2,~0,comp_sz/8);
    sprintf(outpath,"%s/%s.%d",set_opt.insketchpath,skch_prefix,i);
    struct stat s;
    if(stat(outpath, &s) != 0) err(errno,"uniq_sketch_union():%s",outpath);
    size_t size = s.st_size / sizeof(unsigned int);
    if( ( co_stat_fp = fopen(outpath,"rb")) == NULL ) err(errno,"uniq_sketch_union():%s",outpath);
    for(size_t n=0; n< size ; n++){
      fread(&tmpcbdco,sizeof(unsigned int),1,co_stat_fp);
   if ( dict[tmpcbdco/64] & ( 0x8000000000000000LLU >> (tmpcbdco % 64) ))
    dict2[tmpcbdco/64] &= ~( 0x8000000000000000LLU >> (tmpcbdco % 64) );
      dict[tmpcbdco/64] |= ( 0x8000000000000000LLU >> (tmpcbdco % 64) ) ;
    }
    fclose(co_stat_fp);
    sprintf(outpath,"%s/%s.%d",set_opt.outdir,uniq_pan_prefix,i);
    if( ( co_stat_fp = fopen(outpath,"wb")) == NULL ) err(errno,"uniq_sketch_union():%s",outpath);
    for(unsigned int n=0;n< comp_sz/64; n++){
      if(dict[n] & dict2[n]){
        for(int b=0; b< 64; b++){
          if ( (0x8000000000000000LLU >> b) & dict[n] & dict2[n] ){
            unsigned int var = 64*n + b ;
            fwrite(&var,sizeof(unsigned int),1,co_stat_fp);
          }
        }
      }
    }
    fclose(co_stat_fp);
  }
  free(dict);
 free(dict2);
  return 1;
}
int combin_pans()
{
  const char *pan_dstat_fpath = test_get_fullpath(set_opt.remaining_args[0], co_dstat);
 FILE * co_stat_fp = fopen(pan_dstat_fpath,"rb") ;
  if( co_stat_fp == NULL ) err(errno,"combin_pans():%s", pan_dstat_fpath);
 co_dstat_t co_dstat_one, co_dstat_it;
  fread(&co_dstat_one, sizeof(co_dstat_t), 1, co_stat_fp);
 fclose(co_stat_fp) ;
 ctx_obj_ct_t *ctx_ct = calloc( set_opt.num_remaining_args , sizeof(ctx_obj_ct_t) );
 llong all_ctx_ct = 0;
 char tmppath[PATHLEN];
 FILE** com_cofp = malloc( sizeof(FILE*) * co_dstat_one.comp_num);
 FILE** indexfp = malloc( sizeof(FILE*) * co_dstat_one.comp_num);
 size_t *index_offset = calloc(co_dstat_one.comp_num,sizeof(size_t));
 mkdir(set_opt.outdir,0777);
 for(int c = 0; c < co_dstat_one.comp_num; c++){
      sprintf(tmppath,"%s/%s.%d",set_opt.outdir,skch_prefix,c);
      com_cofp[c] = fopen(tmppath,"wb");
      if(com_cofp[c] == NULL) err(errno,"%s",tmppath);
   sprintf(tmppath,"%s/%s.%d",set_opt.outdir,idx_prefix,c);
      indexfp[c] = fopen(tmppath,"wb");
      if(indexfp[c] == NULL) err(errno,"%s",tmppath);
   fwrite(index_offset+c, sizeof(size_t),1,indexfp[c]);
  }
 struct stat file_stat;
 for(int i=0; i<set_opt.num_remaining_args;i++){
  pan_dstat_fpath = test_get_fullpath(set_opt.remaining_args[i], co_dstat);
  if( ( co_stat_fp = fopen(pan_dstat_fpath,"rb")) == NULL ) err(errno,"combin_pans():%s", pan_dstat_fpath);
  fread(&co_dstat_it, sizeof(co_dstat_t), 1, co_stat_fp);
  fclose(co_stat_fp) ;
  if( co_dstat_one.shuf_id != co_dstat_it.shuf_id )
   err(errno,"combin_pans(): %dth shuf_id: %u not match 0th shuf_id: %u\n",i, co_dstat_it.shuf_id, co_dstat_one.shuf_id);
  else if (co_dstat_one.comp_num != co_dstat_it.comp_num )
   err(errno,"combin_pans(): %dth comp_num: %u not match 0th comp_num: %u\n",i, co_dstat_it.comp_num, co_dstat_one.comp_num);
  for(int c = 0; c < co_dstat_one.comp_num; c++){
   sprintf(tmppath,"%s/%s.%d",set_opt.remaining_args[i],pan_prefix,c);
   if(stat(tmppath, &file_stat) == -1) {
    sprintf(tmppath,"%s/%s.%d",set_opt.remaining_args[i],uniq_pan_prefix,c);
    if(stat(tmppath, &file_stat) == -1)
     err(errno,"%s",tmppath);
   }
   unsigned int *tmpco = malloc(file_stat.st_size);
   FILE *pan_fp = fopen(tmppath,"rb");
   fread(tmpco,file_stat.st_size, 1, pan_fp);
   fwrite(tmpco,file_stat.st_size, 1, com_cofp[c]);
   fclose(pan_fp);
   index_offset[c] += file_stat.st_size/ sizeof(unsigned int);
   fwrite(index_offset+c, sizeof(size_t),1,indexfp[c]);
   ctx_ct[i] += file_stat.st_size/ sizeof(unsigned int);
  }
  all_ctx_ct+= ctx_ct[i];
 }
 for(int c = 0; c < co_dstat_one.comp_num; c++) {
  fclose(com_cofp[c]) ;
  fclose(indexfp[c]);
 }
 free(com_cofp);
 free(indexfp);
 free(index_offset);
 co_dstat_one.infile_num = set_opt.num_remaining_args;
 co_dstat_one.all_ctx_ct = all_ctx_ct ;
 sprintf(tmppath,"%s/%s",set_opt.outdir,co_dstat);
 if( (co_stat_fp = fopen(tmppath,"wb") ) == NULL) { err(errno,"%s",tmppath) ;} ;
 fwrite(&co_dstat_one,sizeof(co_dstat_t),1,co_stat_fp);
 fwrite(ctx_ct,sizeof(ctx_obj_ct_t), set_opt.num_remaining_args, co_stat_fp);
 for(int i=0; i<set_opt.num_remaining_args;i++)
   fwrite(set_opt.remaining_args[i],1, PATHLEN,co_stat_fp);
 fclose(co_stat_fp) ;
 free(ctx_ct);
 return 1;
}
void print_gnames(){
  const char* co_dstat_fpath = NULL;
  co_dstat_fpath = test_get_fullpath(set_opt.insketchpath,co_dstat);
  if(co_dstat_fpath == NULL ) err(errno,"cannot find %s under %s ",co_dstat,set_opt.insketchpath);
  FILE *co_stat_fp;
  if( ( co_stat_fp = fopen(co_dstat_fpath,"rb")) == NULL ) err(errno,"sketch_union():%s",co_dstat_fpath);
  co_dstat_t co_dstat_readin;
  fread( &co_dstat_readin, sizeof(co_dstat_t),1,co_stat_fp );
 ctx_obj_ct_t *tmp_ctx_ct = malloc(sizeof(ctx_obj_ct_t)*co_dstat_readin.infile_num);
 fread(tmp_ctx_ct,sizeof(ctx_obj_ct_t),co_dstat_readin.infile_num,co_stat_fp);
 char (*tmpname)[PATHLEN] = malloc( PATHLEN * co_dstat_readin.infile_num );
 fread(tmpname,PATHLEN,co_dstat_readin.infile_num,co_stat_fp);
 for(int i=0; i<co_dstat_readin.infile_num; i++ ){
  printf("%s\n",tmpname[i]);
 }
 free(tmp_ctx_ct);
  free(tmpname);
}
compan_t *organize_taxf(char* taxfile){
 FILE *tf = fopen(taxfile,"r");
 if(tf == NULL) err(errno,"%s",taxfile);
 int ln = 0;
 for (char c = getc(tf); c != EOF; c = getc(tf))
  if (c == '\n') ln++;
 rewind(tf);
#define D_TAXID (-1)
 int hashsz = nextPrime( (int)((double) ln / LD_FCTR) );
 subset_t *tmphs = malloc(hashsz * sizeof(subset_t));
 for(int i = 0 ; i < hashsz ; i++) {
  tmphs[i].taxid = D_TAXID;
  tmphs[i].taxname = NULL;
 }
 char tmpstr[PATHLEN+1];
 const char s[2] = "\t";
 int tax_count = 0 ;
 for(int i= 0; i< ln ; i++) {
  fgets(tmpstr,PATHLEN,tf);
  if( tmpstr[strlen(tmpstr)-1] != '\n')
   err(errno,"organize_taxf(): %dth line %s is not full read, exceed PATHLEN %d ",i,tmpstr,PATHLEN);
  else tmpstr[strlen(tmpstr)-1] = '\0' ;
  int taxid = atoi(strtok(tmpstr,s)) ;
  char *taxname = strtok(NULL,s);
  for(int n = 0; n<hashsz ;n++){
   int hv = HASH(taxid,n,hashsz);
   if(tmphs[hv].taxid == D_TAXID) {
    tmphs[hv].taxid = taxid;
    tax_count++;
    if(taxname != NULL){
     tmphs[hv].taxname = malloc(strlen(taxname)+1);
     strcpy(tmphs[hv].taxname,taxname);
    }
    tmphs[hv].gids = malloc(2*sizeof(int));
    tmphs[hv].gids[0] = 1;
    tmphs[hv].gids[1] = i;
    break;
   }
   else if(tmphs[hv].taxid == taxid){
    if( (tmphs[hv].taxname == NULL && taxname != NULL )
     ||(tmphs[hv].taxname != NULL && taxname == NULL)
     || (taxname != NULL && strcmp(tmphs[hv].taxname,taxname) != 0)
     ) err(errno,"organize_taxf() abort!: taxid %d has different taxnames in %dth and %dth lines",taxid,tmphs[hv].gids[1],i);
     tmphs[hv].gids[0]++;
     tmphs[hv].gids = realloc(tmphs[hv].gids,(tmphs[hv].gids[0]+1) * sizeof(int));
     tmphs[hv].gids[tmphs[hv].gids[0]] = i;
     break;
   }
  }
 }
 fclose(tf);
 compan_t * subset = malloc(sizeof(compan_t));
 subset->taxn = tax_count;
 subset->gn = ln;
 subset->tax = malloc(tax_count * sizeof(subset_t)) ;
 int it = 0;
 for(int n = 0; n <hashsz; n++){
  if(tmphs[n].taxid != D_TAXID){
   subset->tax[it] = tmphs[n];
   it++;
  }
 }
 free(tmphs);
 return subset;
}
int combin_subset_pans(char* taxfile){
 compan_t *subset = organize_taxf(taxfile);
  const char* co_dstat_fpath = NULL;
 co_dstat_fpath = test_get_fullpath(set_opt.insketchpath,co_dstat);
  if(co_dstat_fpath == NULL ) err(errno,"cannot find %s under %s ",co_dstat,set_opt.insketchpath);
 FILE *tmpfh;
 if( ( tmpfh = fopen(co_dstat_fpath,"rb")) == NULL ) err(errno,"combin_subset_pans():%s",co_dstat_fpath);
 co_dstat_t co_dstat_readin;
 fread( &co_dstat_readin, sizeof(co_dstat_t),1,tmpfh );
 fclose(tmpfh);
 if(co_dstat_readin.infile_num != subset->gn)
  err(errno,"combin_subset_pans():%s's genome number %d not matches %s's genome number %d",co_dstat_fpath,co_dstat_readin.infile_num,taxfile,subset->gn);
 free(co_dstat_fpath);
 mkdir(set_opt.outdir,0777);
 char tmppath[PATHLEN]; struct stat s; int outfn = 0; llong all_ctx_ct = 0;
 ctx_obj_ct_t *ctx_ct_list = calloc(subset->taxn,sizeof(ctx_obj_ct_t));
 size_t comp_sz = (1LLU << 4*COMPONENT_SZ);
  llong* dict = (llong*)malloc(comp_sz/8);
  size_t *outcombcoidx = malloc(sizeof(size_t)* (subset->taxn+1));
  for(int c=0; c < co_dstat_readin.comp_num; c++){
  sprintf(tmppath,"%s/%s.%d",set_opt.insketchpath,skch_prefix,c);
  if(stat(tmppath, &s) != 0) err(errno,"combin_subset_pans():%s",tmppath);
  unsigned int *tmpcombco = malloc(s.st_size);
  if((tmpfh = fopen(tmppath,"rb")) == NULL) err(errno,"combin_subset_pans():%s",tmppath);
  fread(tmpcombco,s.st_size, 1,tmpfh);
  fclose(tmpfh);
  sprintf(tmppath,"%s/%s.%d",set_opt.insketchpath,idx_prefix,c);
  if(stat(tmppath, &s) != 0) err(errno,"combin_subset_pans():%s",tmppath);
  size_t *tmpcombcoidx = malloc(s.st_size);
  if((tmpfh = fopen(tmppath,"rb")) == NULL) err(errno,"combin_subset_pans():%s",tmppath);
  fread(tmpcombcoidx,s.st_size, 1,tmpfh);
  fclose(tmpfh);
  sprintf(tmppath,"%s/%s.%d",set_opt.outdir,skch_prefix,c);
  if((tmpfh = fopen(tmppath,"wb")) == NULL) err(errno,"combin_subset_pans():%s",tmppath);
  outfn = 0;
  size_t offset = 0;
  outcombcoidx[0] = 0;
  for(int t = 0; t < subset->taxn; t++){
   if(subset->tax[t].taxid == 0) continue;
    memset(dict,0,comp_sz/8);
    for(int n = 1; n <= subset->tax[t].gids[0];n++){
     int gid = subset->tax[t].gids[n] ;
     for(size_t i= tmpcombcoidx[gid]; i < tmpcombcoidx[gid+1];i++){
      dict[tmpcombco[i]/64] |= ( 0x8000000000000000LLU >> (tmpcombco[i] % 64) ) ;
     }
    }
    for(unsigned int n=0;n< comp_sz/64; n++){
        if(dict[n]){
          for(int b=0; b< 64; b++){
            if ((0x8000000000000000LLU >> b) & dict[n]){
              unsigned int var = 64*n + b ;
              fwrite(&var,sizeof(unsigned int),1,tmpfh);
        offset++; all_ctx_ct++; ctx_ct_list[outfn]++;
            }
          }
        }
    }
    outfn++;
    outcombcoidx[outfn]= offset;
    printf("%d/%d species pangenome grouped\r",t,subset->taxn);
  }
  printf("\n");
  fclose(tmpfh);
  sprintf(tmppath,"%s/%s.%d",set_opt.outdir,idx_prefix,c);
  tmpfh = fopen(tmppath,"wb");
  if(tmpfh == NULL) err(errno,"combin_subset_pans():%s",tmppath);
  fwrite(outcombcoidx,sizeof(size_t),outfn+1,tmpfh);
  fclose(tmpfh);
  free(tmpcombco);
   free(tmpcombcoidx);
 }
 co_dstat_readin.infile_num = outfn;
 co_dstat_readin.koc = 0;
 co_dstat_readin.all_ctx_ct = all_ctx_ct;
 sprintf(tmppath,"%s/%s",set_opt.outdir,co_dstat);
 if((tmpfh = fopen(tmppath,"wb")) == NULL) err(errno,"combin_subset_pans():%s",tmppath);
 fwrite(&co_dstat_readin,sizeof(co_dstat_t),1,tmpfh);
 fwrite(ctx_ct_list,sizeof(ctx_obj_ct_t),outfn,tmpfh);
 char (*tmpfname)[PATHLEN] = malloc(outfn * PATHLEN);
 int idx = 0;
  for(int t=0;t<subset->taxn;t++){
  if(subset->tax[t].taxid != 0) {
   if(subset->tax[t].taxname != NULL)
    sprintf(tmpfname[idx],"%d_%s",subset->tax[t].taxid,subset->tax[t].taxname);
   else sprintf(tmpfname[idx],"%d",subset->tax[t].taxid);
   idx++;
  }
  free(subset->tax[t].gids);
    free(subset->tax[t].taxname);
  }
  free(subset->tax);
  free(subset);
 fwrite(tmpfname,PATHLEN,outfn,tmpfh);
 fclose(tmpfh);
 free(tmpfname);
 free(ctx_ct_list);
 free(dict);
 free(outcombcoidx);
 return 0;
}
int grouping_genomes(char* taxfile){
  compan_t *subset = organize_taxf(taxfile);
  const char* co_dstat_fpath = NULL;
  co_dstat_fpath = test_get_fullpath(set_opt.insketchpath,co_dstat);
  if(co_dstat_fpath == NULL ) err(errno,"cannot find %s under %s ",co_dstat,set_opt.insketchpath);
  FILE *tmpfh;
  if( ( tmpfh = fopen(co_dstat_fpath,"rb")) == NULL ) err(errno,"grouping_genomes():%s",co_dstat_fpath);
  co_dstat_t co_dstat_readin;
  fread( &co_dstat_readin, sizeof(co_dstat_t),1,tmpfh );
  fclose(tmpfh);
  if(co_dstat_readin.infile_num != subset->gn)
    err(errno,"grouping_genomes():%s's genome number %d not matches %s's genome number %d",co_dstat_fpath,co_dstat_readin.infile_num,taxfile,subset->gn);
  free(co_dstat_fpath);
  mkdir(set_opt.outdir,0777);
  char tmppath[PATHLEN]; struct stat s; int outfn = 0; llong all_ctx_ct = 0;
  ctx_obj_ct_t *ctx_ct_list = calloc(subset->taxn,sizeof(ctx_obj_ct_t));
 unsigned int **tax_dict_ar = malloc(subset->taxn * sizeof(unsigned int *));
 int* tax_dict_size = malloc(subset->taxn*sizeof(int));
  size_t *outcombcoidx = malloc(sizeof(size_t)* (subset->taxn+1));
  for(int c=0; c < co_dstat_readin.comp_num; c++){
    sprintf(tmppath,"%s/%s.%d",set_opt.insketchpath,skch_prefix,c);
    if(stat(tmppath, &s) != 0) err(errno,"grouping_genomes():%s",tmppath);
    unsigned int *tmpcombco = malloc(s.st_size);
    if((tmpfh = fopen(tmppath,"rb")) == NULL) err(errno,"grouping_genomes():%s",tmppath);
    fread(tmpcombco,s.st_size, 1,tmpfh);
    fclose(tmpfh);
    sprintf(tmppath,"%s/%s.%d",set_opt.insketchpath,idx_prefix,c);
    if(stat(tmppath, &s) != 0) err(errno,"grouping_genomes():%s",tmppath);
    size_t *tmpcombcoidx = malloc(s.st_size);
    if((tmpfh = fopen(tmppath,"rb")) == NULL) err(errno,"grouping_genomes():%s",tmppath);
    fread(tmpcombcoidx,s.st_size, 1,tmpfh);
    fclose(tmpfh);
    for(int t = 0; t < subset->taxn; t++){
      if(subset->tax[t].taxid == 0) continue;
    int hashsize = 0;
    for(int n = 1; n <= subset->tax[t].gids[0];n++)
     hashsize += (tmpcombcoidx[subset->tax[t].gids[n]+1] - tmpcombcoidx[subset->tax[t].gids[n]]) ;
    int primer_ind = LOG2(hashsize * 1.5);
    tax_dict_size[t] = hashsize = primer_ind > 7 ? primer[primer_ind - 7] : primer[0] ;
    tax_dict_ar[t] = calloc(hashsize,sizeof(unsigned int));
        for(int n = 1; n <= subset->tax[t].gids[0];n++){
          int gid = subset->tax[t].gids[n] ;
          for(size_t i= tmpcombcoidx[gid]; i < tmpcombcoidx[gid+1];i++){
      for (int x = 0 ; x < hashsize; x++){
       int y = HASH(tmpcombco[i],x,hashsize);
       if (tax_dict_ar[t][y] == 0){
        tax_dict_ar[t][y] = tmpcombco[i];
        break;
       }
       else if (tax_dict_ar[t][y] == tmpcombco[i])
        break;
       if (x==hashsize -1){
        printf("grouping_genomes(): hashtable overflow! taxn=%d\tgid=%d\tkmer=%u\n",t,gid,tmpcombco[i]);
       }
      }
          }
        }
  }
  sprintf(tmppath,"%s/%s.%d",set_opt.outdir,skch_prefix,c);
    if((tmpfh = fopen(tmppath,"wb")) == NULL) err(errno,"grouping_genomes():%s",tmppath);
  outfn = 0;
  size_t offset = 0;
  outcombcoidx[0] = 0;
  for(int t = 0; t < subset->taxn; t++){
   if(subset->tax[t].taxid == 0) continue;
   for (int x = 0; x < tax_dict_size[t]; x++){
    if(tax_dict_ar[t][x] != 0){
     fwrite(&tax_dict_ar[t][x],sizeof(unsigned int),1,tmpfh);
     offset++;
     all_ctx_ct++;
     ctx_ct_list[outfn]++;
    }
   }
   outfn++;
   outcombcoidx[outfn]= offset;
   printf("%d/%d species pangenome grouped\r",t,subset->taxn);
   free(tax_dict_ar[t]);
  }
    printf("\n");
    fclose(tmpfh);
    sprintf(tmppath,"%s/%s.%d",set_opt.outdir,idx_prefix,c);
    tmpfh = fopen(tmppath,"wb");
    if(tmpfh == NULL) err(errno,"grouping_genomes():%s",tmppath);
    fwrite(outcombcoidx,sizeof(size_t),outfn+1,tmpfh);
    fclose(tmpfh);
    free(tmpcombco);
    free(tmpcombcoidx);
  }
  co_dstat_readin.infile_num = outfn;
  co_dstat_readin.koc = 0;
  co_dstat_readin.all_ctx_ct = all_ctx_ct;
  sprintf(tmppath,"%s/%s",set_opt.outdir,co_dstat);
  if((tmpfh = fopen(tmppath,"wb")) == NULL) err(errno,"grouping_genomes():%s",tmppath);
  fwrite(&co_dstat_readin,sizeof(co_dstat_t),1,tmpfh);
  fwrite(ctx_ct_list,sizeof(ctx_obj_ct_t),outfn,tmpfh);
  char (*tmpfname)[PATHLEN] = malloc(outfn * PATHLEN);
  int idx = 0;
  for(int t=0;t<subset->taxn;t++){
    if(subset->tax[t].taxid != 0) {
      if(subset->tax[t].taxname != NULL)
        sprintf(tmpfname[idx],"%d_%s",subset->tax[t].taxid,subset->tax[t].taxname);
      else sprintf(tmpfname[idx],"%d",subset->tax[t].taxid);
      idx++;
    }
    free(subset->tax[t].gids);
    free(subset->tax[t].taxname);
  }
  free(subset->tax);
  free(subset);
  fwrite(tmpfname,PATHLEN,outfn,tmpfh);
  fclose(tmpfh);
  free(tmpfname);
  free(ctx_ct_list);
  free(tax_dict_ar);
 free(tax_dict_size);
  free(outcombcoidx);
  return 0;
}
