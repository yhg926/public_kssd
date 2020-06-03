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
 {"outdir",'o',"<path>",0,"specify the output directory.\v",3},
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
 char insketchpath[PATHLEN];
 char pansketchpath[PATHLEN];
 char outdir[PATHLEN];
} set_opt_t ;
set_opt_t set_opt = {
.operation = -1,
.insketchpath[0] = '\0',
.pansketchpath[0]='\0',
.outdir = "./",
};
extern const char co_dstat[];
int sketch_union();
int sketch_operate();
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
  case 'o':
  {
   strcpy(set_opt.outdir, arg);
   break;
  }
  case ARGP_KEY_ARGS:
   strcpy(set_opt.insketchpath, state->argv[state->next]);
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
  if(set_opt.operation == 0 || set_opt.operation == 1 )
   return sketch_operate() ;
  else {
   printf("set operation use : -u, -i or -s\n");
   return -1 ;
  }
 }
 else
  return -1;
}
const char skch_prefix[]="combco";
const char pan_prefix[]="pan";
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
 unsigned int comp_sz = (1 << 4*COMPONENT_SZ);
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
 unsigned int comp_sz = (1 << 4*COMPONENT_SZ);
  llong* dict = (llong*)malloc(comp_sz/8);
  unsigned int tmppanco;
 for(int c=0; c< co_dstat_pan.comp_num; c++ ){
  memset(dict,0,comp_sz/8);
  sprintf(tmppath,"%s/%s.%d",set_opt.pansketchpath, pan_prefix, c);
    if(stat(tmppath, &s) != 0) err(errno,"sketch_operate():%s",tmppath);
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
