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
   
#include "command_dist_wrapper.h"
#include "command_dist.h"
#include "global_basic.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <argp.h>
#include <argz.h>
#include <sys/stat.h>
#include <dirent.h>
#include <err.h>
#include <errno.h>
#include <unistd.h>
#include <stdbool.h>
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif
struct arg_dist
{
  struct arg_global* global;
  char* name;
};
static struct argp_option opt_dist[] =
{
 {"halfKmerlength",'k',"INT",0, "set half Kmer length: 2-15 [8]\v" },
 {"threadN",'p',"INT",0,"set threads number [all threads]\v"},
 {"list",'l',"file",0,"a file contain paths for all query sequences\v"},
 {"DimRdcLevel",'L',"INT",0,"Dimension Reduction Level or provide .shuf file[2]\v"},
 {"maxMemory",'m',"NUM",0,"maximal memory (in G) usage allowed\v"},
 {"LstKmerOcrs",'n',"INT",0,"Specify the Least Kmer occurence in fastq file\v"},
 {"quality",'Q',"INT",0,"Filter Kmer with lowest base quality < q (Phred).[0]\v"},
 {"reference_dir",'r',"<path>",0,"reference genome/database search against.\v"},
 {"outdir",'o',"<path>",0,"folder path for results files.\v" },
 {"neighborN_max",'N',"INT",0,"max number of nearest reference genomes.[1]\v"},
 {"mutDist_max",'D',"FLT",0,"max mutation allowed for distance output.[1]\v"},
 {"metric",'M',"0/1",0,"output metrics: 0: Jaccard/1: Containment [0]\v"},
 {"outfields",'O',"0/1/2",0,"output fields(latter includes former): Distance/Q-values/Confidence Intervels.[2]\v"},
 {"correction",333,"0/1",0,"perform correction for shared k-mer counts or not .[0]\v" },
  {"abundance",'A',0,0,"abundance estimate mode.\v"},
 {"dedup",'u',0,0,"ignore repeated k-mer in reference.\v"},
 {"keepcofile",888,0,0,"keep intermedia .co files.\v"},
 {"pipecmd",'P',"<cmd>",0,"pipe command.\v"},
 {"keepskf",777,0,0,"turn on share_kmer_ct file keep mode.[false]\v"},
 {"skf",'f',"<skfpath>",0,"share_kmer_ct file path.\v"},
 {"byread",555,0,0,"sketch the file by read[false].\v"},
  { 0 }
};
static char doc_dist[] =
  "\n"
  "The dist doc prefix."
  "\v"
  "The dist doc suffix."
  ;
dist_opt_val_t dist_opt_val =
{
.k = 8,
.p = 0,
.dr_level = 2,
.dr_file = "",
.mmry = 0,
.fmt = "mfa",
.refpath = "",
.fpath = "",
.outdir = ".",
.kmerocrs = 1,
.kmerqlty = 0,
.keepco = false,
.stage2 = false,
.num_neigb = 0,
.mut_dist_max = 1,
.metric = Jcd,
.outfields = CI,
.correction = false,
.abundance = false,
.u = false,
.pipecmd = "",
.shared_kmerpath="",
.keep_shared_kmer=false,
.byread=false,
.num_remaining_args = 0,
.remaining_args = NULL
} ;
const char outdir_name[] = "kssd_rslt" ;
static error_t parse_dist(int key, char* arg, struct argp_state* state) {
  struct arg_dist* dist = state->input;
  assert( dist );
  assert( dist->global );
  switch(key)
  {
  case 'k':
   dist_opt_val.k = atoi(arg);
   break;
  case 'p':
  {
#ifdef _OPENMP
    dist_opt_val.p = atoi(arg) ;
#else
     warnx("This version of kssd was built without OpenMP and "
          "thus does not support multi threading. Ignoring -p %d",atoi(arg));
      break;
#endif
  }
   break;
  case 'm':
  {
   double sys_mm = get_sys_mmry();
   double rqst_mm = atoi(arg);
   if ( rqst_mm > sys_mm ){
    warnx("Memory request is larger than system available %f. Ignoring -m %f",sys_mm,rqst_mm);
    dist_opt_val.mmry = sys_mm;
   }
   else
    dist_opt_val.mmry = rqst_mm ;
  }
  break;
  case 'r':
  {
   if(strlen(arg) > PATHLEN) {
        err(errno,"the list path should not longer than %d", PATHLEN);
        exit(EXIT_FAILURE);
      };
   strcpy(dist_opt_val.refpath,arg);
   break;
  }
  case 'l':
  {
   if(strlen(arg) > PATHLEN) {
    err(errno,"the list path should not longer than %d", PATHLEN);
    exit(EXIT_FAILURE);
   };
   strcpy(dist_opt_val.fpath,arg);
   break;
  }
  case 'L':
  {
   struct stat path_stat;
   if( stat(arg,&path_stat) >=0 && S_ISREG(path_stat.st_mode)){
    if(strlen(arg) < PATHLEN )
     strcpy(dist_opt_val.dr_file,arg);
    else
      err(errno,"-L argument path should not longer than %d",PATHLEN);
   }
   else{
    if ( atoi(arg) >= dist_opt_val.k - 2 || atoi(arg) < 0 )
     err(errno,"-L: dimension reduction level should never larger than Kmer length - 2,"
        " which is %d here",dist_opt_val.k - 2 );
    dist_opt_val.dr_level = atoi(arg);
   }
   break;
  }
  case 'n':
  {
   if( atoi(arg) > 7 ) {
    dist_opt_val.kmerocrs = 7;
    warnx("-n argument is larger than Max, it has been set to 7, ignorned -n %d ",atoi(arg));
   }else if( atoi(arg) < 1 ){
    dist_opt_val.kmerocrs = 1;
     warnx("-n argument is smaller than Min, it has been set to 1, ignorned -n %d ",atoi(arg));
   }
   else dist_opt_val.kmerocrs = atoi(arg);
  }
   break;
  case 'Q':
  {
   dist_opt_val.kmerqlty = atoi(arg);
  }
   break;
    case 'o':
  {
      if(strlen(arg) > PATHLEN) {
        err(errno,"the outdir path should not longer than %d", PATHLEN);
        exit(EXIT_FAILURE);
      };
      strcpy(dist_opt_val.outdir,arg);
   break;
    }
  case 'N':
  {
   dist_opt_val.num_neigb = atoi(arg);
   break;
  }
  case 'D':
  {
   dist_opt_val.mut_dist_max = atof(arg);
   break;
  }
  case 'A':
  {
   dist_opt_val.abundance = true;
   break;
  }
  case 'u':
    {
      dist_opt_val.u = true;
      break;
    }
  case 555:
  {
   dist_opt_val.byread = true;
   break;
  }
  case 'P':
  {
   if(strlen(arg) > PATHLEN) {
        err(errno,"the piped command should not longer than %d", PATHLEN);
        exit(EXIT_FAILURE);
      };
      strcpy(dist_opt_val.pipecmd,arg);
      break;
  }
  case 'M':
  {
   dist_opt_val.metric = atoi(arg) ;
   break;
  }
  case 'O':
    {
      dist_opt_val.outfields = atoi(arg) ;
      break;
    }
  case 'f':
  {
   strcpy(dist_opt_val.shared_kmerpath,arg);
   break;
  }
  case 888:
  {
   dist_opt_val.keepco = true;
   break;
  }
  case 999:
  {
   dist_opt_val.stage2 = true;
      break;
  }
  case 333:
    {
      dist_opt_val.correction = atoi(arg) ;
      break;
    }
  case 777:
  {
   dist_opt_val.keep_shared_kmer = true;
   break;
  }
  case ARGP_KEY_ARGS:
    dist_opt_val.num_remaining_args = state->argc - state->next;
        dist_opt_val.remaining_args = state->argv + state->next;
   break;
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
  {
#ifdef _OPENMP
  if(dist_opt_val.p == 0)
     dist_opt_val.p = omp_get_num_procs();
#else
  if(dist_opt_val.p == 0)
     dist_opt_val.p = 1;
#endif
    return ARGP_ERR_UNKNOWN;
  }
 }
  return 0;
}
static struct argp argp_dist =
{
  opt_dist,
  parse_dist,
  "[-r <refernce>] [<query>]",
  doc_dist
};
int cmd_dist(struct argp_state* state)
{
  struct arg_dist dist = { 0, };
  int argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  char* argv0 = argv[0];
  dist.global = state->input;
  argv[0] = malloc(strlen(state->name) + strlen(" dist") + 1);
  if(!argv[0])
    argp_failure(state, 1, ENOMEM, 0);
  sprintf(argv[0], "%s dist", state->name);
  argp_parse(&argp_dist, argc, argv, ARGP_IN_ORDER, &argc, &dist);
  free(argv[0]);
  argv[0] = argv0;
  state->next += argc - 1;
 return dist_dispatch(&dist_opt_val);
}
const char *mk_dist_rslt_dir (const char *parentdirpath, const char * outdirpath )
{
  struct stat dstat;
  const char *outfullpath = malloc(PATHLEN *sizeof(char));
  sprintf((char *)outfullpath,"%s/%s",parentdirpath,outdirpath);
  if(stat(parentdirpath, &dstat) == 0 && S_ISDIR(dstat.st_mode)){
    if( stat(outfullpath, &dstat) == 0 ){
   errno = EEXIST;
      err(errno,"%s",outfullpath);
  }
    else{
      mkdir(outfullpath,0777);
    }
  }
  else {
    mkdir(parentdirpath,0777);
    mkdir(outfullpath,0777);
  }
  return outfullpath;
};
