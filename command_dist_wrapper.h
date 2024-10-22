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
   
#ifndef DIST_WRAPPER
#define DIST_WRAPPER 
#include "global_basic.h"
#include <stdbool.h>
#include <argp.h>
int cmd_dist(struct argp_state* state);
const char *mk_dist_rslt_dir (const char *parentdirpath, const char * outdirpath );
typedef enum { Jcd, Ctm } MTRIC;
typedef enum { Dst, Qv, CI} PFIELD;
typedef struct dist_opt_val
{
  int k ;
  int p ;
 int dr_level;
 char dr_file[PATHLEN];
  double mmry;
  char fmt[10];
  char refpath[PATHLEN];
  char fpath[PATHLEN];
  char outdir[PATHLEN];
  int kmerocrs;
  int kmerqlty;
 bool keepco;
 bool stage2;
 int num_neigb;
 double mut_dist_max;
 MTRIC metric;
 PFIELD outfields;
 bool correction;
 bool abundance;
 bool u;
 char pipecmd[PATHLEN];
 char shared_kmerpath[PATHLEN];
 bool keep_shared_kmer;
 bool byread;
  int num_remaining_args;
  char **remaining_args;
} dist_opt_val_t;
#endif
