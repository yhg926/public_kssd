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
   
#ifndef COMMAND_COMPONENT
#define COMMAND_COMPONENT 
#include <argp.h>
#include "global_basic.h"
extern const char binVec_suffix[];
extern const char abunMtx_suffix[];
extern const char binVec_dirname[];
extern const char abunMtx_idx_suffix[];
extern const char abunMtx_name_suffix[];
extern const char y_l2n_suffix[];
typedef struct binVec
{
 int ref_idx;
 float pct;
} binVec_t;
typedef struct composite_opt
{
 int b;
 int i;
 int s;
 int d;
 int p;
 char refdir[PATHLEN];
 char qrydir[PATHLEN];
 char outdir[PATHLEN];
 int num_remaining_args;
  char **remaining_args;
} composite_opt_t;
int cmd_composite(struct argp_state* state);
int comparator_idx (const void*, const void*);
int comparator_measure (const void *a, const void *b);
int comparator (const void*, const void*);
int get_species_abundance (composite_opt_t *);
int index_abv (composite_opt_t *);
int abv_search (composite_opt_t *composite_opt);
int read_abv (composite_opt_t *composite_opt);
#endif
