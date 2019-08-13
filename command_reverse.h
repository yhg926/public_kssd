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
   
#ifndef COMMAND_REVERSE
#define COMMAND_REVERSE 
#include "global_basic.h"
#include <stdbool.h>
typedef struct reverse_opt_val
{
 char shufile[PATHLEN];
 char outdir[PATHLEN];
 int p;
 int num_remaining_args;
  char **remaining_args;
} reverse_opt_val_t;
int co_reverse2kmer(reverse_opt_val_t *opt_val);
llong core_reverse2unituple(unsigned int kid, int compid, int compbit, int pf_bits, int inner_ctx_bits, int half_outer_ctx_bits, unsigned int *rev_shuf_arr);
#include <argp.h>
int cmd_reverse(struct argp_state* state);
#endif
