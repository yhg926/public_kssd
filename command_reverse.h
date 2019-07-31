# 1 "command_reverse.h"
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
