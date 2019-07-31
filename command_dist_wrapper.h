# 1 "command_dist_wrapper.h"
#ifndef DIST_WRAPPER
#define DIST_WRAPPER 

#include "global_basic.h"
#include <stdbool.h>
#include <argp.h>

int cmd_dist(struct argp_state* state);
const char *mk_dist_rslt_dir (const char *parentdirpath, const char * outdirpath );

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
 bool abundance;
  int num_remaining_args;
  char **remaining_args;

} dist_opt_val_t;

#endif
