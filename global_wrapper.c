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
   
#include "global_wrapper.h"
#include "global_basic.h"
#include "command_shuffle.h"
#include "command_dist_wrapper.h"
#include "command_set.h"
#include "command_reverse.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <argp.h>
#include <argz.h>
static struct argp_option opt_global[] = {
 {"usage",'u',0,OPTION_NO_USAGE,0},
 {"help",'?',0,OPTION_NO_USAGE,0},
 {"license",'l',0, OPTION_NO_USAGE,"license and copyright information.",0},
 { 0 }
};
static char doc_license[] =
"\n"
  "  Copyright 2019 Huiguang Yi. All Rights Reservered.\n\n"
  "Licensed under the Apache License, Version 2.0 (the \"License\");\n"
  "you may not use this file except in compliance with the License.\n"
  "You may obtain a copy of the License at\n\n"
  "	http://www.apache.org/licenses/LICENSE-2.0\n\n"
  "Unless required by applicable law or agreed to in writing, software\n"
  "distributed under the License is distributed on an \"AS IS\" BASIS,\n"
  "WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.\n"
  "See the License for the specific language governing permissions and\n"
  "limitations under the License.\n"
"\v"
;
static char doc_global[] =
"\n"
      "Example of parsing a nested command line."
"\v"
      "Supported subcommands are:\n"
"\n"
      "  shuffle	shuffle/sampling k-mer substring space.\n"
"\n"
      "  dist   	sequences sketching and distance estimation.\n"
"\n"
    "  set    	sketch union/intersection/substraction.\n"
"\n"
   "  reverse	reverse kssd sketch to k-mer set.\n"
"\n"
;
static error_t parse_global(int key, char* arg, struct argp_state* state)
{
  if(key == '?' || key=='u')
   state->name = long_domain;
 else if (key == 'l'){
  printf("%s\n",doc_license);
  return EINVAL;
 }
  else if(key==ARGP_KEY_ARG){
      assert( arg );
      state->name = domain;
      if(strcmp(arg, "shuffle") == 0) {
         cmd_shuffle(state);
      }
   else if(strcmp(arg, "set") == 0) {
    cmd_set(state);
   }
   else if(strcmp(arg, "dist") == 0)
   {
    cmd_dist(state);
   }
   else if(strcmp(arg, "reverse") == 0)
    cmd_reverse(state);
   else if(strcmp(arg, "primer") == 0)
     for(int i = 8;i<52;i++ )
       printf("%llu\n",find_lgst_primer_2pow(i));
   else {
        argp_error(state, "%s is not a valid command", arg);
      }
  }
  else if(key == ARGP_KEY_NO_ARGS){
        state->name = long_domain;
    printf("\n%s\n\n",argp_program_version);
    printf("Unit_space_size = %d\n\n", COMPONENT_SZ);
    printf("Type 'kssd --license' for license and copyright information.\n\n");
        argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
        printf("\v");
        argp_state_help(state,stdout,ARGP_HELP_POST_DOC);
        return EINVAL;
  }
  else if(key == ARGP_KEY_INVALID){
      state->name = state->argv[0] = domain;
      return ARGP_ERR_UNKNOWN;
  }
  else return ARGP_ERR_UNKNOWN;
  return 0;
}
static struct argp argp =
{
  opt_global,
  parse_global,
  "[arguments ...]",
  doc_global,
  0,
  0,
  0
};
void cmd_global(int argc, char**argv)
{
  struct arg_global global = { };
  argp_parse(&argp, argc, argv, ARGP_IN_ORDER, &argc, &global);
}
