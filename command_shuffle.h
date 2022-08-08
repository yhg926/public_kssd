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
   
#ifndef COMMAND_SHUFFLE
#define COMMAND_SHUFFLE 
typedef struct dim_shuffle_stat
{
 int id;
 int k;
 int subk;
 int drlevel;
} dim_shuffle_stat_t;
typedef struct dim_shuffle
{
 dim_shuffle_stat_t dim_shuffle_stat;
 int *shuffled_dim;
} dim_shuffle_t;
#define MIN_SUBCTX_DIM_SMP_SZ 4096
int * shuffle ( int arr[], int len_arr);
int * shuffleN ( int n, int base);
int write_dim_shuffle_file(dim_shuffle_stat_t* dim_shuffle_stat, char *outfile_prefix);
dim_shuffle_t* read_dim_shuffle_file(char *dim_shuffle_file);
int add_len_drlevel2subk(void);
 #include <argp.h>
 int cmd_shuffle(struct argp_state* state);
#endif
