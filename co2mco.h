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
   
#ifndef CO2MCO_H
#define CO2MCO_H 
#include "global_basic.h"
#include <stdint.h>
typedef uint32_t gidobj_t;
typedef struct mco_index
{
  size_t *row_offset;
  unsigned int *row_gnum;
  unsigned int *row_bin_gnum;
} mco_index_t;
typedef struct kmerdb_index
{
  size_t *row_offset;
  unsigned int *row_gnum;
  unsigned int *row_bin_gnum;
} kmerdb_index_t;
typedef struct gid_arr_llist
{
  gidobj_t gidobj[GID_ARR_SZ];
  struct gid_arr_llist *next;
} gid_arr_llist_t ;
typedef struct
{
 gidobj_t g_num ;
 gid_arr_llist_t *next;
} mco_entry_stat_t ;
mco_entry_stat_t** co2unitllmco(const char *codirname, int bin_sz, int bin_id, int component_id);
mco_entry_stat_t** fread_co2unitllmco(const char *codirname, int bin_sz, int bin_id, int component_id);
void combco2mco(const char *mcodirname, const char *codirname, int cofnum, int comp_num, int p_fit_mem);
void cdb_kmerf2kmerdb(const char *mcodirname, const char *codirname, int cofnum, int comp_num, int p_fit_mem);
gidobj_t** llmco2arrmco(mco_entry_stat_t** llmco);
unsigned int write_unit_arrmco_file(const char* unitmcofname, gidobj_t** arrmco);
gidobj_t** read_unit_arrmco_file(const char *mco_fncode);
void free_unit_arrmco(gidobj_t** unit_arrmco);
size_t est_unitllmco_mem(void);
size_t precise_est_unitllmco_mem(const char *co_dstat_fpath);
#endif
