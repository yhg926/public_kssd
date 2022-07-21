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
   
#ifndef GLOBAL_BASIC
#define GLOBAL_BASIC 
#include <stdbool.h>
#define _64MASK 0xffffffffffffffffLLU
#define BIT1MASK 0x0000000000000001LLU
#if ALPHABET == 1
 #define DEFAULT 15
 #define OBJ_ALPH 16
 #define OBJ_BITS 4
 #define BIN_SZ 4096
  extern const bool Objdist[16][16];
  #define IS_DIFF(X,Y) ( Objdist[(X)][(Y)] )
#elif ALPHABET == 2
 #define DEFAULT (-1)
 #define OBJ_BITS 5
 #define BIN_SZ 2048
 #define IS_DIFF(X,Y) ( (X) != (Y) )
#else
 #define DEFAULT (-1)
 #define OBJ_ALPH 4
 #define OBJ_BITS 2
 #define BIN_SZ 65536
 #define IS_DIFF(X,Y) ( (X) != (Y) )
#endif
#define LMAX 4096
#define PATHLEN 256
#define MCO_BUF_S 4096
#ifndef COMPONENT_SZ
#define COMPONENT_SZ 7
#endif
#ifndef CTX_SPC_USE_L
#define CTX_SPC_USE_L 8
#endif
#define CTX_DR_LMT 100
#define LD_FCTR 0.6
#define GID_ARR_SZ 16
#define BBILLION 1073741824
typedef unsigned long long int llong;
 #include <argp.h>
  struct arg_global { int verbosity; };
 void log_printf(struct arg_global* g, int level, const char* fmt, ...);
 #define ARGP_KEY_INVALID 16777219
FILE * fpathopen (const char *dpath, const char *fname, const char *mode );
double get_sys_mmry(void);
#define SWAP2 0x3333333333333333ULL
#define SWAP4 0x0F0F0F0F0F0F0F0FULL
#define SWAP8 0x00FF00FF00FF00FFULL
#define SWAP16 0x0000FFFF0000FFFFULL
#define SWAP32 0x00000000FFFFFFFFULL
static inline llong crvs64bits(llong n) {
  n = ((n >> 2 ) & SWAP2 ) | ((n & SWAP2 ) << 2 );
  n = ((n >> 4 ) & SWAP4 ) | ((n & SWAP4 ) << 4 );
  n = ((n >> 8 ) & SWAP8 ) | ((n & SWAP8 ) << 8 );
  n = ((n >> 16) & SWAP16) | ((n & SWAP16) << 16);
  n = ((n >> 32) & SWAP32) | ((n & SWAP32) << 32);
  return ~n;
}
extern const int Basemap[128];
extern const char Mapbase[];
extern const unsigned int primer[25];
llong find_lgst_primer_2pow(int w);
int nextPrime(int);
typedef struct infile_entry {
 size_t fsize;
 char* fpath;
} infile_entry_t ;
typedef struct infile_tab {
  int infile_num;
  infile_entry_t* organized_infile_tab;
} infile_tab_t ;
#define BASENAME_LEN 128
typedef struct bin_stat {
 llong est_kmc_bf_dr;
 char (*seqfilebasename)[BASENAME_LEN];
 llong AllcoMem;
} bin_stat_t;
infile_tab_t * organize_infile_list(char* list_path,int fmt_ck);
infile_tab_t * organize_infile_frm_arg (int num_remaining_args, char ** remaining_args,int fmt_ck);
bin_stat_t * get_bin_basename_stat(infile_entry_t* organized_infile_tab, int *shuffle_arr,int binsz);
typedef struct co_dirstat
{
  unsigned int shuf_id;
  bool koc;
  int kmerlen;
  int dim_rd_len;
  int comp_num;
  int infile_num;
  llong all_ctx_ct;
} co_dstat_t;
#define ACPT_FMT_SZ 7
#define FAS_FMT_SZ 4
#define FQ_FMT_SZ 2
#define CO_FMT_SZ 1
#define MCO_FMT_SZ 1
#define CMPRESS_FMT_SZ 2
extern const char
*acpt_infile_fmt[ACPT_FMT_SZ],
*fasta_fmt[FAS_FMT_SZ],
*fastq_fmt[FQ_FMT_SZ],
*co_fmt[CO_FMT_SZ],
*mco_fmt[MCO_FMT_SZ],
*compress_fmt[CMPRESS_FMT_SZ];
#include <string.h>
static inline int isCompressfile(char *fname)
{
 int ret = 0;
 for(int i=0; i < CMPRESS_FMT_SZ;i++)
 {
  int basename_len = strlen(fname) - strlen(compress_fmt[i]);
  if( strcmp( ( fname + basename_len ),compress_fmt[i]) == 0 )
   return 1;
 }
 return ret;
}
static inline int isOK_fmt_infile (char* fname, const char *test_fmt[], int test_fmt_arr_sz)
{
  int ret = 0 ;
  char suftmp[10];
 for(int i=0; i < CMPRESS_FMT_SZ;i++ ){
   int basename_len = strlen(fname) - strlen(compress_fmt[i]);
   if( strcmp( ( fname + basename_len ),compress_fmt[i]) == 0 ){
     char cp_fname[PATHLEN];
     strcpy(cp_fname, fname);
     *(cp_fname + basename_len) = '\0';
     fname = cp_fname;
   break;
   };
 };
  for(int i=0; i< test_fmt_arr_sz; i++){
    sprintf(suftmp,".%s",test_fmt[i]);
    if ( strcmp((char *)(fname+strlen(fname) - strlen(suftmp)), suftmp) == 0 ){
      return 1;
    }
  };
  return ret;
};
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
static inline void check (int test, const char * message, ...)
{
    if (test) {
        va_list args;
        va_start (args, message);
        vfprintf (stderr, message, args);
        va_end (args);
        fprintf (stderr, "\n");
        exit (EXIT_FAILURE);
    }
};
typedef struct mmpco
{
 size_t fsize;
 unsigned int *mmpco;
} mmp_uint_t;
static inline mmp_uint_t mmp_uint_arr (char *cofname)
{
 mmp_uint_t cofilemmp;
 int fd;
 struct stat s;
 fd = open (cofname, O_RDONLY);
 check (fd < 0, "open %s failed: %s", cofname, strerror (errno));
 fstat (fd, & s);
 cofilemmp.fsize = s.st_size;
 cofilemmp.mmpco = mmap(NULL, s.st_size , PROT_READ, MAP_PRIVATE, fd, 0);
 check ( cofilemmp.mmpco == MAP_FAILED, "mmap %s failed: %s", cofname, strerror (errno));
 close(fd);
 return cofilemmp;
};
int str_suffix_match(char *str, const char *suf);
const char * get_pathname(const char *fullpath, const char *suf);
const char* test_get_fullpath(const char *parent_path, const char *dstat_f);
typedef struct
{
  int fasta;
  int fastq;
  int co;
 int mco;
} infile_fmt_count_t ;
infile_fmt_count_t *infile_fmt_count ( infile_tab_t * infile_tab );
extern const char co_dstat[];
extern const char skch_prefix[];
extern const char idx_prefix[];
extern const char pan_prefix[];
extern const char uniq_pan_prefix[];
typedef unsigned int ctx_obj_ct_t;
#define H1(K,HASH_SZ) ((K)%(HASH_SZ))
#define H2(K,HASH_SZ) ( 1 + (K) % ( (HASH_SZ) - 1 ) )
#define HASH(K,I,HASH_SZ) ( ( H1(K,HASH_SZ) + I * H2(K,HASH_SZ) ) % HASH_SZ )
#endif
