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
   
#include "global_basic.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <dirent.h>
#include <sys/sysinfo.h>
#include <math.h>
#if ALPHABET == 1
const int Basemap[128] =
{
  [0 ... 127] = DEFAULT, ['z'] = DEFAULT, ['Z'] = DEFAULT ,
  ['a'] = 0, ['A'] = 0, ['c'] = 1, ['C'] = 1, ['g'] = 2, ['G'] = 2, ['t'] = 3, ['T'] = 3,
  ['w'] = 4, ['W'] = 4, ['s'] = 5, ['S'] = 5, ['m'] = 6, ['M'] = 6, ['k'] = 7, ['K'] = 7,
  ['r'] = 8, ['R'] = 8, ['y'] = 9, ['Y'] = 9, ['b'] = 10, ['B'] = 10, ['d'] = 11, ['D'] = 11,
  ['h'] = 12, ['H'] = 12, ['v'] = 13, ['V'] = 13, ['n'] = 14, ['N'] = 14,
};
const char Mapbase[] = {'A','C','G','T','W','S','M','K','R','Y','B','D','H','V','N','Z'};
#include <stdbool.h>
const bool Objdist[16][16] =
{
  [0 ... 15][0 ... 15] = 0,
  [0][1] = 1, [0][2] = 1, [0][3] = 1, [0][5] = 1, [0][7] = 1, [0][9] = 1, [0][10] = 1,
  [1][0] = 1, [2][0] = 1, [3][0] = 1, [5][0] = 1, [7][0] = 1, [9][0] = 1, [10][0] = 1,
  [1][2] = 1, [1][3] = 1, [1][4] = 1, [1][7] = 1, [1][8] = 1, [1][11] = 1,
  [2][1] = 1, [3][1] = 1, [4][1] = 1, [7][1] = 1, [8][1] = 1, [11][1] = 1,
  [2][3] = 1, [2][4] = 1, [2][6] = 1, [2][9] = 1, [2][12] = 1,
  [3][2] = 1, [4][2] = 1, [6][2] = 1, [9][2] = 1, [12][2] = 1,
  [3][5] = 1, [3][6] = 1, [3][8] = 1, [3][13] = 1,
  [5][3] = 1, [6][3] = 1, [8][3] = 1, [13][3] = 1,
  [4][5] = 1, [6][7] = 1, [8][9] = 1,
  [5][4] = 1, [7][6] = 1, [9][8] = 1,
};
#elif ALPHABET == 2
const int Basemap[128] =
{
  [0 ... 127] = DEFAULT, ['a'] = 0, ['A'] = 0, ['c'] = 1, ['C'] = 1, ['d'] = 2, ['D'] = 2, ['e'] = 3, ['E'] = 3,
  ['f'] = 4, ['F'] = 4, ['g'] = 5, ['G'] = 5, ['h'] = 6, ['H'] = 6, ['i'] = 7, ['I'] = 7,
  ['k'] = 8, ['K'] = 8, ['l'] = 9, ['L'] = 9, ['m'] = 10, ['M'] = 10, ['n'] = 11, ['N'] = 11,
  ['p'] = 12, ['P'] = 12, ['q'] = 13, ['Q'] = 13, ['r'] = 14, ['R'] = 14, ['s'] = 15, ['S'] = 15,
  ['t'] = 16, ['T'] = 16, ['v'] = 17, ['V'] = 17, ['w'] = 18, ['W'] = 18, ['y'] = 19, ['Y'] = 19
};
const char Mapbase[] = { 'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' };
#else
const int Basemap[128] =
{
  [0 ... 127] = DEFAULT,
  ['a'] = 0, ['A'] = 0,
  ['c'] = 1, ['C'] = 1,
  ['g'] = 2, ['G'] = 2,
  ['t'] = 3, ['T'] = 3,
};
const char Mapbase[]={'A','C','G','T'};
#endif
const unsigned int primer[25] =
{
  251, 509, 1021, 2039, 4093, 8191, 16381,
  32749, 65521, 131071, 262139, 524287,
  1048573, 2097143, 4194301, 8388593, 16777213,
  33554393, 67108859, 134217689, 268435399,
  536870909, 1073741789, 2147483647, 4294967291
};
double get_sys_mmry(void)
{
  struct sysinfo myinfo;
  double available_bytes;
  sysinfo(&myinfo) ;
 available_bytes = myinfo.mem_unit * myinfo.totalram ;
 return (available_bytes/BBILLION);
};
const char *acpt_infile_fmt[ACPT_FMT_SZ] = {
 "fna",
 "fas",
 "fasta",
 "fq",
 "fastq",
 "fa",
 "co"
};
const char *fasta_fmt[FAS_FMT_SZ] = {
 "fasta",
 "fna",
 "fas",
 "fa"
};
const char *fastq_fmt[FQ_FMT_SZ] = {
 "fq",
 "fastq"
};
const char *co_fmt[CO_FMT_SZ] = {
 "co"
};
const char *mco_fmt[MCO_FMT_SZ] = {
 "mco"
};
const char *compress_fmt[CMPRESS_FMT_SZ] = {
  ".gz",
  ".bz2"
};
void log_printf(struct arg_global* g, int level, const char* fmt, ...)
{
  va_list ap;
  FILE* f = stdout;
  if(g->verbosity < level)
    return;
  if(level == 0)
    f = stderr;
  va_start(ap, fmt);
  vfprintf(f, fmt, ap);
  va_end(ap);
}
FILE * fpathopen (const char *dpath, const char *fname, const char *mode )
{
  char *fullname = malloc(PATHLEN*sizeof(char));
  sprintf(fullname,"%s/%s",dpath,fname);
 struct stat s;
 if(! ( ( stat(dpath, &s) == 0 ) && S_ISDIR(s.st_mode) ) )
  mkdir(dpath, 0777);
  FILE *fp;
  if( (fp = fopen(fullname, mode) ) == NULL )
    err(errno,"fpathopen()::%s",fullname);
  return fp;
}
infile_tab_t * organize_infile_list(char* list_path, int fmt_ck)
{
  infile_tab_t *infile_stat = malloc(sizeof(infile_tab_t));
  int alloc_usize = 1024;
  infile_stat->organized_infile_tab = malloc( sizeof(infile_entry_t) * alloc_usize );
    struct stat path_stat;
    FILE * list;
    list = fopen(list_path,"r");
    if(!list) err(errno,"can't open file %s",list_path);
    char *buf = malloc( LMAX * sizeof(char));
    int file_num = 0;
  if(fmt_ck){
     while ( (fgets(buf,LMAX,list))!=NULL){
    while(isspace(*buf)) buf++;
       buf[strcspn(buf, "\r\n")] = 0;
       if( strlen(buf) < 1 )
         continue;
       if( strlen(buf) > PATHLEN )
         err(errno,"the input list: %s\n %dth line:  %s  has %lu characters exceed the maximal allowed length %d",
            list_path,file_num, buf,strlen(buf),PATHLEN);
       memset(&path_stat, 0, sizeof path_stat);
       stat(buf, &path_stat) ;
       if(!S_ISREG(path_stat.st_mode))
          err(errno,"%dth line: %s",file_num, buf);
    else if(!isOK_fmt_infile(buf,acpt_infile_fmt,ACPT_FMT_SZ)){
      printf ("isOK_fmt_infile(): wrong format %dth line: %s\nSupported format are:\n",file_num, buf);
      for(int i=0; acpt_infile_fmt[i]!=NULL;i++)
       printf(".%s ",acpt_infile_fmt[i]);
      printf("\n");
       err(errno,"program exit");
    }
    else {
     infile_stat->organized_infile_tab[file_num].fsize = path_stat.st_size;
     infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char) );
     strcpy(infile_stat->organized_infile_tab[file_num].fpath, buf);
     file_num++;
     if( file_num >= alloc_usize){
         alloc_usize+=alloc_usize;
          infile_stat->organized_infile_tab
                = realloc(infile_stat->organized_infile_tab, sizeof(infile_entry_t) * alloc_usize );
         }
    }
     };
  }
  else{
   while ( (fgets(buf,LMAX,list))!=NULL){
        while(isspace(*buf)) buf++;
        buf[strcspn(buf, "\r\n")] = 0;
        if( strlen(buf) < 1 )
          continue;
        if( strlen(buf) > PATHLEN )
          err(errno,"the input list: %s\n %dth line:  %s  has %lu characters exceed the maximal allowed length %d",
            list_path,file_num, buf,strlen(buf),PATHLEN);
    infile_stat->organized_infile_tab[file_num].fsize = 0;
    infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char) );
    strcpy(infile_stat->organized_infile_tab[file_num].fpath, buf);
    file_num++;
    if( file_num >= alloc_usize){
            alloc_usize+=alloc_usize;
            infile_stat->organized_infile_tab
                = realloc(infile_stat->organized_infile_tab, sizeof(infile_entry_t) * alloc_usize );
       }
  }
 }
    infile_stat->infile_num = file_num ;
    fclose(list);
  free(buf);
    return infile_stat;
};
infile_tab_t * organize_infile_frm_arg (int num_remaining_args, char ** remaining_args,int fmt_ck)
{
 infile_tab_t *infile_stat = malloc(sizeof(infile_tab_t));
 int file_num = 0;
 struct stat path_stat;
 DIR *dirp;
 struct dirent *dirent;
 char fullpath[PATHLEN];
 int alloc_usize = 1024;
 infile_stat->organized_infile_tab = malloc( sizeof(infile_entry_t) * alloc_usize );
if(fmt_ck) {
 for(int i=0;i<num_remaining_args;i++){
  stat(remaining_args[i],&path_stat);
  if( S_ISDIR(path_stat.st_mode)){
    if (( dirp = opendir(remaining_args[i]) ) == NULL)
     err(errno, "%dth argument: can't open %s",i+1, remaining_args[i] );
    while ((dirent = readdir(dirp)) != NULL){
     if(strlen(remaining_args[i]) + strlen(dirent->d_name) + 1 >PATHLEN)
      err(errno,"path: %s/%s exceed maximal path lenth %d",remaining_args[i], dirent->d_name,PATHLEN);
     sprintf(fullpath, "%s/%s", remaining_args[i], dirent->d_name);
     stat(fullpath,&path_stat);
        if(isOK_fmt_infile(fullpath,acpt_infile_fmt,ACPT_FMT_SZ)){
      infile_stat->organized_infile_tab[file_num].fsize = path_stat.st_size;
      infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char));
      sprintf(infile_stat->organized_infile_tab[file_num].fpath, "%s/%s", remaining_args[i], dirent->d_name);
      file_num++;
      if( file_num >= alloc_usize){
       alloc_usize+=alloc_usize;
       infile_stat->organized_infile_tab
        = realloc(infile_stat->organized_infile_tab,sizeof(infile_entry_t) * alloc_usize );
      }
     }
    }
    closedir(dirp);
  }
  else if(isOK_fmt_infile(remaining_args[i],acpt_infile_fmt,ACPT_FMT_SZ)){
   stat(remaining_args[i],&path_stat);
   infile_stat->organized_infile_tab[file_num].fsize = path_stat.st_size;
   infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char));
   strcpy(infile_stat->organized_infile_tab[file_num].fpath, remaining_args[i]);
   file_num++;
   if( file_num >= alloc_usize){
       alloc_usize+=alloc_usize;
        infile_stat->organized_infile_tab
     = realloc( infile_stat->organized_infile_tab, sizeof(infile_entry_t) * alloc_usize );
      }
  }
  else{
   printf ("wrong format %dth argument: %s\nSupported format are:\n",i+1,remaining_args[i]);
      for(int i=0; acpt_infile_fmt[i]!=NULL;i++)
       printf(".%s ",acpt_infile_fmt[i]);
      printf("\n");
      err(errno,"program exit");
  }
 };
}
else{
 for(int i=0;i<num_remaining_args;i++){
  infile_stat->organized_infile_tab[file_num].fsize = 0;
  infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char));
  strcpy(infile_stat->organized_infile_tab[file_num].fpath, remaining_args[i]);
  file_num++;
  if( file_num >= alloc_usize){
        alloc_usize+=alloc_usize;
        infile_stat->organized_infile_tab
          = realloc( infile_stat->organized_infile_tab, sizeof(infile_entry_t) * alloc_usize );
    }
 }
}
 infile_stat->infile_num = file_num ;
 return infile_stat;
};
infile_fmt_count_t *infile_fmt_count ( infile_tab_t * infile_tab )
{
 infile_fmt_count_t tmp_fmt_count = {0,0,0,0};
 infile_fmt_count_t* fmt_count = (infile_fmt_count_t*)malloc(sizeof(infile_fmt_count_t));
 *fmt_count = tmp_fmt_count;
 for(int i = 0; i < infile_tab->infile_num; i++ )
 {
  if( isOK_fmt_infile(infile_tab->organized_infile_tab[i].fpath,fasta_fmt,FAS_FMT_SZ) )
   fmt_count->fasta++;
  else if (isOK_fmt_infile (infile_tab->organized_infile_tab[i].fpath, fastq_fmt,FQ_FMT_SZ ) )
   fmt_count->fastq++;
  else if(isOK_fmt_infile (infile_tab->organized_infile_tab[i].fpath, co_fmt,CO_FMT_SZ) )
   fmt_count->co++;
  else if(isOK_fmt_infile (infile_tab->organized_infile_tab[i].fpath, mco_fmt,MCO_FMT_SZ))
   fmt_count->mco++;
  else if(infile_tab->organized_infile_tab[i].fsize != 0)
   err(errno,"infile_fmt_count(): %s is not accept format(.fasta,.fastq,.co)",infile_tab->organized_infile_tab[i].fpath);
 }
 return fmt_count;
}
#define GZCOMPRESS_RATE 5
bin_stat_t * get_bin_basename_stat( infile_entry_t* organized_infile_tab, int *shuffle_arr,int binsz)
{
 bin_stat_t *ret = malloc(sizeof(bin_stat_t));
 ret->seqfilebasename = malloc( binsz * BASENAME_LEN );
 llong fsize;
 ret->est_kmc_bf_dr = 0;
 char *fullname,*filename;
 char suftmp[10];
 char cp_filename[BASENAME_LEN];
  int compress_fmt_count = sizeof(compress_fmt)/sizeof(compress_fmt[0]) ;
 int acpt_infile_fmt_count = sizeof(acpt_infile_fmt)/sizeof(acpt_infile_fmt[0]);
 int basename_len;
 for(int i=0; i < binsz; i++){
  fsize = organized_infile_tab[shuffle_arr[i]].fsize;
  fullname = organized_infile_tab[shuffle_arr[i]].fpath;
  (filename = strrchr(fullname,'/') ) ? ++filename : (filename = fullname);
  if( strlen(filename) > BASENAME_LEN)
   err(errno,"input filename:%s excess %d ",filename,BASENAME_LEN);
  strcpy(cp_filename,filename);
  for(int j = 0; j < compress_fmt_count;j++ ){
   basename_len = strlen(filename) - strlen(compress_fmt[j]);
   if( strcmp( ( filename + basename_len ),compress_fmt[j]) == 0 ){
    strcpy( ret->seqfilebasename[i],filename);
    *(cp_filename + basename_len) = '\0';
    fsize *= GZCOMPRESS_RATE;
    break;
   };
  };
  for(int j = 0; j < acpt_infile_fmt_count; j++){
   sprintf(suftmp,".%s",acpt_infile_fmt[j]);
   basename_len = strlen(cp_filename) - strlen(suftmp);
   if( strcmp( (cp_filename + basename_len ), suftmp ) == 0 ){
    if( isOK_fmt_infile(cp_filename,fastq_fmt,FQ_FMT_SZ) )
     fsize = fsize / 2 ;
    else if( isOK_fmt_infile(cp_filename, co_fmt,CO_FMT_SZ ) )
     fsize = fsize / sizeof(llong) ;
    *(cp_filename + basename_len) = '\0';
    break;
   }
  };
  ret->est_kmc_bf_dr += fsize;
  strcpy(ret->seqfilebasename[i],cp_filename);
 };
 return ret;
}
int str_suffix_match(char *str, const char *suf)
{
 int ret = 0;
 if( (strlen(str) > strlen(suf)) && (strcmp( (str + strlen(str) - strlen(suf)),suf) == 0) )
  ret = 1;
 return ret;
};
const char * get_pathname(const char *fullpath, const char *suf)
{
 char *pathcp = malloc( strlen(fullpath) + 1) ;
 strcpy(pathcp,fullpath);
 *(pathcp + strlen(pathcp) - strlen(suf)) = '\0';
 return pathcp;
}
llong find_lgst_primer_2pow(int w)
{
 if( w < 2 || w > 62 ){
  perror("find_1st_primer_after_2pow: argument should between 8 and 62");
  exit(EXIT_FAILURE);
 }
  llong n = ( 1llu << w ) ;
 llong hshsz = (llong) ( n * CTX_SPC_USE_L / LD_FCTR) ;
 printf("w=%d\tspace_sz=%llu\thashsize=%llu\tkmerlimt=%llu\n",w,n,hshsz,(llong)(hshsz*LD_FCTR) ) ;
 llong i = 3, c ; llong prime = 0;
 for(i = n - 1 ; i > (n >> 1) ; i--)
 {
  for ( c = 2 ; c <= (int)pow(i+1,0.5) ; c++ )
  {
   if( i%c == 0 )
    break;
  }
  if( c*c > i ){
   prime = i;
   break;
  }
 }
 printf("nearest prime=%llu\n",prime);
 return prime;
}
