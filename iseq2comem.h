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
   
#ifndef SEQ2CO_H
#define SEQ2CO_H 
#include "global_basic.h"
#define READSEQ_BUFFSZ 65536
#define OCCRC_BIT 16
#define OCCRC_MAX 0xffffLLU
extern void seq2co_global_var_initial(void);
llong * mmpfasta2co(char* seqfname, llong *co);
llong * fasta2co(char* seqfname,llong *co,char * pipecmd);
llong * fastq2co(char* seqfname, llong *co, char * pipecmd,int Q, int M );
llong * fastq2koc (char* seqfname, llong *co, char * pipecmd, int Q);
llong * mt_shortreads2koc (char* seqfname, llong *co, char *pipecmd,int p);
llong write_fqco2file(char* cofilename, llong *co);
llong wrt_co2cmpn_use_inn_subctx(char* cofilename, llong *co);
llong writeco2file(char* cofilename, llong *co);
llong write_fqkoc2file(char* cofilename, llong *co);
unsigned int write_fqkoc2files(char* cofilename, llong *co);
int reads2mco(char* seqfname,const char *co_dir, char * pipecmd);
#endif
