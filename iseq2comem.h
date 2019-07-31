# 1 "iseq2comem.h"
#ifndef SEQ2CO_H
#define SEQ2CO_H 
#include "global_basic.h"
#define READSEQ_BUFFSZ 65536

#define OCCRC_BIT 24
#define OCCRC_MAX 0xffffffLLU

extern void seq2co_global_var_initial(void);


llong * mmpfasta2co(char* seqfname, llong *co);
llong * fasta2co(char* seqfname,llong *co);
llong * fastq2co(char* seqfname, llong *co, int Q, int M );
llong * fastq2koc (char* seqfname, llong *co, int Q);

llong write_fqco2file(char* cofilename, llong *co);
llong wrt_co2cmpn_use_inn_subctx(char* cofilename, llong *co);
llong writeco2file(char* cofilename, llong *co);
llong write_fqkoc2file(char* cofilename, llong *co);
#endif
