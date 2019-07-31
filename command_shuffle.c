# 1 "command_shuffle.c"
#include "command_shuffle.h"
#include "global_basic.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <argp.h>
#include <argz.h>
#include <err.h>
#include <errno.h>
#include <math.h>

struct arg_shuffle
{
  struct arg_global* global;

  char* name;
};

static struct argp_option opt_shuffle[] =
{

  { "genomeSize", 'g', "K/M/G/T",0,"Genome size scale ( with K/M/G). K for single or few genes, M for prokayote genome, G for mammals genome. If specified, will be used for determing context length.[M]\v", 1},
  { "halfKmerLen",'k', "INT", 0, "the length of one side flanking sequence of object, both sides flanking sequence of which catenated in the context.  K should be long enough be nearly unique in genome and is suggested as short as possilbe to keep a good sensitivity, given the uniquess is satisfied, the experical formula is FORMAULA, for proyakat genome,or -g=M, k=8 is suggested, for mamals, -g=G, 11 is suggested.[8]\v" },
 {"halfSubctxLen",'s',"INT", 0,"the length of one side flanking sequence of object, both sides flanking sequence of which catenated in the subcontext. so assumming comparing dna sequences, which have an alphbet size of 4, the dimension of the whole subcontext space is 4^(2*NUM), or 16^NUM, namely, the size of population from which random samples of subspaces are drawed for the purpose of reduced representative of the orignal sequence or for decomposing orignal sequence into blocks, each block is a reduced representative of the orignal sequence using different subcontext space, subcontext should never longer than the context length presetted.[3]\v" },
 {"level",'l',"INT", 0, "the level of dimensionality reduction, this option is used the control the dimensionality reduction rate. The  dimensionality reduction rate is 16^NUM (assuming dna sequence used), e.g. set --reductionLevel=1, the reduction rate will be 16.  The reduction rate is the ratio of the dimension of the whole subcontext space divided by the dimension of drawed subspace, or the ratio of dimension of whole context space divided by dimension of context space after dimensionality reduction, which determine the memory occupy. This value should never larger than the value feed for option -s. To keep the statistical stability of dimensionality reduction, it is suggested the subcontext space dimension after reduction are still large enough( >=256 ), the pragram automatically alert when dimension after reduction are too small.[0]\v",7},
 {"outfile",'o',"STRING",0,"specify the output file name prefix, if not specify default shuffle named 'default.shuf generated'\v"},

 {"usedefault",999,0,0,"All options use default value, which assuming prokaryote genomes, halfCtxLen=8, halfSubctxLen=3, and dimensionality reduction level=2.\v",8},
  { 0 }
};

static char doc_shuffle[] =
  "\n"
  "The shuffle doc prefix."
  "\v"
  "The shuffle doc suffix."
  ;

dim_shuffle_stat_t dim_shuffle_stat = {
0,
8,
3,
2,
};

char shuf_out_file_prefix[PATHLEN]="./default";

static error_t parse_shuffle(int key, char* arg, struct argp_state* state) {
  struct arg_shuffle* shuffle = state->input;

  assert( shuffle );
  assert( shuffle->global );

  switch(key)
  {
    case 'k':
  {
   dim_shuffle_stat.k = atoi(arg);
   break;
  }
  case 's':
  {
   dim_shuffle_stat.subk = atoi(arg);
   break;
  }
  case 'l':
  {
   dim_shuffle_stat.drlevel = atoi(arg);
   break;
  }
  case 'o':
  {
   strcpy(shuf_out_file_prefix, arg);
   break;
  }
    case ARGP_KEY_NO_ARGS:
    {
   if(state->argc<2)
   {
       printf("\v");
    argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
    printf("\v");
       argp_state_help(state,stdout,ARGP_HELP_LONG);
       printf("\v");
       return EINVAL;
   }
    }
    break;
  case 999:
  {
   printf("use default values for all options\n");
  }
  break;
    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp_shuffle =
{
  opt_shuffle,
  parse_shuffle,
 0,
  doc_shuffle
};

int cmd_shuffle(struct argp_state* state)
{
  struct arg_shuffle shuffle = { 0, };
  int argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  char* argv0 = argv[0];

  shuffle.global = state->input;
 argv[0] = malloc(strlen(state->name) + strlen(" shuffle") + 1);

  if(!argv[0])
    argp_failure(state, 1, ENOMEM, 0);
 sprintf(argv[0], "%s shuffle", state->name);
  argp_parse(&argp_shuffle, argc, argv, ARGP_IN_ORDER, &argc, &shuffle);

  free(argv[0]);
  argv[0] = argv0;
  state->next += argc - 1;


  return write_dim_shuffle_file(&dim_shuffle_stat, shuf_out_file_prefix) ;
}



int * shuffle( int arr[], int len_arr)
{
 if(len_arr > RAND_MAX)
  err(errno,"shuffling array length %d longer than RAND_MAX: %d",len_arr,RAND_MAX);
 srand ( time(NULL) );
  int j, temp;
 for (int i = len_arr - 1; i >0; i-- ){
  j = rand() % (i + 1) ;
  temp = arr[i];
  arr[i] = arr[j];
  arr[j] = temp;
 }
 return arr;
};

int * shuffleN (int n, int base)
{
 int *arr;
 arr = malloc(n*sizeof(int)) ;
 if(arr==NULL) err(errno,"shuffleN");
 for(int i = 0; i < n ; i++)
   arr[i] = i + base;
 return shuffle(arr,n);
};



int add_len_drlevel2subk(void)
{
 int min_smp_len = 0, min_subctx_dim_smp_sz;
  min_subctx_dim_smp_sz = MIN_SUBCTX_DIM_SMP_SZ;
  while( min_subctx_dim_smp_sz >>= 1){ min_smp_len++; };
  return ceil((float)min_smp_len/4);
};

int write_dim_shuffle_file(dim_shuffle_stat_t* dim_shuffle_stat, char *outfile_prefix)
{
 if( dim_shuffle_stat->k < dim_shuffle_stat->subk )
  err(errno,"write_dim_shuffle_file(): half-context len: %d"
   "should larger than half-subcontext len (or dimension reduce level + 2) %d",
   dim_shuffle_stat->k,dim_shuffle_stat->subk);
 if(dim_shuffle_stat->subk >= 8)
    err(errno,"write_dim_shuffle_file(): subk shoud smaller than 8");

 int dim_after_reduction = 1<< 4*(dim_shuffle_stat->subk - dim_shuffle_stat->drlevel);
 if( dim_after_reduction < MIN_SUBCTX_DIM_SMP_SZ )
  warnx("dimension after reduction %d is smaller than the suggested minimal"
  " dimension sample size %d, which might cause loss of robutness, -s %d is suggested",dim_after_reduction,MIN_SUBCTX_DIM_SMP_SZ, dim_shuffle_stat->drlevel+2);
 if(strlen(outfile_prefix) + strlen(".shuf") >PATHLEN)
  err(errno,"output path name %s should less than %d characters",outfile_prefix,PATHLEN);

 char outfile[PATHLEN];
 sprintf(outfile,"%s.shuf",outfile_prefix);

 FILE *shuf_out;
 if( (shuf_out = fopen(outfile,"wb")) == NULL)
  err(errno,"write_dim_shuffle_file(): open file %s failed",outfile);

 srand ( time(NULL) );
 dim_shuffle_stat->id = rand();


  int *shuffled_dim = shuffleN(1 << 4*dim_shuffle_stat->subk, 0);
 shuffled_dim = shuffle(shuffled_dim,1 << 4*dim_shuffle_stat->subk);


 int ret = fwrite(dim_shuffle_stat,sizeof(dim_shuffle_stat_t),1,shuf_out)
  + fwrite(shuffled_dim,sizeof(int),1 << 4*dim_shuffle_stat->subk,shuf_out);
 fclose(shuf_out);
 free(shuffled_dim);



 return ret;
};

dim_shuffle_t* read_dim_shuffle_file(char *dim_shuffle_file)
{
 int basename_len = strlen(dim_shuffle_file) - strlen(".shuf") ;
 if( strcmp( (dim_shuffle_file + basename_len),".shuf") !=0 )
  err(errno,"read_dim_shuffle_file(): input file %s is not .shuf file",dim_shuffle_file);

 FILE *shuf_in;
 if( (shuf_in = fopen(dim_shuffle_file,"rb")) == NULL)
  err(errno,"read_dim_shuffle_file(): open file %s failed",dim_shuffle_file);

 dim_shuffle_t *dim_shuffle = malloc( sizeof(dim_shuffle_t) );
 fread(&(dim_shuffle->dim_shuffle_stat),sizeof(dim_shuffle_stat_t),1,shuf_in);

 int shuf_arr_len = 1<< 4*dim_shuffle->dim_shuffle_stat.subk ;

 dim_shuffle->shuffled_dim = malloc( sizeof(int)*shuf_arr_len );
 fread(dim_shuffle->shuffled_dim,sizeof(int)*shuf_arr_len,1,shuf_in);
 fclose(shuf_in);

 return dim_shuffle;
}
