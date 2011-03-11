#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdbool.h>
#include <unistd.h>
#include "utils.h"

#define DEF_OL2MERGE_ADAPTER (10)
#define DEF_OL2MERGE_READS (10)
#define DEF_QCUT (20)
#define DEF_MIN_MATCH_ADAPTER (0.5)
#define DEF_MIN_MATCH_READS (0.5)
#define DEF_MIN_INFERRED_INS_LEN_TO_PRINT (30)
#define DEF_MAX_MISMATCH_ADAPTER (0.1)
#define DEF_MAX_MISMATCH_READS (0.1)
#define DEF_FORWARD_PRIMER ("")
#define DEF_REVERSE_PRIMER ("")
void help ( char *prog_name ) {
  fprintf(stderr, "\n\nUsage:\n%s [Required Args] [Options]\n",prog_name );
  fprintf(stderr, "Required Arguments:\n" );
  fprintf(stderr, "\t-f <forward fastq filename>\n" );
  fprintf(stderr, "\t-r <reverse fastq filename>\n" );
  fprintf(stderr, "\t-1 <forward output fastq filename>\n" );
  fprintf(stderr, "\t-2 <reverse output fastq filename>\n" );
  fprintf(stderr, "Optional General Arguments:\n" );
  fprintf(stderr, "\t-q <Quality score cutoff for mismatches in overlap; default = %d>\n", DEF_QCUT );
  fprintf(stderr, "Optional Arguments for Adapter/Primer Trimming:\n" );
  fprintf(stderr, "\t-q <maximum fraction of good quality, mismatching bases for pimer/adapter; default = %f>\n", DEF_MAX_MISMATCH_ADAPTER );
  fprintf(stderr, "\t-d <maximum fraction of good quality, mismatching bases for pimer/adapter; default = %f>\n", DEF_MAX_MISMATCH_ADAPTER );
  fprintf(stderr, "\t-d <maximum fraction of good quality, mismatching bases for ; default = %f>\n", DEF_MIN_MATCH_ADAPTER );
  fprintf(stderr, "Optional Arguments for Merging:\n" );
  fprintf(stderr, "\t-m <perform merging and output the merged reads to this file>\n" );
  fprintf(stderr, "\t-o <base pair overlap to merge; default = %d>\n", DEF_OL2MERGE );
  fprintf(stderr, "\t-u <min fraction of good quality, matching bases; default = %f>\n", DEF_MIN_MATCH_READS );
  fprintf(stderr, "\n");
  exit( 1 );
}

int main( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;
  bool p64 = false;
  char forward_fn[MAX_FN_LEN];
  char reverse_fn[MAX_FN_LEN];
  char foward_out_fn[MAX_FN_LEN];
  char reverse_out_fn[MAX_FN_LEN];
  char merged_out_fn[MAX_FN_LEN];
  char forward_primer[MAX_SEQ_LEN+1];
  char reverse_primer[MAX_SEQ_LEN+1];
  int ich;
  int ol_req = DEF_OL2MERGE;
  int qcut   = DEF_QCUT;
  SQP sqp = SQP_init();
  /* No args - help!  */
  if ( argc == 1 ) {
    help(argv[0]);
  }

  while( (ich=getopt( argc, argv, "f:r:1:2:m:o:q:6:h" )) != -1 ) {
    switch( ich ) {
    case 'f' :      
    strcpy( forward_fn, optarg );
    break;
    case 'r' :
    strcpy( reverse_fn, optarg );
    break;
    case 'o' :
    ol_req = atoi( optarg );
    break;
    case 'q' :
    qcut = atoi( optarg );
    break;
    case 'h' :
    help(argv[0]);
    break;
    default :
    help(argv[0]);
    }
  }

  gzFile ffq = fileOpen(forward_fn, "r");
  gzFile rfq = fileOpen(reverse_fn, "r");
  while(next_fastqs( ffq, rfq, sqp, p64 )){
    //now first we need to find if a read has probable adapter overlap
    // either trim adapter and return sequence, or merge and print file:
    // if desired, also calculate overlap between reads and merge

  }


}
