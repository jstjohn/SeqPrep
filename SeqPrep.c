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
#define DEF_FORWARD_PRIMER ("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG")
#define DEF_REVERSE_PRIMER ("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
void help ( char *prog_name ) {
  fprintf(stderr, "\n\nUsage:\n%s [Required Args] [Options]\n",prog_name );
  fprintf(stderr, "Required Arguments:\n" );
  fprintf(stderr, "\t-f <forward fastq filename>\n" );
  fprintf(stderr, "\t-r <reverse fastq filename>\n" );
  fprintf(stderr, "\t-1 <forward output fastq filename>\n" );
  fprintf(stderr, "\t-2 <reverse output fastq filename>\n" );
  fprintf(stderr, "General Arguments (Optional):\n" );
  fprintf(stderr, "\t-h Display this help message and exit (also works with no args) \n" );
  fprintf(stderr, "\t-6 Input sequence is in phred+64 rather than phred+33 format, the output will still be phred+33 \n" );
  fprintf(stderr, "\t-q <Quality score cutoff for mismatches in overlap; default = %d>\n", DEF_QCUT );
  fprintf(stderr, "Arguments for Adapter/Primer Trimming (Optional):\n" );
  fprintf(stderr, "\t-A <forward read primer/adapter sequence to trim as it would appear at the end of a read (should validate by grepping a file); default = %d>\n", DEF_FORWARD_PRIMER );
  fprintf(stderr, "\t-B <reverse read primer/adapter sequence to trim as it would appear at the end of a read (should validate by grepping a file); default = %d>\n", DEF_REVERSE_PRIMER );
  fprintf(stderr, "\t-O <minimum overall base pair overlap with adapter sequence to trim; default = %d>\n", DEF_OL2MERGE_ADAPTER );
  fprintf(stderr, "\t-M <maximum fraction of good quality, mismatching bases for primer/adapter overlap; default = %f>\n", DEF_MAX_MISMATCH_ADAPTER );
  fprintf(stderr, "\t-N <minimum fraction of good quality, matching bases for primer/adapter overlap; default = %f>\n", DEF_MIN_MATCH_ADAPTER );
  fprintf(stderr, "Optional Arguments for Merging:\n" );
  fprintf(stderr, "\t-m <perform merging and output the merged reads to this file>\n" );
  fprintf(stderr, "\t-o <minimum overall base pair overlap to merge two reads; default = %d>\n", DEF_OL2MERGE_READS );
  fprintf(stderr, "\t-m <minimum fraction of good quality, matching bases to overlap reads; default = %f>\n", DEF_MIN_MATCH_READS );
  fprintf(stderr, "\t-n <maximum fraction of good quality, mismatching bases to overlap reads; default = %f>\n", DEF_MAX_MISMATCH_READS );
  fprintf(stderr, "\n");
  exit( 1 );
}

int main( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;
  bool p64 = false;
  char forward_fn[MAX_FN_LEN];
  char reverse_fn[MAX_FN_LEN];
  char forward_out_fn[MAX_FN_LEN];
  char reverse_out_fn[MAX_FN_LEN];
  char merged_out_fn[MAX_FN_LEN];
  char forward_primer[MAX_SEQ_LEN+1];
  char reverse_primer[MAX_SEQ_LEN+1];
  int ich;
  int ol_req = DEF_OL2MERGE;
  char qcut   = DEF_QCUT+33;
  SQP sqp = SQP_init();
  /* No args - help!  */
  if ( argc == 1 ) {
    help(argv[0]);
  }
  int req_args = 0;
  while( (ich=getopt( argc, argv, "f:r:1:2:q:A:B:O:M:N:o:m:n:6h" )) != -1 ) {
    switch( ich ) {

    //REQUIRED ARGUMENTS
    case 'f' :
    req_args ++;
    strcpy( forward_fn, optarg );
    break;
    case 'r' :
    req_args ++;
    strcpy( reverse_fn, optarg );
    break;
    case '1' :
    req_args ++;
    strcpy(forward_out_fn, optarg);
    break;
    case '2' :
    req_args ++;
    strcpy(reverse_out_fn, optarg);
    break;

    //OPTIONAL GENERAL ARGUMENTS
    case 'h' :
    help(argv[0]);
    break;
    case '6' :
    p64 = true;
    break;
    case 'q' :
    qcut = atoi(optarg)+33;
    break;

    //TODO: OPTIONAL ADAPTER/PRIMER TRIMMING ARGUMENTS

    //TODO: OPTIONAL MERGING ARGUMENTS


    default :
    help(argv[0]);
    }
  }
  if(req_args < 4){
    fprintf(stderr, "Missing a required argument!\n");
    help(argv[0]);
  }
  gzFile ffq = fileOpen(forward_fn, "r");
  gzFile rfq = fileOpen(reverse_fn, "r");
  while(next_fastqs( ffq, rfq, sqp, p64 )){
    //now first we need to find if a read has probable adapter overlap

    // either trim adapter and return sequence, or merge and print file:
    // if desired, also calculate overlap between reads and merge

  }


}
