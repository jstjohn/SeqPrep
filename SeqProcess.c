#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdbool.h>
#include <unistd.h>
#include "utils.h"

#define DEF_OL2MERGE (10)
#define DEF_QCUT (20)
void help ( void ) {
  printf( "merge_ol_seqs -f <forward fastq filename>\n" );
  printf( "              -r <reverse fastq filename>\n" );
  printf( "              -o <base pair overlap to merge; default = %d>\n",
	  DEF_OL2MERGE );
  printf( "              -q <Quality score cutoff for mismatches in overlap; default = %d>\n",
	  DEF_QCUT );
  exit( 0 );
}

int main( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;
  char forward_fn[MAX_FN_LEN];
  char reverse_fn[MAX_FN_LEN];
  int ich;
  int ol_req = DEF_OL2MERGE;
  int qcut   = DEF_QCUT;
  SQP sqp;
  /* No args - help!  */
  if ( argc == 1 ) {
    help();
  }

  while( (ich=getopt( argc, argv, "f:r:o:q:" )) != -1 ) {
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
    default :
      help();
    }
  }
