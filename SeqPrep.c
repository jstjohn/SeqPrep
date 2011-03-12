#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>
#include "utils.h"

#define DEF_OL2MERGE_ADAPTER (10)
#define DEF_OL2MERGE_READS (10)
#define DEF_QCUT (13)
#define DEF_MIN_MATCH_ADAPTER (0.4)
#define DEF_MIN_MATCH_READS (0.4)
#define DEF_MIN_READ_LEN (30)
#define DEF_MAX_MISMATCH_ADAPTER (0.06)
#define DEF_MAX_MISMATCH_READS (0.02)
//following primer sequences are from:
//http://intron.ccam.uchc.edu/groups/tgcore/wiki/013c0/Solexa_Library_Primer_Sequences.html
//and I validated both with grep, the first gets hits to the forward file only and the second
//gets hits to the reverse file only.
#define DEF_FORWARD_PRIMER ("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG")
#define DEF_REVERSE_PRIMER ("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
void help ( char *prog_name ) {
  fprintf(stderr, "\n\nUsage:\n%s [Required Args] [Options]\n",prog_name );
  fprintf(stderr, "Required Arguments:\n" );
  fprintf(stderr, "\t-f <first read input fastq filename>\n" );
  fprintf(stderr, "\t-r <second read input fastq filename>\n" );
  fprintf(stderr, "\t-1 <first read output fastq filename>\n" );
  fprintf(stderr, "\t-2 <second read output fastq filename>\n" );
  fprintf(stderr, "General Arguments (Optional):\n" );
  fprintf(stderr, "\t-h Display this help message and exit (also works with no args) \n" );
  fprintf(stderr, "\t-6 Input sequence is in phred+64 rather than phred+33 format, the output will still be phred+33 \n" );
  fprintf(stderr, "\t-q <Quality score cutoff for mismatches to be counted in overlap; default = %d>\n", DEF_QCUT );
  fprintf(stderr, "\t-L <Minimum length of a trimmed or merged read to print it; default = %d>\n", DEF_MIN_READ_LEN );
  fprintf(stderr, "Arguments for Adapter/Primer Trimming (Optional):\n" );
  fprintf(stderr, "\t-A <forward read primer/adapter sequence to trim as it would appear at the end of a read\n\t\t (should validate by grepping a file); default = %s>\n", DEF_FORWARD_PRIMER );
  fprintf(stderr, "\t-B <reverse read primer/adapter sequence to trim as it would appear at the end of a read\n\t\t (should validate by grepping a file); default = %s>\n", DEF_REVERSE_PRIMER );
  fprintf(stderr, "\t-O <minimum overall base pair overlap with adapter sequence to trim; default = %d>\n", DEF_OL2MERGE_ADAPTER );
  fprintf(stderr, "\t-M <maximum fraction of good quality, mismatching bases for primer/adapter overlap; default = %f>\n", DEF_MAX_MISMATCH_ADAPTER );
  fprintf(stderr, "\t-N <minimum fraction of good quality, matching bases for primer/adapter overlap; default = %f>\n", DEF_MIN_MATCH_ADAPTER );
  fprintf(stderr, "Optional Arguments for Merging:\n" );
  fprintf(stderr, "\t-s <perform merging and output the merged reads to this file>\n" );
  fprintf(stderr, "\t-o <minimum overall base pair overlap to merge two reads; default = %d>\n", DEF_OL2MERGE_READS );
  fprintf(stderr, "\t-m <minimum fraction of good quality, matching bases to overlap reads; default = %f>\n", DEF_MIN_MATCH_READS );
  fprintf(stderr, "\t-n <maximum fraction of good quality, mismatching bases to overlap reads; default = %f>\n", DEF_MAX_MISMATCH_READS );
  fprintf(stderr, "\n");
  exit( 1 );
}

int main( int argc, char* argv[] ) {
  extern char* optarg;
  bool p64 = false;
  char forward_fn[MAX_FN_LEN];
  char reverse_fn[MAX_FN_LEN];
  char forward_out_fn[MAX_FN_LEN];
  char reverse_out_fn[MAX_FN_LEN];
  char merged_out_fn[MAX_FN_LEN];
  bool do_read_merging = false;
  char forward_primer[MAX_SEQ_LEN+1];
  strcpy(forward_primer, DEF_FORWARD_PRIMER); //set default
  char forward_primer_dummy_qual[MAX_SEQ_LEN+1];
  char reverse_primer[MAX_SEQ_LEN+1];
  strcpy(reverse_primer, DEF_REVERSE_PRIMER); //set default
  char reverse_primer_dummy_qual[MAX_SEQ_LEN+1];
  int i;
  for(i=0;i<MAX_SEQ_LEN+1;i++){
    forward_primer_dummy_qual[i] = 'N';//phred score of 45
    reverse_primer_dummy_qual[i] = 'N';
  }
  int ich;
  int min_ol_adapter = DEF_OL2MERGE_ADAPTER;
  int min_ol_reads = DEF_OL2MERGE_READS;
  unsigned short int min_read_len =DEF_MIN_READ_LEN;
  float min_match_adapter_frac = DEF_MIN_MATCH_ADAPTER;
  float min_match_reads_frac = DEF_MIN_MATCH_READS;
  float max_mismatch_adapter_frac = DEF_MAX_MISMATCH_ADAPTER;
  float max_mismatch_reads_frac = DEF_MAX_MISMATCH_READS;
  unsigned short max_mismatch_adapter[MAX_SEQ_LEN+1];
  unsigned short max_mismatch_reads[MAX_SEQ_LEN+1];
  unsigned short min_match_adapter[MAX_SEQ_LEN+1];
  unsigned short min_match_reads[MAX_SEQ_LEN+1];
  char qcut = (char)DEF_QCUT+33;
  SQP sqp = SQP_init();
  /* No args - help!  */
  if ( argc == 1 ) {
    help(argv[0]);
  }
  int req_args = 0;
  while( (ich=getopt( argc, argv, "f:r:1:2:q:A:s:B:O:M:N:L:o:m:n:6h" )) != -1 ) {
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
    case 'L' :
    min_read_len = atoi(optarg);
    break;

    //OPTIONAL ADAPTER/PRIMER TRIMMING ARGUMENTS
    case 'A':
    strcpy(forward_primer, optarg);
    break;
    case 'B':
    strcpy(reverse_primer, optarg);
    break;
    case 'O':
    min_ol_adapter = atoi(optarg);
    break;
    case 'M':
    max_mismatch_adapter_frac = atof(optarg);
    break;
    case 'N':
    min_match_adapter_frac = atof(optarg);
    break;

    //OPTIONAL MERGING ARGUMENTS
    case 's' :
    do_read_merging = true;
    strcpy( merged_out_fn, optarg );
    break;
    case 'o':
    min_ol_reads = atoi(optarg);
    break;
    case 'm':
    max_mismatch_reads_frac = atof(optarg);
    break;
    case 'n':
    min_match_reads_frac = atof(optarg);
    break;


    default :
    help(argv[0]);
    }
  }
  if(req_args < 4){
    fprintf(stderr, "Missing a required argument!\n");
    help(argv[0]);
  }

  //Calculate table matching overlap length to min matches and max mismatches
  for(i=0;i<MAX_SEQ_LEN+1;i++){
    max_mismatch_reads[i] = floor(((float)i)*max_mismatch_reads_frac);
    max_mismatch_adapter[i] = floor(((float)i)*max_mismatch_adapter_frac);
    min_match_reads[i] = floor(((float)i)*min_match_reads_frac);
    min_match_adapter[i] = floor(((float)i)*min_match_adapter_frac);
  }
  //get length of forward and reverse primers
  int forward_primer_len = strlen(forward_primer);
  int reverse_primer_len = strlen(reverse_primer);


  gzFile ffq = fileOpen(forward_fn, "r");
  gzFile ffqw = fileOpen(forward_out_fn,"w");
  gzFile rfq = fileOpen(reverse_fn, "r");
  gzFile rfqw = fileOpen(forward_out_fn,"w");
  gzFile mfqw = NULL;
  if(do_read_merging)
    mfqw = fileOpen(merged_out_fn,"w");
  int fpos,rpos;
  while(next_fastqs( ffq, rfq, sqp, p64 )){ //returns false when done

    fpos = compute_ol(sqp->fseq,sqp->fqual,sqp->flen,
        forward_primer, forward_primer_dummy_qual, forward_primer_len,
        min_ol_adapter, min_match_adapter, max_mismatch_adapter,
        false, qcut);
    rpos = compute_ol(sqp->rseq,sqp->rqual,sqp->rlen,
        reverse_primer, reverse_primer_dummy_qual, reverse_primer_len,
        min_ol_adapter, min_match_adapter, max_mismatch_adapter,
        false, qcut);
    if(fpos != CODE_NOMATCH || rpos != CODE_NOMATCH){
      //check if reads are long enough to do anything with.
      if(((sqp->flen - fpos) < min_read_len) || ((sqp->rlen - rpos) < min_read_len))
        continue; //ignore these reads and move on.

      // trim adapters
      sqp->fseq[fpos] = '\0';
      sqp->fqual[fpos] = '\0';
      sqp->flen = (sqp->flen - fpos);
      sqp->rseq[rpos] = '\0';
      sqp->rqual[rpos] = '\0';
      sqp->rlen = (sqp->rlen - rpos);
      if(!do_read_merging){ //just print
        write_fastq(ffqw, sqp->fid, sqp->fseq, sqp->fqual);
        write_fastq(rfqw, sqp->rid, sqp->rseq, sqp->rqual);

      }else{ //force merge
        adapter_merge(sqp);
        write_fastq(mfqw,sqp->fid,sqp->merged_seq,sqp->merged_qual);
      }
      //we are done
      continue;
    }
    //To be here we know that there isn't significant adapter overlap
    if(do_read_merging){
      if(read_merge(sqp, min_ol_reads, min_match_reads, max_mismatch_reads, qcut)){
        //print merged output
        write_fastq(mfqw,sqp->fid,sqp->merged_seq,sqp->merged_qual);
      }else{

      }


    }else{ //just write reads to output fastqs
      write_fastq(ffqw, sqp->fid, sqp->fseq, sqp->fqual);
      write_fastq(rfqw, sqp->rid, sqp->rseq, sqp->rqual);
      continue; //done
    }
  }

  return 0;
}
