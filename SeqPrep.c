#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "utils.h"
#include "stdaln.h"

#define DEF_OL2MERGE_ADAPTER (10)
#define DEF_OL2MERGE_READS (20)
#define DEF_QCUT (10)
#define DEF_MIN_MATCH_ADAPTER (0.7)
#define DEF_MIN_MATCH_READS (0.75)
#define DEF_MIN_READ_LEN (30)
#define DEF_MAX_MISMATCH_ADAPTER (0.13)
#define DEF_MAX_MISMATCH_READS (0.02)
#define DEF_MAX_PRETTY_PRINT (10000)
#define DEF_ADAPTER_SCORE_THRES (18)
#define DEF_READ_SCORE_THRES (-500)
//two revolutions of 4 positions = 5000 reads
#define SPIN_INTERVAL (1250)
//following primer sequences are from:
//http://intron.ccam.uchc.edu/groups/tgcore/wiki/013c0/Solexa_Library_Primer_Sequences.html
//and I validated both with grep, the first gets hits to the forward file only and the second
//gets hits to the reverse file only.
#define DEF_FORWARD_PRIMER ("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG")
#define DEF_REVERSE_PRIMER ("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")

void help ( char *prog_name ) {
  fprintf(stderr, "\n\nUsage:\n%s [Required Args] [Options]\n",prog_name );
  fprintf(stderr, "NOTE 1: The output is always gziped compressed.\n");
  fprintf(stderr, "NOTE 2: If the quality strings in the output contain characters less than ascii 33 on an ascii table (they look like lines from a binary file), try running again with or without the -6 option.\n");
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
  fprintf(stderr, "\t-M <maximum fraction of good quality mismatching bases for primer/adapter overlap; default = %f>\n", DEF_MAX_MISMATCH_ADAPTER );
  fprintf(stderr, "\t-N <minimum fraction of matching bases for primer/adapter overlap; default = %f>\n", DEF_MIN_MATCH_ADAPTER );
  fprintf(stderr, "\t-b <adapter alignment band-width; default = %d>\n", aln_param_nt2nt.band_width );
  fprintf(stderr, "\t-Q <adapter alignment gap-open; default = %d>\n", aln_param_nt2nt.gap_open );
  fprintf(stderr, "\t-t <adapter alignment gap-extension; default = %d>\n", aln_param_nt2nt.gap_ext );
  fprintf(stderr, "\t-e <adapter alignment gap-end; default = %d>\n", aln_param_nt2nt.gap_end );
  fprintf(stderr, "\t-Z <adapter local alignment cutoff score ((2*num_matches) - (gap_open*num_gaps) - (gap_close*num_gaps) - (gap_ext*gap_len)) ; default = %d>\n", DEF_ADAPTER_SCORE_THRES );
  fprintf(stderr, "Optional Arguments for Merging:\n" );
  fprintf(stderr, "\t-g <print overhang when adapters are present and stripped (use this if reads are different length)>\n");
  fprintf(stderr, "\t-s <perform merging and output the merged reads to this file>\n" );
  fprintf(stderr, "\t-E <write pretty alignments to this file for visual Examination>\n" );
  fprintf(stderr, "\t-x <max number of pretty alignments to write (if -E provided); default = %d>\n", DEF_MAX_PRETTY_PRINT );
  fprintf(stderr, "\t-o <minimum overall base pair overlap to merge two reads; default = %d>\n", DEF_OL2MERGE_READS );
  fprintf(stderr, "\t-m <minimum fraction of matching bases to overlap reads; default = %f>\n", DEF_MIN_MATCH_READS );
  fprintf(stderr, "\t-n <maximum fraction of good quality mismatching bases to overlap reads; default = %f>\n", DEF_MAX_MISMATCH_READS );
  fprintf(stderr, "\n");
  exit( 1 );
}

static unsigned short spcount = 0;
/**
 * Have a nice spinner to give you a false sense of hope
 */
inline void update_spinner(unsigned long long num_reads){
  if(num_reads == 0){
    fprintf(stderr,"Processing reads... |");
    fflush(stderr);
  }else if (num_reads % SPIN_INTERVAL == 0){
    switch(spcount % 4){
    case 0:
      fprintf(stderr,"\b/");
      fflush(stderr);
      break;
    case 1:
      fprintf(stderr,"\b-");
      fflush(stderr);
      break;
    case 2:
      fprintf(stderr,"\b\\");
      fflush(stderr);
      break;
    default:
      fprintf(stderr,"\b|");
      fflush(stderr);
      break;
    }
    spcount++;
  }
}


int main( int argc, char* argv[] ) {
  unsigned long long num_pairs;
  unsigned long long num_merged;
  unsigned long long num_adapter;
  unsigned long long num_discarded;
  unsigned long long num_too_ambiguous_to_merge;
  unsigned long long max_pretty_print = DEF_MAX_PRETTY_PRINT;
  unsigned long long num_pretty_print = 0;
  unsigned int adapter_thresh = DEF_ADAPTER_SCORE_THRES;
  unsigned int read_thresh = DEF_READ_SCORE_THRES;
  clock_t start, end;
  //init to 0
  num_pairs = num_merged = num_adapter = num_discarded = num_too_ambiguous_to_merge = 0;
  extern char* optarg;
  bool p64 = false;
  char forward_fn[MAX_FN_LEN];
  char reverse_fn[MAX_FN_LEN];
  char forward_out_fn[MAX_FN_LEN];
  char reverse_out_fn[MAX_FN_LEN];
  char merged_out_fn[MAX_FN_LEN];
  bool do_read_merging = false;
  bool print_overhang = false;
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
  bool pretty_print = false;
  char pretty_print_fn[MAX_FN_LEN+1];
  SQP sqp = SQP_init();
  /* No args - help!  */
  if ( argc == 1 ) {
    help(argv[0]);
  }
  int req_args = 0;
  while( (ich=getopt( argc, argv, "f:r:1:2:q:A:s:B:O:E:x:M:N:L:o:m:b:Q:t:e:Z:n:6gh" )) != -1 ) {
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
    case 'b':
      aln_param_nt2nt.band_width = atoi(optarg);
      break;
    case 'Q':
      aln_param_nt2nt.gap_open = atoi(optarg);
      break;
    case 't':
      aln_param_nt2nt.gap_ext = atoi(optarg);
      break;
    case 'e':
      aln_param_nt2nt.gap_end = atoi(optarg);
      break;
    case 'Z':
      adapter_thresh = atoi(optarg);
      break;

      //OPTIONAL MERGING ARGUMENTS
    case 'g' :
      print_overhang = true;
      break;
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
    case 'E':
      pretty_print = true;
      strcpy(pretty_print_fn,optarg);
      break;
    case 'x':
      max_pretty_print = atol(optarg);
      break;


    default :
      help(argv[0]);
    }
  }
  if(req_args < 4){
    fprintf(stderr, "Missing a required argument!\n");
    help(argv[0]);
  }
  start = clock();
  //allocate alignment memory

  //  int min_match = 8;
  //  int ngaps = 1;
  //  int maxglen = 3;

  // AlnParam aln_param_adapter   = {  5, 13, 19, aln_sm_read, 16, 75 };
  //


  //Calculate table matching overlap length to min matches and max mismatches
  for(i=0;i<MAX_SEQ_LEN+1;i++){
    max_mismatch_reads[i] = floor(((float)i)*max_mismatch_reads_frac);
    max_mismatch_adapter[i] = floor(((float)i)*max_mismatch_adapter_frac);
    min_match_reads[i] = ceil(((float)i)*min_match_reads_frac);
    min_match_adapter[i] = ceil(((float)i)*min_match_adapter_frac);
  }
  //get length of forward and reverse primers
  int forward_primer_len = strlen(forward_primer);
  int reverse_primer_len = strlen(reverse_primer);


  gzFile ffq = fileOpen(forward_fn, "r");
  gzFile ffqw = fileOpen(forward_out_fn,"w");
  gzFile rfq = fileOpen(reverse_fn, "r");
  gzFile rfqw = fileOpen(reverse_out_fn,"w");
  gzFile mfqw = NULL;
  gzFile ppaw = NULL;
  if(do_read_merging)
    mfqw = fileOpen(merged_out_fn,"w");
  if(pretty_print)
    ppaw = fileOpen(pretty_print_fn,"w");




  while(next_fastqs( ffq, rfq, sqp, p64 )){ //returns false when done
    update_spinner(num_pairs++);


    //    //see if we can trim based on read-read overlap alone
    //    if(adapter_trim(sqp, min_ol_adapter,
    //        forward_primer, forward_primer_dummy_qual,
    //        forward_primer_len,
    //        reverse_primer, reverse_primer_dummy_qual,
    //        reverse_primer_len,
    //        min_match_adapter,
    //        max_mismatch_adapter,
    //        min_match_reads,
    //        max_mismatch_reads,
    //        qcut)){
    //      //we trimmed the adapter!
    //      num_adapter++;
    //      if((sqp->flen < min_read_len) || (sqp->rlen < min_read_len)){
    //        num_discarded++;
    //        continue;
    //      }
    //      if(!do_read_merging){ //just print
    //        write_fastq(ffqw, sqp->fid, sqp->fseq, sqp->fqual);
    //        write_fastq(rfqw, sqp->rid, sqp->rseq, sqp->rqual);
    //
    //      }else{ //force merge
    //        num_merged++;
    //        adapter_merge(sqp,false);
    //        write_fastq(mfqw,sqp->fid,sqp->merged_seq,sqp->merged_qual);
    //        if(pretty_print && num_pretty_print < max_pretty_print){
    //          num_pretty_print++;
    //          pretty_print_alignment(ppaw,sqp,qcut,true); //true b/c merged input sorted
    //        }
    //      }
    //      continue;
    //    }




    AlnAln *faaln, *raaln, *fraln;

    read_thresh = sqp->flen + sqp->rlen - (sqp->flen/2 * aln_param_rd2rd.gap_ext) - (sqp->rlen/2 * aln_param_rd2rd.gap_ext) - aln_param_rd2rd.gap_open*2 - aln_param_rd2rd.gap_end*2;

    faaln = aln_stdaln_aux(sqp->fseq, forward_primer, &aln_param_nt2nt,
        ALN_TYPE_LOCAL, adapter_thresh , sqp->flen, forward_primer_len);
    raaln = aln_stdaln_aux(sqp->rseq, reverse_primer, &aln_param_nt2nt,
        ALN_TYPE_LOCAL, adapter_thresh, sqp->rlen, reverse_primer_len);

    //check for direct adapter match.
    if(adapter_trim(sqp, min_ol_adapter,
        forward_primer, forward_primer_dummy_qual,
        forward_primer_len,
        reverse_primer, reverse_primer_dummy_qual,
        reverse_primer_len,
        min_match_adapter,
        max_mismatch_adapter,
        min_match_reads,
        max_mismatch_reads,
        qcut) ||
        faaln->score >= adapter_thresh ||
        raaln->score >= adapter_thresh){
      num_adapter++; //adapter present
      //print it if user wants
      if(pretty_print && num_pretty_print < max_pretty_print){
        //void pretty_print_alignment_stdaln(gzFile out, SQP sqp, AlnAln *aln, bool first_adapter, bool second_adapter)
        if(faaln->score >= adapter_thresh){
          num_pretty_print++;
          pretty_print_alignment_stdaln(ppaw,sqp,faaln,true,false,false);
        }
        if(raaln->score >= adapter_thresh){
          num_pretty_print++;
          pretty_print_alignment_stdaln(ppaw,sqp,raaln,false,true,false);
        }
      }

      //do stuff to it
      //assume full length adapter and squish it down to the read with no gaps
      int rpos,fpos;
      rpos = fpos = -MAX_SEQ_LEN;
      if(faaln->score >= adapter_thresh){
        fpos = faaln->start1 - faaln->start2;
      }
      if(raaln->score >= adapter_thresh){
        rpos = raaln->start1 - raaln->start2;
      }
      if(rpos == -MAX_SEQ_LEN && fpos == -MAX_SEQ_LEN){
        rpos = fpos = min(sqp->rlen,sqp->flen);
      }else if (rpos == -MAX_SEQ_LEN){
        rpos = max(0,min(fpos,sqp->rlen));
      }else if (fpos == -MAX_SEQ_LEN){
        fpos = max(0,min(rpos,sqp->flen));
      }else{
        rpos = max(0,min(rpos,min(fpos,sqp->rlen)));
        fpos = max(0,min(fpos,min(rpos,sqp->flen)));
      }
      sqp->rlen = rpos;
      sqp->flen = fpos;

      if(fpos < min_read_len || rpos < min_read_len){
        num_discarded++;
        goto CLEAN_ADAPTERS;
      }else{ //trim the adapters
        sqp->fseq[fpos] = '\0';
        sqp->fqual[fpos] = '\0';
        sqp->rc_rseq[rpos] = '\0';
        sqp->rc_rqual[rpos] = '\0';
        sqp->flen = fpos;
        sqp->rlen = rpos;
        strcpy(sqp->rseq,sqp->rc_rseq); //move RC reads into reg place and reverse them
        strcpy(sqp->rqual,sqp->rc_rqual);
        rev_qual(sqp->rqual,sqp->rlen);
        revcom_seq(sqp->rseq,sqp->rlen);
      }

      //do we want read merging?
      if(do_read_merging){

        //do a nice global alignment between two reads, and print consensus
        fraln = aln_stdaln_aux(sqp->fseq, sqp->rc_rseq, &aln_param_rd2rd,
            ALN_TYPE_GLOBAL, 1, sqp->flen, sqp->rlen);

        //there is some kind of alignment, so lets do the merge
        if(fraln->score > read_thresh){
          //write the merged sequence
          fill_merged_sequence(sqp, fraln, true);
          if(pretty_print && num_pretty_print < max_pretty_print){
            num_pretty_print += 1;
            pretty_print_alignment_stdaln(ppaw,sqp,fraln,false,false,true);
          }
          write_fastq(mfqw,sqp->fid,sqp->merged_seq,sqp->merged_qual);
        }
        goto CLEAN_ALL;
      }else{
        // we know that the adapters are present, trimmed, and the resulting
        // read lengths are both long enough to print.
        // We also know that we aren't doing merging.
        // Now we just need to print.
        write_fastq(ffqw, sqp->fid, sqp->fseq, sqp->fqual);
        write_fastq(rfqw, sqp->rid, sqp->rseq, sqp->rqual);


        //didn't do read alignment
        goto CLEAN_ADAPTERS;
      }
      goto CLEAN_ADAPTERS;
    }

    //check for strong read overlap to assist trimming ends of adapters from end of read
    if(do_read_merging){
      if(read_merge(sqp, min_ol_reads, min_match_reads, max_mismatch_reads, qcut)){
        //print merged output
        num_merged++;
        write_fastq(mfqw,sqp->fid,sqp->merged_seq,sqp->merged_qual);
        if(pretty_print && num_pretty_print < max_pretty_print){
          num_pretty_print++;
          pretty_print_alignment(ppaw,sqp,qcut,false); //false b/c merged input in fixed order
        }
      }else{
        //no significant overlap so just write them
        write_fastq(ffqw, sqp->fid, sqp->fseq, sqp->fqual);
        write_fastq(rfqw, sqp->rid, sqp->rseq, sqp->rqual);
      }
      //done
      goto CLEAN_ADAPTERS;
    }else{ //just write reads to output fastqs
      write_fastq(ffqw, sqp->fid, sqp->fseq, sqp->fqual);
      write_fastq(rfqw, sqp->rid, sqp->rseq, sqp->rqual);
      goto CLEAN_ADAPTERS;
    }


    /**
     * Section for heirarchial cleanup
     *
     * In every case we will at least have to free up the alignment between the adapter and two reads.
     * however in some cases there will be an additional alignment between the two reads. We can do
     * good cleanup in this case with gotos
     */

    CLEAN_ALL:
    aln_free_AlnAln(fraln);

    CLEAN_ADAPTERS:
    aln_free_AlnAln(faaln);
    aln_free_AlnAln(raaln);

    //End the loop over reads
  }
  end = clock();
  double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  fprintf(stderr,"\nPairs Processed:\t%lld\n",num_pairs);
  fprintf(stderr,"Pairs Merged:\t%lld\n",num_merged);
  fprintf(stderr,"Pairs With Adapters:\t%lld\n",num_adapter);
  fprintf(stderr,"Pairs Discarded:\t%lld\n",num_discarded);
  fprintf(stderr,"CPU Time Used (Minutes):\t%lf\n",cpu_time_used/60.0);



  SQP_destroy(sqp);
  gzclose(ffq);
  gzclose(ffqw);
  gzclose(rfq);
  gzclose(rfqw);
  if(mfqw != NULL)
    gzclose(mfqw);
  if(ppaw != NULL)
    gzclose(ppaw);
  return 0;
}
