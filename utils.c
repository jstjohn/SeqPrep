#include <ctype.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "stdaln.h"
#include "utils.h"


SQP SQP_init(){
  //allocate an SQP
  return (SQP) malloc(sizeof(Sqp));
}

/**
 * Calculates the resulting phred 33 score given a mismatch
 */
inline char mismatch_p33_merge(char pA, char pB){
  if(pA > pB){
    return pA-(pB-33);
  }else{
    return pB-(pA-33);
  }
}

bool isXDNA(char c){
  //A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  *  X
  switch(toupper(c)){
  case 'A':
  case 'C':
  case 'G':
  case 'T':
  case 'N':
  case 'X':
  case 'R':
  case 'W':
  case '.':
    return true;
  default:
    return false;
  }
}

/**
 * Calculates the resulting phred 33 score given a match
 */
inline char match_p33_merge(char pA, char pB){
  char res = pA+(pB-33);
  if(res > MAX_QUAL)
    return MAX_QUAL;
  return res;
}


void pretty_print_alignment_stdaln(gzFile out, SQP sqp, AlnAln *aln, bool first_adapter, bool second_adapter, bool print_merged){
  if(!(first_adapter || second_adapter)){
    gzprintf(out,"Read Alignment Score:%d, Suboptimal Score:%d\nID:%s\n",aln->score, aln->subo ,sqp->fid);
    gzprintf(out,"READ1: %s\n",aln->out1);
    gzprintf(out,"       %s\n",aln->outm);
    gzprintf(out,"READ2: %s\n",aln->out2);
    if(print_merged)
      gzprintf(out,"MERGD: %s\n\n",sqp->merged_seq);
    else
      gzprintf(out,"\n");
    return;
  }else if(first_adapter){
    gzprintf(out,"Adapter Alignment Score:%d, Suboptimal Score:%d\nID:%s\n",aln->score, aln->subo ,sqp->fid);
  }else if(second_adapter){
    gzprintf(out,"Adapter Alignment Score:%d, Suboptimal Score:%d\nID:%s\n",aln->score, aln->subo ,sqp->rid);
  }
  gzprintf(out,"READ: %s\n",aln->out1);
  gzprintf(out,"      %s\n",aln->outm);
  gzprintf(out,"ADPT: %s\n\n",aln->out2);
}


/**
 * Print the alignment in the pretty form:
 * Subj:        TGCTAGACTAGCATCG
 *              | |||*|||
 * Quer: ACGTGCATCCTANACT
 * Stars represent good quality mismatches,
 * spaces represent ignored mismatches,
 * vertical bars represent matches of any type
 *
 */
void pretty_print_alignment(gzFile out, SQP sqp, char adj_q_cut, bool sort){
  char *queryseq;
  char *queryqual;
  char *subjseq;
  char *subjqual;
  int querylen = 0;
  int subjlen = 0;
  int i;
  if((!sort) || (sqp->flen >= sqp->rlen)){
    subjseq = sqp->fseq;
    subjqual = sqp->fqual;
    queryseq = sqp->rc_rseq;
    queryqual = sqp->rc_rqual;
    querylen = sqp->rlen;
    subjlen = sqp->flen;
  }else{
    subjseq = sqp->rc_rseq;
    subjqual = sqp->rc_rqual;
    queryseq = sqp->fseq;
    queryqual = sqp->fqual;
    querylen = sqp->flen;
    subjlen = sqp->rlen;
  }
  gzprintf(out, "ID: %s\n",sqp->fid);
  gzprintf(out, "SUBJ: %s\n",subjseq);
  //now print out the bars
  gzprintf(out, "      "); //initial space
  for(i=0;i<sqp->merged_len;i++){
    if(i >= sqp->mpos && i < subjlen && i < (querylen + sqp->mpos)){
      //we are in the overlapping region
      if(subjseq[i] == queryseq[i-sqp->mpos])
        gzputc(out,'|');
      else if(subjqual[i] < adj_q_cut || queryqual[i-sqp->mpos] < adj_q_cut)
        gzputc(out,' ');
      else
        gzputc(out,'*');
    }else{
      gzputc(out,' ');
    }
  }
  gzprintf(out,"\nQUER: ");
  for(i=0;i<sqp->mpos;i++)
    gzputc(out,' '); //spaces before aln
  gzprintf(out,"%s",queryseq);
  gzprintf(out,"\nMERG: %s\n\n",sqp->merged_seq);
}

void fill_merged_sequence(SQP sqp, AlnAln *aln, bool trim_overhang){
  int len = strlen(aln->out1);
  char *out1, *out2;
  out1 = aln->out1;
  out2 = aln->out2;
  int i,p1,p2; //p1,2 store pointers to corresponding pos in original seqs
  p1 = p2 = 0;
  char c1,c2,q1,q2,t1,t2;
  bool end_gaps;
  bool begin_gaps = trim_overhang;
  int j = 0;
  int k;
  for(i=0;i<len;i++){
    c1 = toupper(out1[i]);
    c2 = toupper(out2[i]);
    q1 = sqp->fqual[p1];
    q2 = sqp->rc_rqual[p2];
    if(isXDNA(c1) && isXDNA(c2)){

      //case 1 both are DNA, choose one with best score and subtract
      if (begin_gaps) begin_gaps = false; //switch it off now that we have seen a match
      if(c1 == c2){
        sqp->merged_seq[j] = c1;
        sqp->merged_qual[j] = q1+q2-33;
      }else if(q2 > q1){
        sqp->merged_seq[j] = c2;
        sqp->merged_qual[j] = q2 - q2 + 33;
      }else{
        sqp->merged_seq[j] = c1;
        sqp->merged_qual[j] = q1 - q2 + 33;
      }
      //increment both positions of the reads
      p1++;
      p2++;
      j++;
    }else if(isXDNA(c1)){
      // c2 is a gap
      if (!begin_gaps){
        sqp->merged_seq[j] = c1;
        sqp->merged_qual[j] = ((q1-33)>>1)+33; //divide score by 2
        //now check to see if we are done:
        if(trim_overhang){
          end_gaps = true;
          for(k=i;k<len;k++){
            t2 = out2[k];
            if(t2 != '-'){
              end_gaps = false;
              break;
            }
          }
          if(end_gaps){
            //everything after this is a gap
            break;
          }
        }
        j++;
      }
      //increment the first
      p1++;
    }else if(isXDNA(c2)){
      //c1 is a gap
      if(!begin_gaps){
        sqp->merged_seq[j] = c2;
        sqp->merged_qual[j] = ((q2-33)>>1)+33; //divide score by 2
        if(trim_overhang){
          end_gaps = true;
          for(k=i;k<len;k++){
            t1 = out1[k];
            if(t1 != '-'){
              end_gaps = false;
              break;
            }
          }
          if(end_gaps){
            //everything after this is a gap
            break;
          }
        }
        j++;
      }
      //increment the second
      p2++;
    }
  }
  sqp->merged_seq[j] = '\0';
  sqp->merged_qual[j] = '\0';
  sqp->merged_len = j;
}



void SQP_destroy(SQP sqp){
  //free up an SQP
  free(sqp);
}


/**
 * adapter_trim:
 *
 *
 */
bool adapter_trim(SQP sqp, size_t min_ol_adapter,
    char *forward_primer, char *forward_primer_dummy_qual,
    int forward_primer_len,
    char *reverse_primer, char *reverse_primer_dummy_qual,
    int reverse_primer_len,
    unsigned short min_match_adapter[MAX_SEQ_LEN+1],
    unsigned short max_mismatch_adapter[MAX_SEQ_LEN+1],
    unsigned short min_match_reads[MAX_SEQ_LEN+1],
    unsigned short max_mismatch_reads[MAX_SEQ_LEN+1],
    char qcut){
  //adapters on reads if the insert size is less than the read length, the adapter
  // appears at the end of the sequence.


  /**
   * First check for adapter match before the first position of the read
   */
  int pfpos = compute_ol(
      forward_primer, forward_primer_dummy_qual, forward_primer_len,
      sqp->fseq,sqp->fqual,sqp->flen,
      max(min(forward_primer_len,sqp->flen)-5,0), min_match_adapter, max_mismatch_adapter,
      false, qcut);

  int prpos = compute_ol(
      reverse_primer, reverse_primer_dummy_qual, reverse_primer_len,
      sqp->rseq,sqp->rqual,sqp->rlen,
      max(min(reverse_primer_len,sqp->rlen)-5,0), min_match_adapter, max_mismatch_adapter,
      false, qcut);

  if(pfpos >= 0 || prpos >= 0){
    //yikes, a match to the adapter at the first position!
    sqp->fseq[0] = '\0';
    sqp->fqual[0] = '\0';
    sqp->flen = 0;
    sqp->rseq[0] = '\0';
    sqp->rqual[0] = '\0';
    sqp->rlen = 0;
    sqp->rc_rqual[0] = '\0';
    sqp->rc_rseq[0] = '\0';
    return true;
  }

  /**
   * now check for the adapter after the first position of the read
   */
  int fpos = compute_ol(sqp->fseq,sqp->fqual,sqp->flen,
      forward_primer, forward_primer_dummy_qual, forward_primer_len,
      min_ol_adapter, min_match_adapter, max_mismatch_adapter,
      false, qcut);
  int rpos = compute_ol(sqp->rseq,sqp->rqual,sqp->rlen,
      reverse_primer, reverse_primer_dummy_qual, reverse_primer_len,
      min_ol_adapter, min_match_adapter, max_mismatch_adapter,
      false, qcut);
  if(fpos != CODE_NOMATCH || rpos != CODE_NOMATCH){
    //check if reads are long enough to do anything with.
    // trim adapters
    sqp->fseq[fpos] = '\0';
    sqp->fqual[fpos] = '\0';
    sqp->flen = fpos;
    sqp->rseq[rpos] = '\0';
    sqp->rqual[rpos] = '\0';
    sqp->rlen = rpos;
    // now re-reverse complement the sequences
    strncpy(sqp->rc_rseq,sqp->rseq,sqp->rlen+1);
    strncpy(sqp->rc_rqual,sqp->rqual,sqp->rlen+1);
    rev_qual(sqp->rc_rqual, sqp->rlen);
    revcom_seq(sqp->rc_rseq, sqp->rlen);
    //adapters present
    return true;
  }

  return read_olap_adapter_trim(sqp, min_ol_adapter,
      min_match_adapter, max_mismatch_adapter,
      min_match_reads, max_mismatch_reads,
      qcut);
}



/**
 * look for adapters by read overlap
 *
 */
bool read_olap_adapter_trim(SQP sqp, size_t min_ol_adapter,
    unsigned short min_match_adapter[MAX_SEQ_LEN+1],
    unsigned short max_mismatch_adapter[MAX_SEQ_LEN+1],
    unsigned short min_match_reads[MAX_SEQ_LEN+1],
    unsigned short max_mismatch_reads[MAX_SEQ_LEN+1],
    char qcut){
  ////////////
  // Look at the adapter overhang
  // Starting from our minimum adapter overlap
  // check to see if there is total overlap with
  //Round1:
  //       ---------- Subj
  //       ---------- Query
  //Round2:
  //      ----------  Subj
  //     ----------   Query
  //...
  //we can get this effect by swapping the query and subj, and then have a high minimum
  //overlap
  char *queryseq= sqp->rc_rseq;
  char *queryqual= sqp->rc_rqual;
  char *subjseq= sqp->fseq;
  char *subjqual= sqp->fqual;
  int querylen = sqp->rlen;
  int subjlen = sqp->flen;

  int ppos = compute_ol(
      queryseq, queryqual, querylen,
      subjseq, subjqual, subjlen,
      min(querylen,subjlen)-min_ol_adapter, min_match_reads, max_mismatch_reads,
      true, qcut ); //pass true here so ambiguous matches are avoided
  if(ppos != CODE_NOMATCH && ppos != CODE_AMBIGUOUS){
    //we have a match, trim the adapter!
    if(ppos == 0){
      //no adapter
      return false;
    }else{
      //ppos gives us the shift to the left of the query
      // One case:
      //   ----X------- fread
      // -X----         rread
      // Another case:
      //   ---X-        fread
      // -X---          rread
      // Another case:
      //   ----         fread
      // -X----X-       rread


      //first calc rlen after the first clip
      sqp->rlen -= ppos;

      //now in the first two cases shown above, the other cut point is just the
      //new rlen
      if(sqp->rlen <= sqp->flen)
        sqp->flen = sqp->rlen;
      //otherwise leave sqp->flen alone
      else if(sqp->rlen > sqp->flen){
        // Another case:
        //   ----         fread
        // -X----X---     rread
        // make initial cut to rc read
        sqp->rc_rqual[ppos + sqp->flen] = '\0';
        sqp->rc_rseq[ppos + sqp->flen] = '\0';
        strncpy(sqp->rseq,sqp->rc_rseq,ppos + sqp->flen+1); //move RC reads into reg place and reverse them
        strncpy(sqp->rqual,sqp->rc_rqual,ppos + sqp->flen+1);
        rev_qual(sqp->rqual, ppos + sqp->flen);
        revcom_seq(sqp->rseq, ppos + sqp->flen);

        //now we have our end cut in place in the regular reads
        sqp->rlen = sqp->flen;

      }

      //now cases have been handled and length has been determined
      sqp->fseq[sqp->flen] = '\0';
      sqp->fqual[sqp->flen] = '\0';
      sqp->rseq[sqp->rlen] = '\0';
      sqp->rqual[sqp->rlen] = '\0';
      // now re-reverse complement the sequences
      strncpy(sqp->rc_rseq,sqp->rseq,sqp->rlen+1);
      strncpy(sqp->rc_rqual,sqp->rqual,sqp->rlen+1);
      rev_qual(sqp->rc_rqual, sqp->rlen);
      revcom_seq(sqp->rc_rseq, sqp->rlen);
      return true;
    }
  }
  return false;
}


/**
 * read_merge:
 *    Computes the potential overlap between two reads,
 *    fills the merged_seq items in sqp
 *    return true if a merging was done, and false otherwise
 */
bool read_merge(SQP sqp, size_t min_olap,
    unsigned short min_match[MAX_SEQ_LEN+1],
    unsigned short max_mismatch[MAX_SEQ_LEN+1],
    char adj_q_cut){
  //now compute overlap
  int i;

  char *queryseq;
  char *queryqual;
  char *subjseq;
  char *subjqual;
  int querylen = 0;
  int subjlen = 0;
  char c,q;
  //  if(sqp->rlen <= sqp->flen){
  subjseq = sqp->fseq;
  subjqual = sqp->fqual;
  queryseq = sqp->rc_rseq;
  queryqual = sqp->rc_rqual;
  querylen = sqp->rlen;
  subjlen = sqp->flen;
  //  }else{
  //    queryseq = sqp->fseq;
  //    queryqual = sqp->fqual;
  //    subjseq = sqp->rc_rseq;
  //    subjqual = sqp->rc_rqual;
  //    subjlen = sqp->rlen;
  //    querylen = sqp->flen;
  //  }
  ////////////
  // Now calculate the other cases
  //Round1:
  //   ---------- Subj
  //   ---------- Query
  //Round2:
  //  ----------  Subj
  //   ---------- Query
  //Round3:
  // ----------   Subj
  //   ---------- Query
  //...
  int mpos = compute_ol(
      subjseq, subjqual, subjlen,
      queryseq, queryqual, querylen,
      min_olap, min_match, max_mismatch,
      true, adj_q_cut );
  if(mpos == CODE_NOMATCH || mpos == CODE_AMBIGUOUS){
    return false;
  }else{
    //part where subj is non-overlapping
    // ---------
    //  -----
    int pos = 0;
    for(i=0;i<mpos;i++){
      sqp->merged_seq[pos] = subjseq[i];
      sqp->merged_qual[pos] = subjqual[i];
      pos++;
    }
    //overlapping section
    int end = min(subjlen,querylen+mpos);
    for(i=mpos;i<end ;i++){
      if(subjseq[i] == queryseq[i-mpos]){
        c = subjseq[i];
        q = match_p33_merge(subjqual[i],queryqual[i-mpos]);
      }else{
        q = mismatch_p33_merge(subjqual[i],queryqual[i-mpos]);
        if(subjqual[i] > queryqual[i-mpos]){
          c = subjseq[i];
        }else{
          c = queryseq[i-mpos];
        }
      }
      sqp->merged_seq[pos] = c;
      sqp->merged_qual[pos] = q;
      pos++;
    }
    //now print either the subject or query that is remaining at the end.
    if(subjlen >= querylen+mpos){ //subject is really long so it is overhanging
      for(i=end;i<subjlen;i++){
        sqp->merged_seq[pos] = subjseq[i];
        sqp->merged_qual[pos] = subjqual[i];
        pos++;
      }
    }else{ //normal case
      for(i=end-mpos;i<querylen;i++){
        sqp->merged_seq[pos] = queryseq[i];
        sqp->merged_qual[pos] = queryqual[i];
        pos++;
      }
    }
    sqp->merged_len = pos;
    sqp->merged_seq[pos] = '\0';
    sqp->merged_qual[pos] = '\0';
    sqp->mpos = mpos;
    return true; //successfull merge complete!
  }
  return false;
}


void adapter_merge(SQP sqp, bool print_overhang){
  //first RC reverse read so we can do direct overlapping
  int i = 0;
  int j = 0;
  char c,q;
  if(sqp->rlen == sqp->flen){
    //easy.. peezy.. lemaon.... squeezy..
    for(i=0; i< sqp->rlen; i++){
      if(sqp->rc_rseq[i] == sqp->fseq[i]){
        c = sqp->rc_rseq[i];
        q = match_p33_merge(sqp->rc_rqual[i],sqp->fqual[i]);
      }else{
        q = mismatch_p33_merge(sqp->rc_rqual[i],sqp->fqual[i]);
        if(sqp->rc_rqual[i]>sqp->fqual[i]){
          c = sqp->rc_rseq[i];
        }else{
          c = sqp->fseq[i];
        }
      }
      sqp->merged_seq[i] = c;
      sqp->merged_qual[i] = q;
    }
    sqp->merged_len = sqp->rlen;
    sqp->merged_qual[i] = '\0';
    sqp->merged_seq[i] = '\0';
    sqp->mpos = 0;
  }else{
    int num_match = 0;
    int max_match = 0;
    int max_offset = 0;
    char *queryseq;
    char *queryqual;
    char *subjseq;
    char *subjqual;
    int querylen = 0;
    int subjlen = 0;
    int ndiff = 0;
    if(sqp->flen > sqp->rlen){
      ndiff = sqp->flen - sqp->rlen;
      subjseq = sqp->fseq;
      subjqual = sqp->fqual;
      queryseq = sqp->rc_rseq;
      queryqual = sqp->rc_rqual;
      querylen = sqp->rlen;
      subjlen = sqp->flen;
    }else{
      ndiff = sqp->rlen - sqp->flen;
      subjseq = sqp->rc_rseq;
      subjqual = sqp->rc_rqual;
      queryseq = sqp->fseq;
      queryqual = sqp->fqual;
      querylen = sqp->flen;
      subjlen = sqp->rlen;
    }
    for(i=0;i<=ndiff;i++){
      num_match = 0;
      for(j=0;j<querylen;j++){
        if(subjseq[i+j] == queryseq[j])
          num_match++;
      }
      if(num_match>max_match){
        max_offset = i;
        max_match = num_match;
      }
    }//end for
    //now we have our best offset, use that and build our merged_sequence
    int pos = 0;

    if(print_overhang){
      for(i=0;i<max_offset;i++){
        sqp->merged_seq[pos] = subjseq[i];
        sqp->merged_qual[pos] = subjqual[i];
        pos++;
      }
    }
    //now do the part that overlaps
    for(i=max_offset;i<querylen+max_offset;i++){
      if(subjseq[i] == queryseq[i-max_offset]){
        c = subjseq[i];
        q = match_p33_merge(subjqual[i],queryqual[i-max_offset]);
      }else{
        q = mismatch_p33_merge(subjqual[i],queryqual[i-max_offset]);
        if(subjqual[i] > queryqual[i-max_offset]){
          c = subjseq[i];
        }else{
          c = queryseq[i-max_offset];
        }
      }
      sqp->merged_seq[pos] = c;
      sqp->merged_qual[pos] = q;
      pos++;
    }
    //finish off the subject sequence
    if(print_overhang){
      for(i=(querylen+max_offset); i<subjlen; i++){
        sqp->merged_seq[pos] = subjseq[i];
        sqp->merged_qual[pos] = subjqual[i];
        pos++;
      }
    }
    sqp->merged_seq[pos]='\0';
    sqp->merged_qual[pos]='\0';
    sqp->merged_len = pos;
    sqp->mpos = max_offset;
  }//end case where we need to find the best merging
}

/* next_fastqs
   Read the next forward and reverse fastq sequences.
   Check to make sure their ID's are compatible and
   put the results in the next SQP of SQPDB. Grow
   this, if necessary.
 */
inline bool next_fastqs( gzFile ffq, gzFile rfq, SQP curr_sqp, bool p64 ) {
  int frs; // forward fastq read status
  int rrs; // reverse fastq read status
  size_t id1len = 0;
  size_t id2len = 0;
  /* Read the next fastq record from the forward and reverse
     pair of each */

  //make sure everything is fresh...
  memset(curr_sqp->fid,'\0',MAX_SEQ_LEN);
  memset(curr_sqp->rid,'\0',MAX_SEQ_LEN);
  memset(curr_sqp->fseq,'\0',MAX_SEQ_LEN);
  memset(curr_sqp->rseq,'\0',MAX_SEQ_LEN);
  memset(curr_sqp->rc_rseq,'\0',MAX_SEQ_LEN);
  memset(curr_sqp->fqual,'\0',MAX_SEQ_LEN);
  memset(curr_sqp->rqual,'\0',MAX_SEQ_LEN);
  memset(curr_sqp->merged_seq,'\0',MAX_SEQ_LEN+MAX_SEQ_LEN);
  memset(curr_sqp->merged_qual,'\0',MAX_SEQ_LEN+MAX_SEQ_LEN);
  memset(curr_sqp->rc_rqual,'\0',MAX_SEQ_LEN);
  curr_sqp->flen = curr_sqp->rlen = 0;


  //

  frs = read_fastq( ffq, curr_sqp->fid, curr_sqp->fseq, 
      curr_sqp->fqual, &id1len, &(curr_sqp->flen), p64 );
  rrs = read_fastq( rfq, curr_sqp->rid, curr_sqp->rseq, 
      curr_sqp->rqual, &id2len, &(curr_sqp->rlen), p64 );

  //  //reverse comp the second read for overlapping and everything.
  //  strcpy(curr_sqp->rc_rseq,curr_sqp->rseq);
  //  strcpy(curr_sqp->rc_rqual,curr_sqp->rqual);
  //  revcom_seq(curr_sqp->rc_rseq,curr_sqp->rlen);
  //  rev_qual(curr_sqp->rc_rqual,curr_sqp->rlen);

  if ( (frs == 1) &&
      (rrs == 1) &&
      f_r_id_check( curr_sqp->fid, id1len, curr_sqp->rid, id2len ) ) {
    strncpy(curr_sqp->rc_rseq,curr_sqp->rseq,curr_sqp->rlen+1);
    strncpy(curr_sqp->rc_rqual,curr_sqp->rqual,curr_sqp->rlen+1);
    rev_qual(curr_sqp->rc_rqual, curr_sqp->rlen);
    revcom_seq(curr_sqp->rc_rseq, curr_sqp->rlen);
    return true;
  } else {
    return false;
  }
}

inline int write_fastq(gzFile out, char id[], char seq[], char qual[]){
  return gzprintf(out,"@%s\n%s\n+\n%s\n", id, seq, qual);
}


inline bool f_r_id_check( char fid[], size_t fid_len, char rid[], size_t rid_len ) {
  if(fid_len != rid_len){
    goto bad_read;
  }else if (strncmp( fid, rid, fid_len - 2) == 0 ) {
    return true;
  }

  bad_read:
  fprintf(stderr,"ERROR: Fastq id lines do not match: %s vs %s \n", fid, rid);
  return false;
}

/* read_fastq
   Return 1 => more sequence to be had
          0 => EOF
 */
int read_fastq( gzFile* fastq, char id[], char seq[], char qual[], size_t *id_len, size_t *seq_len, bool p64 ) {
  char c;
  size_t i;
  c = gzgetc( fastq );
  if ( c == EOF ) return 0;
  if ( c != '@' ) {
    fprintf( stderr, "fastq record not beginning with @\n" );
    return 0;
  }

  /* get identifier */
  i = 0;
  while( (!isspace(c=gzgetc( fastq ) ) &&
      (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return 0;
    }
    id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      /* Id is too long - truncate it now */
      id[i] = '\0';
    }
  }
  id[i] = '\0';
  *id_len = i;
  /* Now, everything else on the line is description (if anything)
     although fastq does not appear to formally support description */
  while ( (c != '\n') &&
      (c != EOF) ) {
    c = gzgetc( fastq );
  }

  /* Now, read the sequence. This should all be on a single line */
  i = 0;
  c = gzgetc( fastq );
  while ( (c != '\n') &&
      (c != EOF) &&
      (i < MAX_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      seq[i++] = (c=='.' ? 'N':c);
    }
    c = gzgetc( fastq );
  }
  seq[i] = '\0';
  *seq_len = i;
  /* If the reading stopped because the sequence was longer than
     INIT_ALN_SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_SEQ_LEN ) {
    while ( (c != '\n') &&
        (c != EOF) ) {
      c = gzgetc( fastq );
    }
  }

  /* Now, read the quality score header */
  c = gzgetc( fastq );
  if ( c != '+' ) {
    fprintf( stderr, "Problem reading quality line for %s\n", id );
    return 1;
  }
  /* Zip through the rest of the line, it should be the same identifier
     as before or blank */
  c = gzgetc( fastq );
  while( (c != '\n') &&
      (c != EOF) ) {
    c = gzgetc( fastq );
  }

  /* Now, get the quality score line */
  c = gzgetc( fastq );
  i = 0;
  while( (c != '\n') &&
      (c != EOF) &&
      (i < MAX_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      qual[i++] = (p64)?((c=='B')?'!':c-31):c;
    }
    c = gzgetc( fastq );
  }
  qual[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     INIT_ALN_SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_SEQ_LEN ) {
    while ( (c != '\n') &&
        (c != EOF) ) {
      c = gzgetc( fastq );
    }
  }

  if ( c == EOF ) {
    fprintf(stderr,"\nWarning: Your last read may have been discarded because you are missing a new line at the end of the file.\n\n");
    return 0;
  }
  return 1;
}


/** fileOpen **/
gzFile * fileOpen(const char *name, char access_mode[]) {
  gzFile * f;
  f = gzopen(name, access_mode);
  if (f == NULL) {
    fprintf( stderr, "%s\n", name);
    perror("Cannot open file");
    return NULL;
  }
  return f;
}






/*
   Supply two sequences in the proper orientation for overlap
   Ie in this example give compute_ol the reversed sequence and quality


   then QUERY:       TGCTAGACTAGCATCG
                     |*||| |||
     SUBJECT: ACGTGCATCCTANACT

   Therefore, the overlap would be 9. Ignore any
   base that has Quality score less than adj_q_cut

 */

int compute_ol(
    char subjectSeq[], char subjectQual[], size_t subjectLen,
    char querySeq[], char queryQual[], size_t queryLen,
    size_t min_olap, unsigned short min_match[MAX_SEQ_LEN+1],
    unsigned short max_mismatch[MAX_SEQ_LEN+1],
    bool check_unique, char adj_q_cut ) {

  size_t  pos;  
  /* Try each possible starting position 
     on the forward sequence */
  int best_hit = CODE_NOMATCH;
  int subject_len = subjectLen;
  for( pos = 0; pos < subjectLen - min_olap + 4; pos++ ) {
    subject_len = subjectLen - pos;
    //Round1:
    //   ------     Subj
    //   ---------- Query
    //Round2:
    //  ------      Subj
    //   ---------- Query
    //Round3:
    // ------       Subj
    //   ---------- Query
    //...
    if ( k_match( &(subjectSeq[pos]), &(subjectQual[pos]),
        subject_len, querySeq, queryQual, queryLen,
        min_match[subject_len>queryLen?queryLen:subject_len],
        max_mismatch[subject_len>queryLen?queryLen:subject_len],
        adj_q_cut ) ) {

      if(check_unique && best_hit != CODE_NOMATCH){
        return CODE_AMBIGUOUS;
      }
      if(best_hit == CODE_NOMATCH){
        if(!check_unique)
          return pos;
        best_hit = pos;
      }
    }
  }
  return best_hit;
}

/* k_match
   Args: pointer to forward seq,
         pointer to forward qual scores,
	 length of the forward seq
	 pointer to rev seq,
	 pointer to rev qual scores
         length of the reverse seq
   This is the comparison function for finding the overlap
   between the forward and reverse reads. It's called at 
   all possible overlapping positions, from longest to
   shortest, until it finds one. It doesn't require a match
   if either read has quality score less that QCUT.
   Returns: true if it's a match, false if it's not
 */
bool k_match( const char* s1, const char* q1, size_t len1, 
    const char* s2, const char* q2, size_t len2,
    unsigned short min_match, unsigned short max_mismatch,
    char adj_q_cut) {
  size_t i;
  size_t mismatch = 0;
  size_t match = 0;
  for( i = 0; ((i < len1) && (i <len2)); i++ ) {
    //if we have a match, or at least bad quality bases...
    if ( s1[i] == s2[i] || ((q1[i] >= adj_q_cut) &&
        (q2[i] >= adj_q_cut))){
      if (s1[i] != s2[i]) {
        mismatch++;
        if(mismatch > max_mismatch)
          return false;
      }else{
        match++;
      }
    }
  }
  if (match >= min_match)
    return true;
  return false;
}


void revcom_seq( char seq[], int len) {
  //int len = strlen(seq);
  char tmp_base;
  int  i;

  for (i = 0; i < len/2; i++) {
    tmp_base = seq[i];
    seq[i] = revcom_char(seq[len-(i+1)]);
    seq[len-(i+1)] = revcom_char(tmp_base);
  }

  /* If sequence length is even, we're done, otherwise there is
     the base right in the center to revcom */
  if (len%2 == 1) {
    seq[i] = revcom_char(seq[len-(i+1)]);
  }
}

inline char revcom_char(const char base) {
  static bool warned = false;
  switch (base) {
  case 'A':
    return 'T';
  case 'a' :
    return 't';

  case 'C':
    return 'G';
  case 'c' :
    return 'g';

  case 'G':
    return 'C';
  case 'g' :
    return 'c';

  case 'T':
    return 'A';
  case 't' :
    return 'a';

    //ignore special characters
  case '-':
    return '-';
  case '.':
    return '.';
  case 'N':
    return 'N';
  case 'n':
    return 'n';

  case 'X':
    return 'X';
  case 'x':
    return 'x';

  default:
    if(!warned){
      warned = true;
      fprintf( stderr, "WARNING: Non standard DNA character in sequence: \"%c\"\n", base);
    }
    return base;
  }
}

inline void rev_qual( char q[], int len ) {
  //int len = strlen(q);
  char tmp_q;
  int  i;

  for (i = 0; i < len/2; i++) {
    tmp_q = q[i];
    q[i] = q[len-(i+1)];
    q[len-(i+1)] = tmp_q;
  }
}
