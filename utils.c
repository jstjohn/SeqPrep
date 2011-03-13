#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
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

/**
 * Calculates the resulting phred 33 score given a match
 */
inline char match_p33_merge(char pA, char pB){
  char res = pA+(pB-33);
  if(res > MAX_QUAL)
    return MAX_QUAL;
  return res;
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
void pretty_print_alignment(gzFile out, SQP sqp, char adj_q_cut){
  char *queryseq;
  char *queryqual;
  char *subjseq;
  char *subjqual;
  int querylen = 0;
  int subjlen = 0;
  int i;
  if(sqp->flen >= sqp->rlen){
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


void SQP_destroy(SQP sqp){
  //free up an SQP
  free(sqp);
}

/**
 * read_merge:
 *    Computes the potential overlap between two reads,
 *    fills the merged_seq items in sqp
 *    return true if a merging was done, and false otherwise
 *
 *  WARNING: this function currently works by assuming that
 *  the second read is not a subset of the first read. This can be
 *  guarenteed if the second read is the same length as the first.
 */
bool read_merge(SQP sqp, size_t min_olap,
    unsigned short min_match[MAX_SEQ_LEN+1],
    unsigned short max_mismatch[MAX_SEQ_LEN+1],
    char adj_q_cut){
  //first reverse complement reads
  strcpy(sqp->rc_rseq,sqp->rseq);
  strcpy(sqp->rc_rqual,sqp->rqual);
  rev_qual(sqp->rc_rqual,sqp->rlen);
  revcom_seq(sqp->rc_rseq,sqp->rlen);
  //now compute overlap
  int i;

  char *queryseq;
  char *queryqual;
  char *subjseq;
  char *subjqual;
  int querylen = 0;
  int subjlen = 0;
  char c,q;
  subjseq = sqp->fseq;
  subjqual = sqp->fqual;
  queryseq = sqp->rc_rseq;
  queryqual = sqp->rc_rqual;
  querylen = sqp->rlen;
  subjlen = sqp->flen;

  int mpos = compute_ol(
      subjseq, subjqual, subjlen,
      queryseq, queryqual, querylen,
      min_olap, min_match, max_mismatch,
      true, adj_q_cut );
  if(mpos == CODE_NOMATCH || mpos == CODE_AMBIGUOUS){
    return false;
  }else{
    //part where subj is non-overlapping
    int pos = 0;
    for(i=0;i<mpos;i++){
      sqp->merged_seq[pos] = subjseq[i];
      sqp->merged_qual[pos] = subjqual[i];
      pos++;
    }
    //overlapping section
    for(i=mpos;i<subjlen;i++){
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
    //part where query is non-overlapping
    for(i=subjlen-mpos;i<querylen;i++){
      sqp->merged_seq[pos] = queryseq[i];
      sqp->merged_qual[pos] = queryqual[i];
      pos++;
    }
    sqp->merged_len = pos;
    sqp->merged_seq[pos] = '\0';
    sqp->merged_qual[pos] = '\0';
    sqp->mpos = mpos;
    return true; //successfull merge complete!
  }
  return false;
}


void adapter_merge(SQP sqp){
  //first RC reverse read so we can do direct overlapping
  strcpy(sqp->rc_rseq, sqp->rseq);
  strcpy(sqp->rc_rqual, sqp->rqual);
  revcom_seq(sqp->rc_rseq,sqp->rlen);
  rev_qual(sqp->rc_rqual,sqp->rlen);
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
    for(i=0;i<max_offset;i++){
      sqp->merged_seq[pos] = subjseq[i];
      sqp->merged_qual[pos] = subjqual[i];
      pos++;
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
    for(i=(querylen+max_offset); i<subjlen; i++){
      sqp->merged_seq[pos] = subjseq[i];
      sqp->merged_qual[pos] = subjqual[i];
      pos++;
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
    return true;
  }

  else {
    return false;
  }
}

inline int write_fastq(gzFile out, char id[], char seq[], char qual[]){
  return gzprintf(out,"@%s\n%s\n+\n%s\n", id, seq, qual);
}


inline bool f_r_id_check( char fid[], size_t fid_len, char rid[], size_t rid_len ) {
  if(fid_len != rid_len){
    return false; //trivial case
  }

  //expect last two characters are not equal
  if (strncmp( fid, rid, fid_len - 2) == 0 ) {
    return true;
  }
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
      seq[i++] = c;
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
                     | |||-|||
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
  for( pos = 0; pos < subjectLen - min_olap; pos++ ) {
    subject_len = subjectLen - pos;
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
        if(mismatch >=max_mismatch)
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


void revcom_seq( char seq[], int len ) {
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

  case '-':
  return '-';

  case 'N':
  return 'N';
  case 'n':
  return 'n';

  case 'X':
  return 'X';
  case 'x':
  return 'x';

  default:
  fprintf( stderr, "Do not know how to revcom \"%c\"\n", base);
  return base;
  }
}

inline void rev_qual( char q[], int len ) {
  char tmp_q;
  int  i;

  for (i = 0; i < len/2; i++) {
    tmp_q = q[i];
    q[i] = q[len-(i+1)];
    q[len-(i+1)] = tmp_q;
  }
}
