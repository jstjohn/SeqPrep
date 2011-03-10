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

void SQP_destroy(SQP sqp){
  //free up an SQP
  free(sqp);
}

/* next_fastqs
   Read the next forward and reverse fastq sequences.
   Check to make sure their ID's are compatible and
   put the results in the next SQP of SQPDB. Grow
   this, if necessary.
*/
inline bool next_fastqs( gzFile ffq, gzFile rfq, SQP curr_sqp ) {
  int frs; // forward fastq read status
  int rrs; // reverse fastq read status
  int id1len = 0;
  int id2len = 0;
  /* Read the next fastq record from the forward and reverse
     pair of each */
  frs = read_fastq( ffq, curr_sqp->fid, curr_sqp->fseq, 
                    curr_sqp->fqual, &id1len, &(curr_sqp->flen) );
  rrs = read_fastq( rfq, curr_sqp->rid, curr_sqp->rseq, 
                    curr_sqp->rqual, &id2len, &(curr_sqp->rlen) );

  if ( (frs == 1) &&
       (rrs == 1) &&
       f_r_id_check( fid, id1len, rid, id2len ) ) {
    return true;
  }

  else {
    return false;
  }
}

inline bool f_r_id_check( char fid[], size_t fid_len, char rid[], size_t rid_len ) {
  if(fid_len != rid_len){
    return false; //trivial case
  }

  //expect last two characters are not equal
  if (strcmpi( fid, rid, fid_len - 2) == 0 ) {
    return true;
  }
  return false;
}

/* read_fastq
   Return 1 => more sequence to be had
          0 => EOF
 */
int read_fastq( gzFile* fastq, char id[], char seq[], char qual[], size_t *id_len, size_t *seq_len, int p64 ) {
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
     char subjectSeq[], char subjectQual[]; size_t subjectLen; 
     char querySeq[]; char queryQual[]; size_t queryLen; 
     size_t min_olap, size_t min_match, size_t max_miss, 
     bool check_unique, char adj_q_cut) {
  
  size_t  pos;  
  /* Try each possible starting position 
     on the forward sequence */
  int best_hit = CODE_NOMATCH;
  for( pos = 0; pos < subjectLen - min_olap; pos++ ) {
    if ( k_match( &(subjectSeq[pos]),
	       &(subjectQual[pos]),
	       subjectLen,
               querySeq, queryQual,
               queryLen, min_match 
               max_miss, adj_q_cut ) ) {
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
	size_t minmatch, size_t maxmiss, char adj_q_cut ) {
  size_t i;
  size_t mismatch = 0;
  size_t match = 0;
  for( i = 0; ((i < len1) && (i <len2)); i++ ) {
    if ( (q1[i] >= adj_q_cut) &&
	 (q2[i] >= adj_q_cut){
      if (s1[i] != s2[i]) ) {
        mismatch++;
        if(mismatch >=max_miss)
          return false;
      }else{
        match++;
      }
    }
  }
  if (match >= minmatch)
    return true;
  else
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
