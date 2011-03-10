#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "seqp.h"

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
inline int next_fastqs( gzFile ffq, gzFile rfq, SQP curr_sqp ) {
  int frs; // forward fastq read status
  int rrs; // reverse fastq read status
 
  /* Read the next fastq record from the forward and reverse
     pair of each */
  frs = read_fastq( ffq, curr_sqp->fid, curr_sqp->fseq, curr_sqp->fqual );
  rrs = read_fastq( rfq, curr_sqp->rid, curr_sqp->rseq, curr_sqp->rqual );

  if ( (frs == 1) &&
       (rrs == 1) &&
       f_r_id_check( fid, rid ) ) {
    return 1;
  }

  else {
    return 0;
  }
}

inline int f_r_id_check( char fid[], size_t fid_len, char rid[], size_t rid_len ) {
  if(fid_len != rid_len){
    return 0; //trivial case
  }

  //expect last two characters are not equal
  if (strcmpi( fid, rid, fid_len - 2) == 0 ) {
    return 1;
  }
  return 0;
}

/* read_fastq
   Return 1 => more sequence to be had
          0 => EOF
 */
int read_fastq( gzFile* fastq, char id[], char seq[], char qual[], size_t *id_len, size_t *seq_len ) {
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
	  (i < SEQ_LEN) ) {
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
  if ( i == SEQ_LEN ) {
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
	 (i < SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      qual[i++] = c;
    }
    c = gzgetc( fastq );
  }
  qual[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     INIT_ALN_SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == SEQ_LEN ) {
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

/* Takes two SQP, a cutoff for base quality score,
   a shortcut number of mismatches beyond which we don't care
   anymore and pointers to the number of matches and mismatches to
   be calculated by comparing these two sequences.
   The forward and reverse fragments are compared. Just considers
   out to the end of the shorter of any two sequences being
   compared.
*/
void comp_seq( const SQP s1, const SQP s2, const int qcut, 
	       const int mmc, int* match, int* mismatch ) {
  size_t len1, len2, i;
  *match = 0;
  *mismatch = 0;

  len1 = strlen( s1->fseq );
  len2 = strlen( s2->fseq );

  for ( i = 0; ((i < len1) && (i < len2)); i++ ) {
    if ( (s1->fqual[i] >= qcut) &&
	 (s2->fqual[i] >= qcut) ) {
      if ( s1->fseq[i] == s2->fseq[i] ) {
	(*match)++;
      }
      else {
	(*mismatch)++;
	if ( *mismatch > mmc ) {
	  return;
	}
      }
    }
  }

  len1 = strlen( s1->rseq );
  len2 = strlen( s2->rseq );
  for ( i = 0; ((i < len1) && (i < len2)); i++ ) {
    if ( (s1->rqual[i] >= qcut) &&
	 (s2->rqual[i] >= qcut) ) {
      if ( s1->rseq[i] == s2->rseq[i] ) {
	(*match)++;
      }
      else {
	(*mismatch)++;
	if ( *mismatch > mmc ) {
	  return;
	}
      }
    }
  }
}


/* 
   FOR: ACGTGCATGCTAGACT
   REV: CGATGCTAGTCTAGCA

   then REVCOM(REV): TGCTAGACTAGCATCG
                     |||||||||
         FOR: ACGTGCATGCTAGACT

   Therefore, the overlap would be 9. Ignore any
   base that has Quality score less than QCUT

*/
void compute_ol(const SQP sqp) {
  size_t i, len1, len2, pos;
  char rc_rev[ SEQ_LEN + 1 ];
  char r_rev_q[ SEQ_LEN + 1 ];
  for( i = 0; i < sqpdb->num_reads; i++ ) {
    len1 = strlen( sqpdb->sqps[i]->fseq );
    len2 = strlen( sqpdb->sqps[i]->rseq );
    strcpy( rc_rev, sqpdb->sqps[i]->rseq );
    strcpy( r_rev_q, sqpdb->sqps[i]->rqual );
    revcom_seq( rc_rev );
    rev_qual( r_rev_q );
    sqpdb->sqps[i]->for_rev_ol = 0;
    /* Try each possible starting position 
       on the forward sequence */
    for( pos = 0; pos < len1; pos++ ) {
      if ( perf_match( &(sqpdb->sqps[i]->fseq[pos]),
		       &(sqpdb->sqps[i]->fqual[pos]),
		       rc_rev, r_rev_q ) ) {
	sqpdb->sqps[i]->for_rev_ol = (len1 - pos);
	break;
      }
    }
  }
  return;
}

/* perf_match
   Args: pointer to forward seq,
         pointer to forward qual scores
	 pointer to rev seq,
	 pointer to rev qual scores
   This is the comparison function for finding the overlap
   between the forward and reverse reads. It's called at 
   all possible overlapping positions, from longest to
   shortest, until it finds one. It doesn't require a match
   if either read has quality score less that QCUT.
   Returns: 1 if it's a match, 0 if it's not
*/
int perf_match( const char* s1, const char* q1, size_t len1, 
		const char* s2, const char* q2, size_t len2 ) {
  size_t i, adj_q_cut;
  adj_q_cut = 33 + QCUT;
  for( i = 0; ((i < len1) && (i <len2)); i++ ) {
    if ( (q1[i] >= adj_q_cut) &&
	 (q2[i] >= adj_q_cut) &&
	 (s1[i] != s2[i]) ) {
      return 0;
    }
  }
  return 1;
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
