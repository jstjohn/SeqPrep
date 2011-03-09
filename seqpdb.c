#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "seqpdb.h"

/* init_SQPDB
   Initialize and return a SQPDB
   Args: (1) size_t size - how big the array of sqps should be
*/
SQPDB init_SQPDB( size_t size ) {
  size_t i;
  SQPDB sqpdb;
  SQP first_seq;

  /* Try to allocate memory for SQPDB */
  sqpdb            = (SQPDB)malloc(sizeof(Sqpdb));
  if ( sqpdb == NULL ) {
    fprintf( stderr, "Not enough memories for database" );
    return NULL;
  }

  /* Try to allocate memory for sqps array */
  first_seq = (SQP)malloc( size * sizeof(Sqp) );
  sqpdb->sqps = (SQP*)malloc( size * sizeof(SQP) );
  
  if ( (first_seq   == NULL) ||
       (sqpdb->sqps == NULL) ) {
    fprintf( stderr, "Not enough memories for database" );
    return NULL;
  } 
    
  for( i = 0; i < size; i++ ) {
    sqpdb->sqps[i] = &first_seq[i];
  }

  sqpdb->seq_len   = SEQ_LEN;
  sqpdb->num_reads = 0;
  sqpdb->size      = size;

  return sqpdb;
}

/* grow_sqpdb
 */
SQPDB grow_sqpdb( SQPDB sqpdb ) {
  size_t new_size, i;
  SQP first_new_seq;
  SQP* new_sqps;

  new_size = sqpdb->size * 2;
  
  /* Try to allocate new memory for 2x sqps array */
  new_sqps = (SQP*)malloc(new_size * sizeof(SQP) );
  first_new_seq = (SQP)malloc(sqpdb->size * sizeof(Sqp) );

  if ( (first_new_seq == NULL) ||
       (new_sqps == NULL) ) {
    fprintf( stderr, "Not enough memories for database\n" );
    return NULL;
  }

  for( i = 0; i < sqpdb->size; i++ ) {
    new_sqps[i              ] = sqpdb->sqps[i];
    new_sqps[i + sqpdb->size] = &first_new_seq[i];
  }
  free( sqpdb->sqps );
  sqpdb->sqps = new_sqps;
  sqpdb->size = new_size;
  return sqpdb;
}

/* read_fastqs
   Takes the filenames of two fastq files, these are the forward
   and reverse reads in the same order. Also takes the sqpdb, which
   may be caused to grow while reading in the fastq data.
   Returns the final read status. This is 0 if everything went
   smoothly to completion.
*/
int read_fastqs( char* ffqfn, char* rfqfn, SQPDB sqpdb ) {
  char fid[MAX_ID_LEN + 1];
  char rid[MAX_ID_LEN + 1];
  int read_status;
  FILE* ffq;
  FILE* rfq;

  ffq = fileOpen( ffqfn, "r" );
  rfq = fileOpen( rfqfn, "r" );
  
  read_status = next_fastqs( ffq, rfq, sqpdb );
  while( read_status == 1 ) {
    read_status = next_fastqs( ffq, rfq, sqpdb );
  }

  return read_status;
}

/* next_fastqs
   Read the next forward and reverse fastq sequences.
   Check to make sure their ID's are compatible and
   put the results in the next SQP of SQPDB. Grow
   this, if necessary.
*/
int next_fastqs( FILE* ffq, FILE* rfq, SQPDB sqpdb ) {
  char fid[MAX_ID_LEN + 1];
  char rid[MAX_ID_LEN + 1];
  int frs; // forward fastq read status
  int rrs; // reverse fastq read status
  SQP curr_sqp;

  /* Check if sqpdb is big enough to accomadate another read
     pair. If not, grow it. */
  if ( sqpdb->size == sqpdb->num_reads ) {
    sqpdb = grow_sqpdb( sqpdb );
  }
  
  curr_sqp = sqpdb->sqps[sqpdb->num_reads];
 
  /* Read the next fastq record from the forward and reverse
     pair of each */
  frs = read_fastq( ffq, fid, curr_sqp->fseq, curr_sqp->fqual );
  rrs = read_fastq( rfq, rid, curr_sqp->rseq, curr_sqp->rqual );

  if ( (frs == 1) &&
       (rrs == 1) &&
       f_r_id_check( fid, rid ) ) {
    /* Populate the qual_sum field */
    calc_qual_sum( curr_sqp );
    sqpdb->num_reads++;
    return 1;
  }

  else {
    return 0;
  }
}

int f_r_id_check( char fid[], char rid[] ) {
  size_t fid_len, rid_len;

  fid_len = strlen( fid );
  rid_len = strlen( rid );

  fid[fid_len - 2] = '\0'; //don't care what the last 2 characters are
  rid[rid_len - 2] = '\0'; //don't care what the last 2 characters are

  if (strcmp( fid, rid ) == 0 ) {
    return 1;
  }
  return 0;
}

/* Takes an SQP and calculates the sum of quality
   scores for all bases on the forward and reverse
   reads. Ass-u-mes a ASCII offset of 33.
   Sets the qual_sum field and returns nothing */
void calc_qual_sum( SQP s ) {
  size_t i, len;
  s->qual_sum = 0;
  len = strlen( s->fqual );
  for(i = 0; i < len; i++ ) {
    s->qual_sum += s->fqual[i] - 33;
  }
  len = strlen( s->rqual );
  for(i = 0; i < len; i++ ) {
    s->qual_sum += s->rqual[i] - 33;
  }
}      

/* read_fastq
   Return 1 => more sequence to be had
          0 => EOF
 */
int read_fastq( FILE* fastq, char id[], char seq[], char qual[] ) {
  char c;
  size_t i;
  c = fgetc( fastq );
  if ( c == EOF ) return 0;
  if ( c != '@' ) {
    fprintf( stderr, "fastq record not beginning with @\n" );
    return 0;
  }

  /* get identifier */
  i = 0;
  while( (!isspace(c=fgetc( fastq ) ) &&
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

  /* Now, everything else on the line is description (if anything)
     although fastq does not appear to formally support description */
  while ( (c != '\n') &&
	  (c != EOF) ) {
    c = fgetc( fastq );
  }
  i = 0;

  /* Now, read the sequence. This should all be on a single line */
  i = 0;
  c = fgetc( fastq );
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
    c = fgetc( fastq );
  }
  seq[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     INIT_ALN_SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == SEQ_LEN ) {
    while ( (c != '\n') &&
	    (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  /* Now, read the quality score header */
  c = fgetc( fastq );
  if ( c != '+' ) {
    fprintf( stderr, "Problem reading quality line for %s\n", id );
    return 1;
  }
  /* Zip through the rest of the line, it should be the same identifier
     as before or blank */
  c = fgetc( fastq );
  while( (c != '\n') &&
	 (c != EOF) ) {
    c = fgetc( fastq );
  }

  /* Now, get the quality score line */
  c = fgetc( fastq );
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
    c = fgetc( fastq );
  }
  qual[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     INIT_ALN_SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == SEQ_LEN ) {
    while ( (c != '\n') &&
	    (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  if ( c == EOF ) {
    return 0;
  }
  return 1;
}


/** fileOpen **/
FILE * fileOpen(const char *name, char access_mode[]) {
  FILE * f;
  f = fopen(name, access_mode);
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

/* sort_sqps
   Sorts the sqps in an sqpdb by their
   quality sum. Uses qual_sum_comp
*/
void sort_sqps( SQPDB sqpdb ) {
  qsort( &sqpdb->sqps[0], sqpdb->num_reads,
	 sizeof(SQP), qual_sum_comp );
}

/* Takes a sorted sqpdb and a quality-sum cutoff value.
   Sets the num_reads field such that all the reads 
   below the quality-sum cutoff value are excluded */
void set_qsum_cut( SQPDB sqpdb, int cutoff ) {
  size_t i = 0;
  while( sqpdb->sqps[i]->qual_sum >= cutoff ) {
    i++;
  }
  sqpdb->num_reads = (i-1);
  return;
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
void compute_ol( const SQPDB sqpdb ) {
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
int perf_match( const char* s1, const char* q1,
		const char* s2, const char* q2 ) {
  size_t len1, len2, i, adj_q_cut;
  adj_q_cut = 33 + QCUT;
  len1 = strlen( s1 );
  len2 = strlen( s2 );
  for( i = 0; ((i < len1) && (i <len2)); i++ ) {
    if ( (q1[i] >= adj_q_cut) &&
	 (q2[i] >= adj_q_cut) &&
	 (s1[i] != s2[i]) ) {
      return 0;
    }
  }
  return 1;
}

/* Computes and populates the num_skf and num_skr fields
   in for all entries in the input sqpdb, given the kmer
   counts in ks and the cutoffs specified by lo_k and hi_k.
*/
void compute_sk( const SQPDB sqpdb, const KSP ks, 
		 const int lo_k, const int hi_k ) {
  size_t i = 0;
  size_t len, pos, inx;
  int k;

  k = ks->k;
  /* Go through each sequence pair */
  for( i = 0; i < sqpdb->num_reads; i++ ) {
    /* Initialize sweet-spot kmers to 0 */
    sqpdb->sqps[i]->num_skf = 0;
    sqpdb->sqps[i]->num_skr = 0;

    len = strlen( sqpdb->sqps[i]->fseq );
    /* Look at each kmer. Is it in the sweet spot? */
    for( pos = 0; pos < (len - k + 1); pos++ ) {
      if (kmer2inx( &sqpdb->sqps[i]->fseq[pos], k, &inx )) {
	if ( (ks->ka[inx]->num_occur >= lo_k) &&
	     (ks->ka[inx]->num_occur <= hi_k) ) {
	  sqpdb->sqps[i]->num_skf += 1;
	}
      }
    }

    len = strlen( sqpdb->sqps[i]->rseq );
    /* Look at each kmer. Is it in the sweet spot? */
    for( pos = 0; pos < (len - k + 1); pos++ ) {
      if (kmer2inx( &sqpdb->sqps[i]->rseq[pos], k, &inx )) {
	if ( (ks->ka[inx]->num_occur >= lo_k) &&
	     (ks->ka[inx]->num_occur <= hi_k) ) {
	  sqpdb->sqps[i]->num_skr += 1;
	}
      }
    }
  }
  return;
}

/* Takes an sqpdb with valid data in the sqps[]->num_skf|r fields
   Calculates and prints a two dimensional matrix of the number
   of sweet-spot kmers in forward and reverse reads (the row and
   column position) over all reads. */
void output_sweet_matrix( const SQPDB sqpdb ) {
  size_t i, j;
  size_t sm[SEQ_LEN][SEQ_LEN];
  int max_dim = 0;
  /* Zero out matrix */
  for( i = 0; i < SEQ_LEN; i++ ) {
    for( j = 0; j < SEQ_LEN; j++ ) {
      sm[i][j] = 0;
    }
  }

  /* Add in correct values */
  for( i = 0; i < sqpdb->num_reads; i++ ) {
    sm[ sqpdb->sqps[i]->num_skf ][ sqpdb->sqps[i]->num_skr ]++;
    if ( sqpdb->sqps[i]->num_skf > max_dim ) {
      max_dim = sqpdb->sqps[i]->num_skf;
    }
    if ( sqpdb->sqps[i]->num_skr > max_dim ) {
      max_dim = sqpdb->sqps[i]->num_skr;
    }
  }
  
  /* Print it out */
  for( i = 0; i < max_dim; i++ ) {
    for( j = 0; j < max_dim; j++ ) {
      printf( "%d ", (int)sm[i][j] );
    }
    printf( "\n" );
  }
}

int qual_sum_comp( const void* sqp1, const void* sqp2 ) {
  SQP* s1p = (SQP*) sqp1;
  SQP* s2p = (SQP*) sqp2;
  SQP s1 = *s1p;
  SQP s2 = *s2p;
  if ( s1->qual_sum < s2->qual_sum ) {
    return 1;
  }
  if ( s1->qual_sum > s2->qual_sum ) {
    return -1;
  }
  return 0;
}

int qual_p_cut( SQPDB sqpdb, int p ) {
  unsigned int qsum_hist[2 * SEQ_LEN * 40];
  size_t i;
  size_t total = 0;
  /* Initialize */
  for( i = 0; i < (SEQ_LEN*40); i++ ) {
    qsum_hist[i] = 0;
  }

  for( i = 0; i < sqpdb->num_reads; i++ ) {
      qsum_hist[ sqpdb->sqps[i]->qual_sum ]++;
  }
  
  for( i = 0; i < (SEQ_LEN*40*2); i++ ) {
    total += qsum_hist[i];

    if ( (((double)total / (double)sqpdb->num_reads) * 100.0) >=
	 (100 - p) ) {
      return i;
    }
  }
}

void output_qual_sum_hist( SQPDB sqpdb ) {
  unsigned int qsum_hist[SEQ_LEN * 40];
  size_t i;

  /* Initialize */
  for( i = 0; i < (SEQ_LEN*40); i++ ) {
    qsum_hist[i] = 0;
  }

  for( i = 0; i < sqpdb->num_reads; i++ ) {
    if ( sqpdb->sqps[i]->qual_sum >= (SEQ_LEN * 40) ) {
      printf( "Holy crap.\n" );
    }
    else {
      qsum_hist[ sqpdb->sqps[i]->qual_sum ]++;
    }
  }
  printf( "# QSUM NUM_OCCURANCES\n" );
  for( i = 0; i < (SEQ_LEN*40); i++ ) {
    printf( "%d %u\n", (int)i, qsum_hist[i] );
  }
  return;
}

/* Initialize a new KSP
   Populate all ka array elements with NULL
*/
KSP init_KSP( int k ) {
  KSP ks;
  size_t i, len;
  len = 1<<(k*2);
  ks = (KSP)malloc(sizeof(Kmers));

  ks->k = k;
  ks->ka = (KOP*)malloc(sizeof(KOP) * len);
  for( i = 0; i < len; i++ ) {
    ks->ka[i] = NULL;
  }
  return ks;
}

KOP init_KOP( void ) {
  KOP ko;
  ko = (KOP)malloc(sizeof(Kmer_occur));
  ko->num_occur = 0;
  return ko;
}

/* Takes a pointer to the kmers struct, the index for
   the kmer we've observed and the sequence number (sn)
   for the sequence it occured in. If this number is
   negative, it denotes occurance in the reverse read.
   Adds this occurance to the occur field and increments
   the num_occur field. If we've already seen the
   MAX_KMER_OCCUR for this kmer, we do nothing.
   Also, do nothing if the last sequence added was this one.
   This happens if a kmer occurs multiple times in a sequence.
   We only care that it occurs once.
 */
int add_to_ksp( KSP ks, size_t inx, int sn ) {
  size_t pos;
  /* Ever see this kmer (inx) before? */
  if ( ks->ka[inx] == NULL ) {
    ks->ka[inx] = (KOP)init_KOP();
    ks->ka[inx]->occur[0] = sn;
    ks->ka[inx]->num_occur = 1;
    return 1;
  }

  pos = ks->ka[inx]->num_occur;
  if ( pos >= MAX_KMER_OCCUR ) {
    return -1;
  }

  if ( ks->ka[inx]->occur[pos - 1] == sn ) {
    return -2;
  }
  ks->ka[inx]->occur[pos] = sn;
  ks->ka[inx]->num_occur += 1;
  return 1;
}

/* kmer2inx
   Args: (1) a pointer to a character string;
             the kmer to find the corresponding index of;
	     might not be null-terminated
	 (2) length of the kmer
	 (3) pointer to size_t to put the index
   Returns: TRUE if the index was set, FALSE if it could not
            be set because of some non A,C,G,T character
   Uses the formula A=>00, C=>01, G=>11, T=>11 to make a
   bit string for the kmer. Any other character is not allowed
   and will cause an error
   The bit string is constructed by reading the kmer from left
   to right. This bit-string is then interpreted as a variable
   of type size_t and is appropriate as an array index
*/
int kmer2inx( const char* kmer,
	      const unsigned int kmer_len,
	      size_t* inx ) {
  size_t l_inx = 0;
  int i = 0;
  char curr_char;

  while( i < kmer_len ) {
    l_inx = l_inx << 2;
    curr_char = toupper(kmer[i]); // Upper case it in case it is not
    switch( curr_char ) {
    case 'A' :
      l_inx += 0;
      break;
    case 'C' :
      l_inx += 1;
      break;
    case 'G' :
      l_inx += 2;
      break;
    case 'T' :
      l_inx += 3;
      break;
    default :
      return 0; // not valid!
    }
    i++;
  }
  *inx = l_inx;
  return 1; // valid!
}

void output_read_kmer_cont_dist( KSP ks ) {
  size_t read_dist[MAX_KMER_OCCUR + 1];
  size_t i;
  size_t num_kmers;
  num_kmers = 1<<(2 * ks->k);

  for( i = 0; i < (MAX_KMER_OCCUR + 1); i++ ) {
    read_dist[i] = 0;
  }
  for( i = 0; i < num_kmers; i++ ) {
    if ( ks->ka[i] == NULL ) {
      ;
    }
    else {
      read_dist[ ks->ka[i]->num_occur ] += ks->ka[i]->num_occur;
    }
  }
  printf( "# kmer_occurance_number number_of_seqs\n" );
  for( i = 0; i < (MAX_KMER_OCCUR+1); i++ ) {
    printf( "%d %d\n", (int)i, (int)read_dist[i] );
  }
}
	

void output_kmer_dist( KSP ks ) {
  /* dist[] => index is the number of occurances for a kmer
               value is the number of kmers that occur that
	       number of times */
  size_t dist[MAX_KMER_OCCUR + 1];
  size_t i;
  size_t num_kmers;
  num_kmers = 1<<(2 * ks->k);

  for( i = 0; i < (MAX_KMER_OCCUR + 1); i++ ) {
    dist[i] = 0;
  }
  for( i = 0; i < num_kmers; i++ ) {
    
    if ( ks->ka[i] == NULL ) {
      dist[0]++;
    }

    else {
      dist[ ks->ka[i]->num_occur ]++;
    }
  }

  for( i = 0; i < (MAX_KMER_OCCUR+1); i++ ) {
    printf( "%d %d\n", (int)i, (int)dist[i] );
  }
  return;
}
  
OLP init_OLP( int max_ol ) {
  OLP olp;

  olp = (OLP)malloc(sizeof(Overlappers));
  
  olp->ols = (int*)malloc(sizeof(int) * max_ol);
  olp->pos = (int*)malloc(sizeof(int) * max_ol);
  olp->n = 0;

  return olp;
}


void revcom_seq( char seq[] ) {
  char tmp_base;
  int len, i;
  len = strlen(seq);

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
    return 'N';
  }
}

void rev_qual( char q[] ) {
  char tmp_q;
  int len, i;
  len = strlen(q);
  
  for (i = 0; i < len/2; i++) {
    tmp_q = q[i];
    q[i] = q[len-(i+1)];
    q[len-(i+1)] = tmp_q;
  }
  
  /* If sequence length is even, we're done, otherwise there is
     the quality score right in the center, But we leave it
     anyway here...
  */
}
