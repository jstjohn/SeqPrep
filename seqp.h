#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <zlib.h>

#define MAX_ID_LEN (256)
#define SEQ_LEN (64)
#define MAX_KMER_OCCUR (128)
#define QCUT (20) // Quality score cutoff for mismatching bases to determine if there is overlap in forward and reverse reads
#define MAX_ALLOW_FOR_REV_OL (4) // Maximum number of base-pairs allowed to overlap to be used in making scontigs

/* Type to hold the forward and reverse read
   of a sequence pair with quality scores */
typedef struct sqp {
  char fseq[SEQ_LEN+1];
  char fqual[SEQ_LEN+1];
  char rseq[SEQ_LEN+1];
  char rqual[SEQ_LEN+1];
  int qual_sum;
  int num_skf; // number of sweep-spot kmers on forward
  int num_skr; // number of sweep-spot kmers on reverse
  int for_rev_ol; // number of bases that (might) overlap in the
                  // forward and reverse reads
} Sqp;
typedef struct sqp* SQP;

/* Type to hold an array of sqp, i.e, a database
   of sequences pairs and their quality scores */
typedef struct sqpdb {
  SQP* sqps; // pointer to array of seq & qual pairs
  unsigned int seq_len; // the length of the sequences
  size_t num_reads; // The current length of the database, 0-indexed
  size_t size; // The currently allocated size of the database
} Sqpdb;
typedef struct sqpdb* SQPDB;

/* Type to hold an array of kmers. The array index
   specifies what the kmer is using this scheme:
   A=>00, C=>01, G=>10, T=>11
   A bit string is this composed, from 5' to 3' of
   the kmer. The bitstring is interpreted as a size_t
   to index the array
   Example: ACTCG => 0001110110b => 118d
*/
typedef struct kmer_occur {
  int occur[MAX_KMER_OCCUR];
  size_t num_occur;
} Kmer_occur;
typedef struct kmer_occur* KOP;
typedef struct kmers {
  int k; // length of kmers
  KOP* ka; // pointer to size_t; the kmer array
} Kmers;
typedef struct kmers* KSP;

typedef struct overlappers {
  size_t n; // number of overlappers
  int* ols; // pointer to list of overlappers
  int* pos; // alignment positions for overlappers (offset
            // relative to the overlappee => 0 means they
            // start at the same place
} Overlappers;
typedef struct overlappers* OLP;



/* Prototypes */
SQPDB init_SQPDB( size_t size );
SQPDB grow_sqpdb( SQPDB sqpdb );
int read_fastqs( char* ffqfn, char* rfqfn, SQPDB sqpdb );
int next_fastqs( FILE* ffq, FILE* rfq, SQPDB sqpdb );
int f_r_id_check( char fid[], char rid[] );
void calc_qual_sum( SQP s );
int read_fastq( FILE* fastq, char id[], char seq[], char qual[] );
FILE * fileOpen(const char *name, char access_mode[]);
int qual_p_cut( SQPDB sqpdb, int p );
void output_qual_sum_hist( SQPDB sqpdb );
void sort_sqps( SQPDB sqpdb );
int qual_sum_comp( const void* sqp1, const void* sqp2 );
void comp_seq( const SQP s1, const SQP s2, const int qcut, 
	       const int mmc, int* match, int* mismatch );
KSP init_KSP( int k );
KOP init_KOP( void );
int add_to_ksp( KSP ks, size_t inx, int sn );
int kmer2inx( const char* kmer,
	      const unsigned int kmer_len,
	      size_t* inx );
void output_kmer_dist( KSP ks );
void output_read_kmer_cont_dist( KSP ks );
void compute_sk( const SQPDB sqpdb, const KSP ks, 
		 const int lo_k, const int hi_k );
void set_qsum_cut( SQPDB sqpdb, int cutoff );
void revcom_seq( char seq[] );
inline char revcom_char(const char base);
OLP init_OLP( int max_ol );
void rev_qual( char q[] );
int perf_match( const char* s1, const char* q1,
		const char* s2, const char* q2 );
void compute_ol( const SQPDB sqpdb );
