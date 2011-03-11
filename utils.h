#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <zlib.h>

#define MAX_ID_LEN (256)
#define MAX_SEQ_LEN (256)
#define CODE_AMBIGUOUS (-2)
#define CODE_NOMATCH (-1)
#define CODE_NOADAPT (9999)

/* Type to hold the forward and reverse read
   of a sequence pair with quality scores */
typedef struct sqp {
  char fid[MAX_ID_LEN+1];
  char fseq[MAX_SEQ_LEN+1];
  char fqual[MAX_SEQ_LEN+1];
  size_t flen;
  char rid[MAX_ID_LEN+1];
  char rseq[MAX_SEQ_LEN+1];
  char rqual[MAX_SEQ_LEN+1];
  char rc_rseq[MAX_SEQ_LEN+1];
  char rc_rqual[MAX_SEQ_LEN+1];
  char merged_seq[MAX_SEQ_LEN+MAX_SEQ_LEN+1];
  char merged_qual[MAX_SEQ_LEN+MAX_SEQ_LEN+1];
  size_t rlen;
  size_t insert_len; // only valid if guessing that primer is present in seq
  size_t for_rev_ol; // number of bases that (might) overlap in the
                     // forward and reverse reads, unless insert_len
} Sqp;
typedef struct sqp* SQP;


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
