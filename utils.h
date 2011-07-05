#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <zlib.h>
#include "stdaln.h"

#define MAX_ID_LEN (256)
#define MAX_FN_LEN (512)
#define MAX_SEQ_LEN (256)
//60+33 = 93 = '[' (was 83='S')
#define MAX_QUAL (93)
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
  size_t merged_len;
  size_t rlen;
  size_t mpos;
} Sqp;
typedef struct sqp* SQP;

SQP SQP_init();
void SQP_destroy(SQP sqp);
void adapter_merge(SQP sqp, bool print_overhang);
void fill_merged_sequence(SQP sqp, AlnAln *aln, bool include_overhang);
void pretty_print_alignment(gzFile out, SQP sqp, char adj_q_cut, bool sort);
void pretty_print_alignment_stdaln(gzFile out, SQP sqp, AlnAln *aln, bool first_adapter, bool second_adapter, bool print_merged);
inline char mismatch_p33_merge(char pA, char pB);
inline char match_p33_merge(char pA, char pB);
void make_blunt_ends(SQP sqp, AlnAln *aln);
bool read_olap_adapter_trim(SQP sqp, size_t min_ol_adapter,
    unsigned short min_match_adapter[MAX_SEQ_LEN+1],
    unsigned short max_mismatch_adapter[MAX_SEQ_LEN+1],
    unsigned short min_match_reads[MAX_SEQ_LEN+1],
    unsigned short max_mismatch_reads[MAX_SEQ_LEN+1],
    char qcut);
bool read_merge(SQP sqp, size_t min_olap,
    unsigned short min_match[MAX_SEQ_LEN+1],
    unsigned short max_mismatch[MAX_SEQ_LEN+1],
    char adj_q_cut);
inline bool next_fastqs( gzFile ffq, gzFile rfq, SQP curr_sqp, bool p64 );
inline int write_fastq(gzFile out, char id[], char seq[], char qual[]);
inline bool f_r_id_check( char fid[], size_t fid_len, char rid[], size_t rid_len );
int read_fastq( gzFile* fastq, char id[], char seq[], char qual[],
    size_t *id_len, size_t *seq_len, bool p64 );
gzFile * fileOpen(const char *name, char access_mode[]);
int compute_ol(
    char subjectSeq[], char subjectQual[], size_t subjectLen,
    char querySeq[], char queryQual[], size_t queryLen,
    size_t min_olap,
    unsigned short min_match[MAX_SEQ_LEN+1],
    unsigned short max_mismatch[MAX_SEQ_LEN+1],
    bool check_unique, char adj_q_cut );
bool k_match( const char* s1, const char* q1, size_t len1,
    const char* s2, const char* q2, size_t len2,
    unsigned short min_match,
    unsigned short max_mismatch, char adj_q_cut);
void revcom_seq( char seq[], int len);
inline char revcom_char(const char base);
inline void rev_qual( char q[], int len );
bool adapter_trim(SQP sqp, size_t min_ol_adapter,
    char *forward_primer, char *forward_primer_dummy_qual,
    int forward_primer_len,
    char *reverse_primer, char *reverse_primer_dummy_qual,
    int reverse_primer_len,
    unsigned short min_match_adapter[MAX_SEQ_LEN+1],
    unsigned short max_mismatch_adapter[MAX_SEQ_LEN+1],
    unsigned short min_match_reads[MAX_SEQ_LEN+1],
    unsigned short max_mismatch_reads[MAX_SEQ_LEN+1],
    char adj_q_cut);

#ifndef max
  #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
  #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif
