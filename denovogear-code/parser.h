#ifndef PARSER_H_
#define PARSER_H_
#include "bcf.h"
#include "pedParser.h"
#define MAX_QCALL_LINE 2048
#define ID_LENGTH 100

#define MIN_READ_DEPTH 0
#define MIN_MAPQ 40 

//int bcf_2qcall(bcf_hdr_t *h, bcf1_t *b, Trio t, vector<vector<string > >  & tgt, lookup_t & lookup);

typedef struct { //New struct for Q calls
	char chr[3]; /* Chromosome Number. */
	long pos; /* position, the first base in a chromosome has offset zero. */
	char ref_base; /* Either A, C, G or T */
	char alt[20]; /* ALT string */
	int depth; /* number of mapped reads */
	int rms_mapQ; /* RMS mapping quality */
	int min_lk; /* minimum lk capped at 255 */	
	int lk[10];   /* log likelihood ratio, capped at 255 */
	char id[ID_LENGTH]; /* string for sample ID */
} qcall_t;

typedef struct { //New struct for Indels
	char chr[3]; /* chromosome Number. */
	long pos; /* position, the first base in a chromosome has offset zero. */
	char ref_base[2048]; /* Either A, C, G or T */
	char alt[2048]; /* ALT string */
	int depth; /* number of mapped reads */
	int rms_mapQ; /* RMS mapping quality */
	int min_lk; /* minimum lk capped at 255 */	
	int lk[3];   /* log likelihood ratio, capped at 255, ORDER = ref/ref, ref/alt & alt/alt */
	char id[ID_LENGTH]; /* string for sample ID */
	int length; /* to be filled in by Ruth */
} indel_t;

int parse(char* ped_file, char* bcf_file);
int bcf_2qcall(bcf_hdr_t *h, bcf1_t *b, Trio t, qcall_t* mom, qcall_t* dad, qcall_t* child, indel_t* mom_indel, indel_t* dad_indel, indel_t* child_indel, int& flag);

#endif