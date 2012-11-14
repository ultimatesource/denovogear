#include <math.h>

#ifndef GLF_H_
#define GLF_H_


// imported from old glfTools

#define expPhred(x) (double)exp((double)(-(x))/4.343)
#define logPhred(x) (int)((x) < 1 ? (0.5-4.343*log(x)) : (-0.5-4.343*log(x)))

static char iupac[16] = {'-','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'} ;
static int isHom[16] = {0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0} ;
static int isHet[16] = {0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0} ;

/* glz genotype order is: AA/AC/AG/AT/CC/CG/CT/GG/GT/TT, or AMRWCSYGKT in IUPAC */
static int glfBase[10] = { 1, 3, 5, 9, 2, 6, 10, 4, 12, 8 } ; /* mapping from 10 genotypes to 4 bit base coding */
static int baseGlf[16] = { -1, 0, 4, 1, 7, 2, 5, -1, 9, 3, 6, -1, 8, -1, -1, -1 } ; /* inverse of glfBase */

/* utility functions in glfMain.c */

typedef struct {
	unsigned char ref_base:4, dummy:4; /** "XACMGRSVTWYHKDBN"[ref_base] gives the reference base */
	unsigned char max_mapQ; /** maximum mapping quality */
	unsigned char lk[10];   /** log likelihood ratio, capped at 255 */
	unsigned min_lk:8, depth:24; /** minimum lk capped at 255, and the number of mapped reads */
} glf1_t;

#include <stdint.h>
#include "bgzf.h"
typedef BGZF *glfFile;

#define GLF3_RTYPE_END   0
#define GLF3_RTYPE_SUB   1
#define GLF3_RTYPE_INDEL 2

typedef struct {
	uint8_t ref_base:4, rtype:4; /** "XACMGRSVTWYHKDBN"[ref_base] gives the reference base */
	uint8_t rms_mapQ; /** RMS mapping quality */
	uint8_t lk[10];   /** log likelihood ratio, capped at 255 */
	uint32_t min_lk:8, depth:24; /** minimum lk capped at 255, and the number of mapped reads */
	int32_t offset; /** the first base in a chromosome has offset zero. */
	// for indel (lkHom1, lkHom2 and lkHet are the first three elements in lk[10])
	int16_t indel_len[2];
	int32_t max_len; // maximum indel len; will be modified by glf3_read1()
	char *indel_seq[2];
} glf3_t;

typedef struct {
	int32_t l_text;
	uint8_t *text;
} glf3_header_t;

#ifdef __cplusplus
extern "C" {
#endif

#define glf3_init1() ((glf3_t*)calloc(1, sizeof(glf3_t)))
#define glf3_destroy1(g3) do { free((g3)->indel_seq[0]); free((g3)->indel_seq[1]); free(g3); } while (0)

	glf3_header_t *glf3_header_init();
	glf3_header_t *glf3_header_read(glfFile fp);
	void glf3_header_write(glfFile fp, const glf3_header_t *h);
	void glf3_header_destroy(glf3_header_t *h);
	char *glf3_ref_read(glfFile fp, int *len);
	void glf3_ref_write(glfFile fp, const char *name, int len);
	int glf3_write1(glfFile fp, const glf3_t *g3);
	int glf3_read1(glfFile fp, glf3_t *g3);

#ifdef __cplusplus
}
#endif

#endif

/* some common functions: */

extern int qAddTable[1024] ;
#define qAdd(x,y)  (x - qAddTable[512+y-x])
int quality (glf3_t *g);
extern int qsumbit ;

/* commands */

int glfDump (int argc, char *argv[]) ;
int glfExtract (int argc, char *argv[]) ;
int glfSnpCall (int argc, char *argv[]) ;
int glfSubCall (int argc, char *argv[]) ;
int glfSoloPrior (int argc, char *argv[]) ;
int glfStats (int argc, char *argv[]) ;
int glfCheckGenotype (int argc, char *argv[]) ;
