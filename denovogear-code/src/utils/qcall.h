#define ID_LENGTH 100

//New struct for Q calls
typedef struct {
	char chr[3]; /* Chromosome Number. */
	long pos; /* position, the first base in a chromosome has offset zero. */
	char ref_base; /* Either A, C, G or T */
	int depth; /* number of mapped reads */
	int rms_mapQ; /* RMS mapping quality */
	int min_lk; /* minimum lk capped at 255 */	
	int lk[10];   /* log likelihood ratio, capped at 255 */
	char id[ID_LENGTH]; /* string for sample ID */
} qcall_t;