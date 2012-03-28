#include <math.h>
#define WANT_STREAM       // include iostream and iomanipulators
#include "newmatap.h"
#include "newmatio.h"

#define MIN_READ_DEPTH_SNP 10


using namespace std;

// SNP Lookup Table
typedef struct {
  Matrix aref; /*priors for "A" reference allele */
  Matrix cref; /*priors for "C" reference allele */
  Matrix gref; /*priors for "G" reference allele */
  Matrix tref; /*priors for "T" reference allele */
  Matrix snpcode; /*code for identifying a biallelic locus */
  Matrix tp; /*code for identifying a biallelic locus */
  Matrix code; /*code for identifying a biallelic locus */
  Matrix mrate; /*mutation probability*/
  Matrix denovo; /*code for identifying de novo events + mutation rate*/
  Matrix norm; /*code for identifying de novo events + mutation rate*/
} lookup_snp_t;

// Indel Lookup Table
typedef struct {
  Matrix priors; /*priors for Indel */
  Matrix snpcode; /*code for identifying a biallelic locus */
  Matrix tp; /*code for identifying a biallelic locus */
  Matrix code; /*code for identifying a biallelic locus */
  Matrix mrate; /*mutation probability*/
  Matrix denovo; /*code for identifying de novo events + mutation rate*/
  Matrix norm; /*code for identifying de novo events + mutation rate*/
  Matrix hit; /* multiplication factor for indel mu */
} lookup_indel_t;

// Paired-sample Lookup Table
typedef struct {
  Matrix snpcode; /*code for identifying a biallelic locus */
  Matrix priors; /*priors */
  Matrix denovo; /*code for identifying de novo events + mutation rate*/
  Matrix norm; /*code for identifying de novo events + mutation rate*/
} lookup_pair_t;
