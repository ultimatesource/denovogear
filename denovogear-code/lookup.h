#include <math.h>
#define WANT_STREAM       // include iostream and iomanipulators
#include "newmatap.h"
#include "newmatio.h"

using namespace std;

// Struct to hold Lookup Table
typedef struct {
  Matrix aref; /*priors for "A" reference allele */
  Matrix cref; /*priors for "C" reference allele */
  Matrix gref; /*priors for "G" reference allele */
  Matrix tref; /*priors for "T" reference allele */
  Matrix priors; /*priors for Indel */
  Matrix snpcode; /*code for identifying a biallelic locus */
  Matrix tp; /*code for identifying a biallelic locus */
  Matrix code; /*code for identifying a biallelic locus */
  Matrix mrate; /*mutation probability*/
  Matrix denovo; /*code for identifying de novo events + mutation rate*/
  Matrix norm; /*code for identifying de novo events + mutation rate*/
  Matrix hit; /* multiplication factor for indel mu */
} lookup_t;

