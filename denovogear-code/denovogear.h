#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "bcf.h"
#include "bgzf.h"
#include "lookup.h"
#include "parser.h"

#define WANT_STREAM       // include iostream and iomanipulators
#include "newmatap.h"
#include "newmatio.h"

using namespace std;

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#ifdef use_namespace
using namespace RBD_LIBRARIES;
#endif


#ifndef VERSION
#define VERSION "dummy"
#endif

#define MRATE 5e-7
#define MIN_MAPQ 40 

// Calculate SNP DNM PP
void trio_like_snp(qcall_t child, qcall_t mom, qcall_t dad, int flag, 
				   vector<vector<string> > &tgt, lookup_t & lookup);
// Calculate INDEL DNM PP
void trio_like_indel(indel_t *child, indel_t *mom, indel_t *dad, int flag, 
					 vector<vector<string> > &tgtIndel, lookup_t & lookupIndel, double mu_scale);
// Make SNP and INDEL lookup tables
void makeLookup(string table_type, double Mrate, double IndelMrate, 
				double PolyRate, vector<vector<string > > & tgt, 
				lookup_t & lookup);
