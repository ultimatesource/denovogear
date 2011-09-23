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
#define MIN_READ_DEPTH 0
#define MIN_READ_DEPTH_INDEL 10
#define MIN_MAPQ 40 

Trio trios[ MAX_TRIOS ];

int autLength[24] = { 0, 247249719, 242951149, 199501827, 191263066, 180837866, 170896992, 158821424, 146274818, 140273252, 135374737, 134452384, 132289542, 114127979, 106360583, 100338915, 88822253, 78654741, 76117153, 63806651, 62435958, 46944323, 49591432, 154913754} ;
int X_length = 154913754 ;


void trio_like_snp(qcall_t child, qcall_t mom, qcall_t dad, int flag, vector<vector<string> > &tgt, lookup_t & lookup);
void trio_like_indel(indel_t *child, indel_t *mom, indel_t *dad, int flag, vector<vector<string> > &tgtIndel, lookup_t & lookupIndel);
int read_lookup (vector<vector<string> > &tgt, lookup_t & lookup);
int read_indelLookup (vector<vector<string> > &tgtIndel, lookup_t & lookupIndel);

