/*
 * Copyright (c) 2010, 2011 Genome Research Ltd.
 * Copyright (c) 2012, 2013 Donald Conrad and Washington University in St. Louis
 * Authors: Donald Conrad <dconrad@genetics.wustl.edu>,
 * Avinash Ramu <aramu@genetics.wustl.edu>
 * This file is part of DeNovoGear.
 *
 * DeNovoGear is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "bcf.h"
#include "bgzf.h"
#include "parser.h"
#include "makeLookup.h"

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


// Calculate SNP DNM PP
void trio_like_snp( qcall_t child, qcall_t mom, qcall_t dad, int flag,
  vector<vector<string > > & tgt, lookup_snp_t & lookup,
  string op_vcf_f, ofstream& fo_vcf, double pp_cutoff, int RD_cutoff, int& n_site_pass);

// Calculate INDEL DNM PP
void trio_like_indel(indel_t *child,indel_t *mom, indel_t *dad, int flag,
                     vector<vector<string > > & tgtIndel,
                     lookup_indel_t & lookupIndel, double mu_scale,
                     string op_vcf_f, ofstream& fo_vcf, double pp_cutoff,
                     int RD_cutoff, int& n_site_pass, double user_indel_mrate);

// Calculate Pair PP
void pair_like(pair_t tumor, pair_t normal, vector<vector<string> > &tgtPair,
	       lookup_pair_t & lookupPair, int flag, string op_vcf_f, ofstream& fo_vcf,
         double pp_cutoff, int RD_cutoff, int& n_site_pass);

