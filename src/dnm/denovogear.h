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

#ifndef DENOVOGEAR_H_
#define DENOVOGEAR_H_

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "parser.h"
#include "lookup.h"
#include "makeLookup.h"


//#ifdef use_namespace
//using namespace RBD_LIBRARIES;
//#endif


// Calculate SNP DNM PP
int trio_like_snp(qcall_t &child, qcall_t &mom, qcall_t &dad, int flag,
                   lookup_table_t &tgt, lookup_snp_t &lookup,
                   std::vector<hts::bcf::File> &vcfout, double pp_cutoff, int RD_cutoff);

// Calculate INDEL DNM PP
int trio_like_indel(indel_t &child, indel_t &mom, indel_t &dad, int flag,
                     lookup_table_t &tgtIndel, lookup_indel_t &lookupIndel, double mu_scale,
                     std::vector<hts::bcf::File> &vcfout, double pp_cutoff,
                     int RD_cutoff, double user_indel_mrate);

// Calculate Pair PP
void pair_like(pair_t &tumor, pair_t &normal,
               lookup_table_t &tgtPair, lookup_pair_t &lookupPair,
               int flag, std::vector<hts::bcf::File> &vcfout,
               double pp_cutoff, int RD_cutoff, int &n_site_pass);


#endif

