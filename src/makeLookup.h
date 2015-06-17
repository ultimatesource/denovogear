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

#ifndef MAKE_LOOKUP_H
#define MAKE_LOOKUP_H

#include <vector>
#include "lookup.h"


void printUsage();
void writeLookup(std::ofstream &fout, bool snp);
void getIndelPriors(std::string g_gts1, int n_uniqa);
void getSNPPriors(std::string g_gts1, int n_uniqa);

// Lookup for paired sample analysis
void makePairedLookup(double pairMrate, lookup_table_t &tgt, 
		      lookup_pair_t &lookup);

void setIndelLines(lookup_table_t &tgt, float lines[][27]);

void setSNPLines(lookup_table_t &tgt,
                 float lines[][1000]);

// SNP Lookup for the autosome model
void makeSNPLookup(double SNPMrate, double PolyRate,
                   lookup_table_t &tgt, lookup_snp_t &lookup);

// Indel Lookup for the autosome model
void makeIndelLookup(double PolyRate,
                     lookup_table_t &tgt, lookup_indel_t &lookupIndel);

// SNP Lookup for the XS model
void makeXSSNPLookup(double SNPMrate, double PolyRate,
                     lookup_table_t &tgt, lookup_snp_t &lookup);

// Indel Lookup for the XS model
void makeXSIndelLookup(double PolyRate,
                       lookup_table_t &tgt, lookup_indel_t &lookupIndel);

// SNP Lookup for the XD model
void makeXDSNPLookup(double SNPMrate, double PolyRate,
                     lookup_table_t &tgt, lookup_snp_t &lookup);

// Indel Lookup for the XD model
void makeXDIndelLookup(double PolyRate,
                       lookup_table_t &tgt, lookup_indel_t &lookupIndel);

#endif
