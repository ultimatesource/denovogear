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

#include <vector>
using namespace std;

#ifndef LOOKUP_H
#include "lookup.h"
#endif

void printUsage();
void writeLookup(ofstream& fout, bool snp);
void getIndelPriors(string g_gts1, int n_uniqa);
void getSNPPriors(string g_gts1, int n_uniqa);

// Lookup for paired sample analysis
void makePairedLookup(double pairMrate, vector<vector<string > > & tgt, lookup_pair_t & lookup);

void setIndelLines(vector<vector<string > > & tgt,
	float lines[][27]);

void setSNPLines(vector<vector<string > > & tgt,
	 float lines[][1000]);

// SNP Lookup for the autosome model
void makeSNPLookup(double SNPMrate, double PolyRate,
	vector<vector<string > > & tgt, lookup_snp_t & lookup);

// Indel Lookup for the autosome model
void makeIndelLookup(double IndelMrate, double PolyRate,
	vector<vector<string > > & tgt, lookup_indel_t & lookupIndel);

// SNP Lookup for the XS model
void makeXSSNPLookup(double SNPMrate, double PolyRate,
	vector<vector<string > > & tgt, lookup_snp_t & lookup);

// Indel Lookup for the XS model
void makeXSIndelLookup(double IndelMrate, double PolyRate,
	vector<vector<string > > & tgt, lookup_indel_t & lookupIndel);

// SNP Lookup for the XD model
void makeXDSNPLookup(double SNPMrate, double PolyRate,
	vector<vector<string > > & tgt, lookup_snp_t & lookup);

// Indel Lookup for the XD model
void makeXDIndelLookup(double IndelMrate, double PolyRate,
	vector<vector<string > > & tgt, lookup_indel_t & lookupIndel);
