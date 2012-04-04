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