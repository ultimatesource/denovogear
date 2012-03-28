#include <vector>
using namespace std;

#ifndef LOOKUP_H
#include "lookup.h"
#endif

void printUsage();
void writeLookup(ofstream& fout, bool snp);
void getIndelPriors(string g_gts1, int n_uniqa);
void getSNPPriors(string g_gts1, int n_uniqa);
void makeSNPLookup(double SNPMrate, double PolyRate, 
	vector<vector<string > > & tgt, lookup_snp_t & lookup);
void makeIndelLookup(double IndelMrate, double PolyRate, 
	vector<vector<string > > & tgt, lookup_indel_t & lookupIndel);
void makePairedLookup(double pairMrate, vector<vector<string > > & tgt, lookup_pair_t & lookup);
void readIndelLookup(ofstream& fout, vector<vector<string > > & tgt, 
	float lines[][27]);
void readSNPLookup(ofstream& fout, vector<vector<string > > & tgt,
	 float lines[][1000]); 