#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>


#include "bcf.h"

#include "bgzf.h"
#include "denovogear.h"

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
#define MIN_READ_DEPTH 10
#define MIN_READ_DEPTH_INDEL 10
#define MIN_MAPQ 40 


int parse(char* ped_file, char* bcf_file, double snp_mrate, double indel_mrate,
		  double poly_rate, double mu_scale) 
{
	printf("\nPED file : %s, BCF file : %s", ped_file, bcf_file);   
	// Parse PED files and read in trios
    Trio trios[ MAX_TRIOS ];
    int trio_count = parse_ped(ped_file, trios); 
	printf ("\nThe number of trios in the ped file : %d\n\n", trio_count);
	
    // Create SNP lookup
	vector<vector<string > > tgt;
	vector<string> tmpV;
	for (int l=0;l<10;l++)
		tmpV.push_back("NA");
	for (int l=0;l<100;l++)
		tgt.push_back(tmpV);
	lookup_t lookup;
	makeLookup("point", snp_mrate, indel_mrate, poly_rate, tgt, lookup);	
	cout<<"\nRead SNP lookup table\n";
	cerr <<" First mrate "<< lookup.mrate(1,1) <<" Last "<< lookup.mrate(100,10)<<endl;
	cerr <<" First code "<< lookup.code(1,1) <<" Last "<< lookup.code(100,10)<<endl;
	cerr <<" First tgt "<< tgt[0][0] <<" Last "<< tgt[99][9] <<endl;
	cerr <<" First tref "<< lookup.tref(1,1) <<" Last "<< lookup.tref(100,10) <<endl;

    // Create INDEL lookup
	vector<vector<string > > tgtIndel;
	vector<string> tmpIndel;
    for (int l=0;l<3;l++)
        tmpIndel.push_back("NA");
    for (int l=0;l<9;l++)
        tgtIndel.push_back(tmpIndel);
	lookup_t lookupIndel;
	makeLookup("indel", snp_mrate, indel_mrate, poly_rate, tgtIndel, lookupIndel);
	cout<<"\nRead indel lookup table";
	cerr <<endl<<" First mrate "<< lookupIndel.mrate(1,1) <<" Last "<< lookupIndel.mrate(9,3)<<endl;
    cerr <<" First code "<< lookupIndel.code(1,1) <<" Last "<< lookupIndel.code(9,3)<<endl;
    cerr <<" First tgt "<< tgtIndel[0][0] <<" Last "<< tgtIndel[8][2] <<endl;
    cerr <<" First prior "<< lookupIndel.priors(1,1) <<" Last "<< lookupIndel.priors(9,3) <<endl;

	// Iterate each position of BCF file
    qcall_t mom_snp, dad_snp, child_snp;
	indel_t mom_indel, dad_indel, child_indel;	
	bcf_hdr_t *hout, *hin;
	bcf_t *bp = vcf_open(bcf_file, "rb");
	hin = hout = vcf_hdr_read(bp);
	bcf1_t *b;
	b = static_cast<bcf1_t *> (calloc(1, sizeof(bcf1_t)));  
	while ( vcf_read(bp, hin, b) > 0 ) {
		int j = 0, flag =0;
		//printf("Position Number %d", pos++);
		for ( j=0; j<trio_count; j++) {
			int is_indel = bcf_2qcall(hout, b, trios[j],  &mom_snp, &dad_snp, 
									  &child_snp, &mom_indel, &dad_indel, 
									  &child_indel, flag);
			if ( is_indel == 0 ) {
				trio_like_snp(child_snp, mom_snp, dad_snp, flag, tgt, lookup);   			
			}   
			else if ( is_indel == 1 ) {
				trio_like_indel(&child_indel, &mom_indel, &dad_indel, flag, 
								tgtIndel, lookupIndel, mu_scale);  
			}
			else if ( is_indel < 0 ) {
				printf("\n BCF PARSING ERROR !  %d", is_indel);
				printf("\n Exiting !");
				exit(1);
			}
		}   		  
	}

	bcf_hdr_destroy(hin);
	bcf_destroy(b);
	bcf_close(bp);
	return 0;
   
}


int main(int argc, char *argv[])
{
	char ped_f[50] = "EMPTY", bcf_f[50] = "EMPTY";
	double indel_mrate = 1e-9, snp_mrate = 1e-8, poly_rate = 1e-3, mu_scale = 1.0;
    
    // Read in Command Line arguments
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {{"ped", 1, 0, 0}, 
			{"bcf", 1, 0, 1}, {"snp_mrate", 1, 0, 2}, {"indel_mrate", 1, 0, 3}, 
			{"poly_rate", 1, 0, 4}, {"indel_mu_scale", 1, 0, 5}};
		int c = getopt_long (argc, argv, "", long_options, &option_index);
		if (c == -1)
			break;
		switch(c) 
		{      
			case 0:
				strcpy(ped_f, optarg);
				break;
			case 1:
				strcpy(bcf_f, optarg);
				break;
			case 2:
				snp_mrate = atof(optarg);
				break;
			case 3:
				indel_mrate = atof(optarg);
				break;
			case 4:
				poly_rate = atof(optarg);
				break;
			case 5:
				mu_scale = atof(optarg);
				break;
		}
	}
	if(!strcmp(bcf_f, "EMPTY") || !strcmp(ped_f, "EMPTY")) {
		cout<<"ERROR ! Please specify both the PED file and BCF file ! Exiting!\n";
		exit(1);
	}
  
    // Create lookup table and read ped, BCF files.
	parse(ped_f, bcf_f, snp_mrate, indel_mrate, poly_rate, mu_scale);
	cerr<<"\nDone !\n";
	return 0;
}
