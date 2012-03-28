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

#include "phaser.h"

#ifdef use_namespace
using namespace RBD_LIBRARIES;
#endif


#ifndef VERSION
#define VERSION "dummy"
#endif

#define MRATE 5e-7

void usage()
{
  	cout<<"DeNovoGear - Identify denovo mutations from next-gen sequencing data.\n";
  	cout<<"Usage \n\t./denovogear dnm --bcf bcf_f --ped ped_f \nOR";
	cout<<"\n\t./denovogear phaser --dnm dnm_F --pgt pgts_f --bam bam_f --window INT[1000]";
	cout<<endl;
}

int findDenovo(char* ped_file, char* bcf_file, double snp_mrate, 
               double indel_mrate, double poly_rate, double pair_mrate, double mu_scale) 
{
	printf("\nPED file : %s, BCF file : %s", ped_file, bcf_file); 
  
	// Parse PED files and read in trios
    Trio trios[ MAX_TRIOS ];
    Pair pairs[ MAX_PAIRS ];
    int trio_count = 0, pair_count = 0;
    parse_ped (ped_file, trios, pairs, trio_count, pair_count); 
	printf ("\nThe number of trios in the ped file : %d", trio_count);
    printf ("\nThe number of paired samples in the ped file : %d\n\n", pair_count);

	
    // Create SNP lookup
	vector<vector<string > > tgt;
	vector<string> tmpV;
	for (int l=0; l<10; l++)
		tmpV.push_back("NA");
	for (int l=0; l<100; l++)
		tgt.push_back(tmpV);
	lookup_snp_t lookup;
	makeSNPLookup(snp_mrate, poly_rate, tgt, lookup);	
	cout<<"\nCreated SNP lookup table\n";
	cerr <<" First mrate "<< lookup.mrate(1,1) <<" Last "<< lookup.mrate(100,10)<<endl;
	cerr <<" First code "<< lookup.code(1,1) <<" Last "<< lookup.code(100,10)<<endl;
	cerr <<" First tgt "<< tgt[0][0] <<" Last "<< tgt[99][9] <<endl;
	cerr <<" First tref "<< lookup.tref(1,1) <<" Last "<< lookup.tref(100,10) <<endl;

    // Create INDEL lookup
	vector<vector<string > > tgtIndel;
	vector<string> tmpIndel;
    for (int l=0; l<3; l++)
        tmpIndel.push_back("NA");
    for (int l=0; l<9; l++)
        tgtIndel.push_back(tmpIndel);
	lookup_indel_t lookupIndel;
	makeIndelLookup(indel_mrate, poly_rate, tgtIndel, lookupIndel);
	cout<<"\nCreated indel lookup table";
	//cerr <<endl<<" First mrate "<< lookupIndel.mrate(1,1) <<" Last "<< lookupIndel.mrate(9,3)<<endl;
    cerr <<" First code "<< lookupIndel.code(1,1) <<" Last "<< lookupIndel.code(9,3)<<endl;
    cerr <<" First tgt "<< tgtIndel[0][0] <<" Last "<< tgtIndel[8][2] <<endl;
    cerr <<" First prior "<< lookupIndel.priors(1,1) <<" Last "<< lookupIndel.priors(9,3) <<endl;

    // Create Pair INDEL lookup
    vector<vector<string > > tgtPair;
	vector<string> tmpPair;
    for (int l=0; l<10; l++)
        tmpPair.push_back("NA");
    for (int l=0; l<10; l++)
        tgtPair.push_back(tmpPair);
    lookup_pair_t lookupPair;
    makePairedLookup(pair_mrate, tgtPair, lookupPair);
    cerr<<"\nCreated paired lookup table"<<endl;
    cerr <<" First tgt "<< tgtPair[0][0] <<" Last "<< tgtPair[9][9] <<endl;
    cerr <<" First prior "<< lookupPair.priors(1,1) <<" Last "<< lookupPair.priors(10,10) <<endl;

	// Iterate each position of BCF file
    qcall_t mom_snp, dad_snp, child_snp;
	indel_t mom_indel, dad_indel, child_indel;
	pair_t tumor, normal;	
	bcf_hdr_t *hout, *hin;
	bcf_t *bp = NULL;
	bp = vcf_open(bcf_file, "rb");
	hin = hout = vcf_hdr_read(bp);
	bcf1_t *b;
	b = static_cast<bcf1_t *> (calloc(1, sizeof(bcf1_t)));  
	while ( vcf_read(bp, hin, b) > 0 ) {
		int j = 0, flag =0;
		//printf("Position Number %d", pos++);
		// PROCESS ALL TRIOS
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
		// PROCESS ALL PAIRS
		for ( j=0; j<pair_count; j++) {
			int is_indel = bcf2Paired (hout, b, pairs[j], &tumor, &normal, flag);
			if ( is_indel == 0 ) 
				pair_like (tumor, normal, tgtPair, lookupPair, flag);	
			else if ( is_indel < 0 ) {
				printf("\n BCF PARSING ERROR - Paired Sample!  %d", is_indel);
				printf("\n Exiting !\n");
				exit(1);
			}	
		}  	
		
			  
	}

	bcf_hdr_destroy(hin);
	bcf_destroy(b);
	bcf_close(bp);
	return 0;
   
}


int mainDNG(int argc, char *argv[])
{
	char ped_f[500] = "EMPTY", bcf_f[500] = "EMPTY";

	double indel_mrate = 1e-9;// indel mutation rate
	double snp_mrate = 1e-8;// snp mutation rate
	double poly_rate = 1e-3;// polymorphism rate - used in prior calculations
	double pair_mrate = 1e-9;// mutation rate in paired samples
	double mu_scale = 1.0;// scaling factor for indel priors
    
    // Read in Command Line arguments
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {{"ped", 1, 0, 0}, 
			{"bcf", 1, 0, 1}, {"snp_mrate", 1, 0, 2}, {"indel_mrate", 1, 0, 3}, 
			{"poly_rate", 1, 0, 4}, {"pair_mrates", 1, 0, 5}, {"indel_mu_scale", 1, 0, 6}};
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
				cout<<"\nThe new point mutation rate is : "<<snp_mrate;
				break;
			case 3:
				indel_mrate = atof(optarg);
				cout<<"\nThe new indel mutation rate is : "<<indel_mrate;
				break;
			case 4:
				poly_rate = atof(optarg);
				cout<<"\nThe new polymorphism rate is : "<<poly_rate;
				break;
			case 5:
				pair_mrate = atof(optarg);
				cout<<"\nThe new paired-sample mutation rate is : "<<pair_mrate;
				break;
			case 6:
				mu_scale = atof(optarg);
				cout<<"\nThe new indel mutation rate scaling factor rate is : "<<mu_scale;
				break;
		}
	}
  
	if(!strcmp(bcf_f, "EMPTY") || !strcmp(ped_f, "EMPTY")) {
		cout<<"ERROR ! Please specify both the PED file and BCF file !"
        "Exiting!\n";
		exit(1);
	}
  
    // Create lookup table and read ped, BCF files.
	findDenovo(ped_f, bcf_f, snp_mrate, indel_mrate, poly_rate, pair_mrate, mu_scale);
	cerr<<"\nDone !\n";
	return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 2) {
    usage();
    exit(1);
  }
  else if (strcmp(argv[1], "dnm") == 0) return mainDNG(argc-1, argv+1);
  else if (strcmp(argv[1], "phaser") == 0) return mainPhaser(argc-1, argv+1);
  else {
    usage();
    exit(1);
  }
  return 0;
}


