
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>


using namespace std;
#include "bcf.h"
#include "bgzf.h"
#include "denovogear.h"

#define WANT_STREAM       // include iostream and iomanipulators
#include "newmatap.h"
#include "newmatio.h"
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
  	cout<<"Usage\n\tAutosomes\n\n\t\t./denovogear dnm auto --bcf bcf_f --ped ped_f ";
  	cout<<"\n\n\tX in male offspring\n\n\t\t./denovogear dnm XS --bcf bcf_f --ped ped_f ";
  	cout<<"\n\n\tX in female offspring\n\n\t\t./denovogear dnm XD --bcf bcf_f --ped ped_f ";
	cout<<"\n\n\tPhaser\n\n\t\t./denovogear phaser --dnm dnm_F --pgt pgts_f --bam bam_f --window INT[1000]";
	cout<<endl;
}

int findDenovo_auto(string ped_file, string bcf_file, double snp_mrate, 
               double indel_mrate, double poly_rate, double pair_mrate, 
               double mu_scale, string op_vcf_f) 
{
	cout<<"\nPED file : "<<ped_file<<", BCF file : "<<bcf_file; 
  
	// Parse PED files and read in trios
	Trio* trios;
    Pair* pairs;
    int trio_count = 0, pair_count = 0;
    parse_ped (ped_file, &trios, &pairs, trio_count, pair_count); 
	cout<<"\nThe number of trios in the ped file : "<<trio_count<<endl;
    cout<<"The number of paired samples in the ped file : "<<pair_count<<endl<<endl;
    
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

    // Create Paired lookup
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

    //create output vcf
    ofstream fo_vcf;
    if(op_vcf_f != "EMPTY") {
    	fo_vcf.open(op_vcf_f.c_str());
    	if(fo_vcf == NULL) { 
      		cerr<<"Unable to open vcf file for writing output. Exiting !"<<endl;
    	}
    	fo_vcf<<"##fileformat=VCFv4.1\n";
    	fo_vcf<<"##source=DeNovoGear\n";
    	fo_vcf<<"##INPUT - BCF:"<<bcf_file<<" PED:"<<ped_file<<"\n";
    	fo_vcf<<"##INFO=<ID=RD_MOM,Number=1,Type=Integer,Description=\"Read depth for MOM\">\n";
    	fo_vcf<<"##INFO=<ID=RD_DAD,Number=1,Type=Integer,Description=\"Read depth for DAD\">\n";
    	fo_vcf<<"##INFO=<ID=RD_NORMAL,Number=1,Type=Integer,Description=\"Read depth for Normal sample in Paired Sample Analysis\">\n";
    	fo_vcf<<"##INFO=<ID=MQ_MOM,Number=1,Type=Integer,Description=\"Mapping quality for MOM\">\n";
    	fo_vcf<<"##INFO=<ID=MQ_DAD,Number=1,Type=Integer,Description=\"Mapping quality for DAD\">\n";
    	fo_vcf<<"##INFO=<ID=MQ_NORMAL,Number=1,Type=Integer,Description=\"Mapping quality for Normal sample in Paired Sample Analysis\">\n";
    	fo_vcf<<"##INFO=<ID=NULL_CONFIG,Number=1,Type=String,Description=\"NULL trio configuration\">\n";
    	fo_vcf<<"##INFO=<ID=PP_NULL,Number=1,Type=Float,Description=\"Posterior probability for the NULL configuration\">\n";
    	fo_vcf<<"##INFO=<ID=ML_NULL,Number=1,Type=Float,Description=\"Maximum Likelihood for the NULL configuration\">\n";
    	fo_vcf<<"##INFO=<ID=ML_DNM,Number=1,Type=Float,Description=\"Maximum Likelihood for the DNM configuration\">\n";
    	fo_vcf<<"##INFO=<ID=SNPcode,Number=1,Type=Integer,Description=\"Code that shows whether DNM configuration shown is monomorphic or contains variation\">\n";
    	fo_vcf<<"##INFO=<ID=INDELcode,Number=1,Type=Integer,Description=\" \">\n"; //// Needs Detail
    	fo_vcf<<"##INFO=<ID=PairSNPcode,Number=1,Type=Integer,Description=\" \">\n"; //// Needs Detail
    	fo_vcf<<"##INFO=<ID=code,Number=1,Type=Integer,Description=\" \">\n"; //// Needs Detail
    	fo_vcf<<"##FORMAT=<ID=DNM_CONFIG,Number=1,Type=String,Description=\"DNM trio configuration\">\n";
    	fo_vcf<<"##FORMAT=<ID=PP_DNM,Number=1,Type=Float,Description=\"Posterior probability for the DNM configuration\">\n";
    	fo_vcf<<"##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Read Depth of the child\">\n";
    	fo_vcf<<"##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality of the child\">\n";
    	fo_vcf<<"##FORMAT=<ID=RD_T,Number=1,Type=Integer,Description=\"Read Depth of the tumor sample\">\n";
    	fo_vcf<<"##FORMAT=<ID=MQ_T,Number=1,Type=Integer,Description=\"Mapping quality of the normal sample\">\n";
  	}

	// Iterate each position of BCF file
    qcall_t mom_snp, dad_snp, child_snp;
	indel_t mom_indel, dad_indel, child_indel;
	pair_t tumor, normal;	
	bcf_hdr_t *hout, *hin;
	bcf_t *bp = NULL;
	bp = vcf_open(bcf_file.c_str(), "rb");
	hin = hout = vcf_hdr_read(bp);
	bcf1_t *b;
	b = static_cast<bcf1_t *> (calloc(1, sizeof(bcf1_t))); 
	int firstLine = 0; 
	while ( vcf_read(bp, hin, b) > 0 ) {
		int j = 0, flag =0;
		//printf("Position Number %d", pos++);
		// PROCESS TRIOS
		for ( j=0; j<trio_count; j++) {
			int is_indel = bcf_2qcall(hout, b, trios[j],  &mom_snp, &dad_snp, 
									  &child_snp, &mom_indel, &dad_indel, 
									  &child_indel, flag);
			if ( is_indel == 0 ) {				
				if(firstLine == 0) {
					fo_vcf<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"<<child_snp.id<<"\n";
					firstLine = 1;
				}
				trio_like_snp(child_snp, mom_snp, dad_snp, flag, tgt, lookup, op_vcf_f, fo_vcf);   			
			}   
			else if ( is_indel == 1 ) {
				if(firstLine == 0) {
					fo_vcf<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"<<child_indel.id<<"\n";
					firstLine = 1;
				}
				trio_like_indel(&child_indel, &mom_indel, &dad_indel, flag, 
								tgtIndel, lookupIndel, mu_scale, op_vcf_f, fo_vcf);   
			}
			else if ( is_indel < 0 ) {
				printf("\nBCF PARSING ERROR !  %d", is_indel);
				printf("\nExiting !\n\n");
				exit(1);
			}
		}
		  
		// PROCESS  PAIRS
		for ( j=0; j<pair_count; j++) {
			int is_indel = bcf2Paired (hout, b, pairs[j], &tumor, &normal, flag);
			if ( is_indel == 0 ) {
				if(firstLine == 0) {
					fo_vcf<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"<<tumor.id<<"\n";
					firstLine = 1;
				}
				pair_like (tumor, normal, tgtPair, lookupPair, flag, op_vcf_f, fo_vcf);
			}	
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
	fo_vcf.close();
	return 0;
   
}

int findDenovo_XS(string ped_file, string bcf_file, double snp_mrate, 
               double indel_mrate, double poly_rate, double pair_mrate, 
               double mu_scale, string op_vcf_f)
{
	cout<<"\n\n\tXS Model";
	cout<<"\nPED file : "<<ped_file<<", BCF file : "<<bcf_file; 
  
	// Parse PED files and read in trios
	Trio* trios;
    Pair* pairs;
    int trio_count = 0, pair_count = 0;
    parse_ped (ped_file, &trios, &pairs, trio_count, pair_count); 
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
	makeXSSNPLookup(snp_mrate, poly_rate, tgt, lookup);	
	cout<<"\nCreated SNP lookup table - XS\n";
	//Check to see if lookup looks OK
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
	makeXSIndelLookup(indel_mrate, poly_rate, tgtIndel, lookupIndel);
	cout<<"\nCreated indel lookup table - XS";
	//cerr <<endl<<" First mrate "<< lookupIndel.mrate(1,1) <<" Last "<< lookupIndel.mrate(9,3)<<endl;
    cerr <<endl<<" First code "<< lookupIndel.code(1,1) <<" Last "<< lookupIndel.code(9,3)<<endl;
    cerr <<" First tgt "<< tgtIndel[0][0] <<" Last "<< tgtIndel[8][2] <<endl;
    cerr <<" First prior "<< lookupIndel.priors(1,1) <<" Last "<< lookupIndel.priors(9,3) <<endl;

    //create output vcf
    ofstream fo_vcf;

    // Iterate each position of BCF file
    qcall_t mom_snp, dad_snp, child_snp;
	indel_t mom_indel, dad_indel, child_indel;
	bcf_hdr_t *hout, *hin;
	bcf_t *bp = NULL;
	bp = vcf_open(bcf_file.c_str(), "rb");
	hin = hout = vcf_hdr_read(bp);
	bcf1_t *b;
	b = static_cast<bcf1_t *> (calloc(1, sizeof(bcf1_t))); 
	int firstLine = 0; 
	while ( vcf_read(bp, hin, b) > 0 ) {
		int j = 0, flag =0;
		//printf("Position Number %d", pos++);
		// PROCESS TRIOS
		for ( j=0; j<trio_count; j++) {
			int is_indel = bcf_2qcall(hout, b, trios[j],  &mom_snp, &dad_snp, 
									  &child_snp, &mom_indel, &dad_indel, 
									  &child_indel, flag);
			if ( is_indel == 0 ) {
				if(firstLine == 0) {
					fo_vcf<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"<<child_snp.id<<"\n";
					firstLine = 1;
				}
				trio_like_snp(child_snp, mom_snp, dad_snp, flag, tgt, lookup, op_vcf_f, fo_vcf);   			
			}   
			else if ( is_indel == 1 ) {
				if(firstLine == 0) {
					fo_vcf<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"<<child_indel.id<<"\n";
					firstLine = 1;
				}
				trio_like_indel(&child_indel, &mom_indel, &dad_indel, flag, 
								tgtIndel, lookupIndel, mu_scale, op_vcf_f, fo_vcf);  
			}
			else if ( is_indel < 0 ) {
				printf("\nBCF PARSING ERROR !  %d", is_indel);
				printf("\nExiting !\n\n");
				exit(1);
			}
		}
	}

    bcf_hdr_destroy(hin);
	bcf_destroy(b);
	bcf_close(bp);
	fo_vcf.close();
	return 0;

    
}	         

int findDenovo_XD(string ped_file, string bcf_file, double snp_mrate, 
               double indel_mrate, double poly_rate, double pair_mrate, 
               double mu_scale, string op_vcf_f)
{
	cout<<"\n\n\tXD Model";
	cout << "\nPED file : "<<ped_file<<", BCF file : "<<bcf_file; 
  
	// Parse PED files and read in trios
	Trio* trios;
    Pair* pairs;
    int trio_count = 0, pair_count = 0;
    parse_ped (ped_file, &trios, &pairs, trio_count, pair_count); 
	cout << "\nThe number of trios in the ped file : "<<trio_count<<endl;
    cout << "The number of paired samples in the ped file : "<<pair_count<<endl<<endl;

    // Create SNP lookup
	vector<vector<string > > tgt;
	vector<string> tmpV;
	for (int l=0; l<10; l++)
		tmpV.push_back("NA");
	for (int l=0; l<100; l++)
		tgt.push_back(tmpV);
	lookup_snp_t lookup;
	makeXDSNPLookup(snp_mrate, poly_rate, tgt, lookup);	
	cout<<"\nCreated SNP lookup table - XD\n";
	//Check to see if lookup looks OK
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
	makeXDIndelLookup(indel_mrate, poly_rate, tgtIndel, lookupIndel);
	cout<<"\nCreated indel lookup table - XD";
	//cerr <<endl<<" First mrate "<< lookupIndel.mrate(1,1) <<" Last "<< lookupIndel.mrate(9,3)<<endl;
    cerr <<endl<<" First code "<< lookupIndel.code(1,1) <<" Last "<< lookupIndel.code(9,3)<<endl;
    cerr <<" First tgt "<< tgtIndel[0][0] <<" Last "<< tgtIndel[8][2] <<endl;
    cerr <<" First prior "<< lookupIndel.priors(1,1) <<" Last "<< lookupIndel.priors(9,3) <<endl;

    //create output vcf
    ofstream fo_vcf;

    // Iterate each position of BCF file
    qcall_t mom_snp, dad_snp, child_snp;
	indel_t mom_indel, dad_indel, child_indel;
	bcf_hdr_t *hout, *hin;
	bcf_t *bp = NULL;
	bp = vcf_open(bcf_file.c_str(), "rb");
	hin = hout = vcf_hdr_read(bp);
	bcf1_t *b;
	b = static_cast<bcf1_t *> (calloc(1, sizeof(bcf1_t)));
	int firstLine = 0;  
	while ( vcf_read(bp, hin, b) > 0 ) {
		int j = 0, flag =0;
		//printf("Position Number %d", pos++);
		// PROCESS TRIOS
		for ( j=0; j<trio_count; j++) {
			int is_indel = bcf_2qcall(hout, b, trios[j],  &mom_snp, &dad_snp, 
									  &child_snp, &mom_indel, &dad_indel, 
									  &child_indel, flag);
			if ( is_indel == 0 ) {
				trio_like_snp(child_snp, mom_snp, dad_snp, flag, tgt, lookup, op_vcf_f, fo_vcf); 
				if(firstLine == 0) {
					fo_vcf<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"<<child_snp.id<<"\n";
					firstLine = 1;
				}   			
			}   
			else if ( is_indel == 1 ) {
				trio_like_indel(&child_indel, &mom_indel, &dad_indel, flag, 
								tgtIndel, lookupIndel, mu_scale, op_vcf_f, fo_vcf);   
				if(firstLine == 0) {
					fo_vcf<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"<<child_indel.id<<"\n";
					firstLine = 1;
				}
			}
			else if ( is_indel < 0 ) {
				printf("\nBCF PARSING ERROR !  %d", is_indel);
				printf("\nExiting !\n\n");
				exit(1);
			}
		}
	}

    bcf_hdr_destroy(hin);
	bcf_destroy(b);
	bcf_close(bp);
	fo_vcf.close();
	return 0;
}


int mainDNG(int argc, char *argv[])
{
	string ped_f = "EMPTY"; // ped file
	string bcf_f = "EMPTY"; // bcf file
	string op_vcf_f = "EMPTY"; // vcf output -- optional

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
			{"poly_rate", 1, 0, 4}, {"pair_mrates", 1, 0, 5}, {"indel_mu_scale", 1, 0, 6}, {"output", 1, 0, 7}};
		int c = getopt_long (argc-1, argv+1, "", long_options, &option_index);
		if (c == -1)
			break;
		switch(c) 
		{      
			case 0:
				ped_f = optarg;
				break;
			case 1:
				bcf_f = optarg;
				break;
			case 2:
				snp_mrate = atof(optarg);
				cerr<<"\nThe new point mutation rate is : "<<snp_mrate;
				break;
			case 3:
				indel_mrate = atof(optarg);
				cerr<<"\nThe new indel mutation rate is : "<<indel_mrate;
				break;
			case 4:
				poly_rate = atof(optarg);
				cerr<<"\nThe new polymorphism rate is : "<<poly_rate;
				break;
			case 5:
				pair_mrate = atof(optarg);
				cerr<<"\nThe new paired-sample mutation rate is : "<<pair_mrate;
				break;
			case 6:
				mu_scale = atof(optarg);
				cerr<<"\nThe new indel mutation rate scaling factor rate is : "<<mu_scale;
				break;
			case 7:
				op_vcf_f = optarg;
				cerr<<"\nThe output vcf file is : "<<op_vcf_f;
				break;
			default:
				cerr<<"\nUnrecognized option ! "<<argv[c]<<" Exiting !";
				usage();
    			exit(1);
		}
	}
  
	if((bcf_f == "EMPTY") || (ped_f == "EMPTY")) {
		cerr<<"ERROR ! Please specify both the PED file and BCF file !"
        "Exiting!\n";
		exit(1);
	}
  
    // Create lookup table and read ped, BCF files. Separate for autosomes, X in males and X in females.
    if (strcmp(argv[1], "auto") == 0)
		findDenovo_auto(ped_f, bcf_f, snp_mrate, indel_mrate, poly_rate, pair_mrate, mu_scale, op_vcf_f);
	else if (strcmp(argv[1], "XS") == 0)
		findDenovo_XS(ped_f, bcf_f, snp_mrate, indel_mrate, poly_rate, pair_mrate, mu_scale, op_vcf_f);
	else if (strcmp(argv[1], "XD") == 0)
		findDenovo_XD(ped_f, bcf_f, snp_mrate, indel_mrate, poly_rate, pair_mrate, mu_scale, op_vcf_f);
	else {
    	usage();
    	exit(1);
  	}

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


