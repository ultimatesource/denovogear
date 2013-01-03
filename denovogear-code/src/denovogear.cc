
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
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
#define VERSION "Denovogear0.5.2"
#endif


void usage()
{
  	cout<<VERSION<<endl;
  	cout<<"Usage\n\tAutosomes\n\n\t\t./denovogear dnm auto --bcf bcf_f --ped ped_f ";
  	cout<<"\n\n\tX in male offspring\n\n\t\t./denovogear dnm XS --bcf bcf_f --ped ped_f ";
  	cout<<"\n\n\tX in female offspring\n\n\t\t./denovogear dnm XD --bcf bcf_f --ped ped_f ";
	cout<<"\n\n\tPhaser\n\n\t\t./denovogear phaser --dnm dnm_F --pgt pgts_f --bam bam_f --window INT[1000]";
	cout<<endl;
	exit(1);
}

int callMakeSNPLookup(vector<vector<string > >& tgtSNP, lookup_snp_t& lookupSNP, const double& snp_mrate, const double& poly_rate, string model)
{
	vector<string> tmp;
	for (int l=0; l<10; l++)
		tmp.push_back("NA");
	for (int l=0; l<100; l++)
		tgtSNP.push_back(tmp);
	

	if (model == "auto") makeSNPLookup(snp_mrate, poly_rate, tgtSNP, lookupSNP);
	else if (model == "XS") makeXSSNPLookup(snp_mrate, poly_rate, tgtSNP, lookupSNP);
	else if (model == "XD") makeXDSNPLookup(snp_mrate, poly_rate, tgtSNP, lookupSNP);
	else {
		cerr<<endl<<"Invalid model for SNP lookup, exiting.";
		exit(1);
	}
	
	cerr<<"\nCreated SNP lookup table\n";
	cerr <<" First mrate: "<< lookupSNP.mrate(1,1) <<" last: "<< lookupSNP.mrate(100,10)<<endl;
	cerr <<" First code: "<< lookupSNP.code(1,1) <<" last: "<< lookupSNP.code(100,10)<<endl;
	cerr <<" First target string: "<< tgtSNP[0][0] <<" last: "<< tgtSNP[99][9] <<endl;
	cerr <<" First tref: "<< lookupSNP.tref(1,1) <<" last: "<< lookupSNP.tref(100,10) <<endl;
	return 0;
}

int callMakeINDELLookup(vector<vector<string > >& tgtIndel, lookup_indel_t& lookupIndel, const double& indel_mrate, const double& poly_rate, string model)
{
	vector<string> tmp;
	for (int l=0; l<3; l++)
        tmp.push_back("NA");
    for (int l=0; l<9; l++)
        tgtIndel.push_back(tmp);
	

	if (model == "auto") makeIndelLookup(indel_mrate, poly_rate, tgtIndel, lookupIndel);
	else if (model == "XS")	makeXSIndelLookup(indel_mrate, poly_rate, tgtIndel, lookupIndel);
	else if (model == "XD") makeXDIndelLookup(indel_mrate, poly_rate, tgtIndel, lookupIndel);
	else {
			cerr<<endl<<"Invalid model for INDEL lookup, exiting.";
			exit(1);
	}
	
	cerr<<"\nCreated indel lookup table";
	cerr <<" First code: "<< lookupIndel.code(1,1) <<" last: "<< lookupIndel.code(9,3)<<endl;
    cerr <<" First target string: "<< tgtIndel[0][0] <<" last: "<< tgtIndel[8][2] <<endl;
    cerr <<" First prior: "<< lookupIndel.priors(1,1) <<" last: "<< lookupIndel.priors(9,3) <<endl;	
    return 0;
}

int callMakePairedLookup(vector<vector<string > >& tgtPair, lookup_pair_t lookupPair, const double& pair_mrate, const double& poly_rate)
{
	vector<string> tmp;
	for (int l=0; l<10; l++)
        tmp.push_back("NA");
    for (int l=0; l<10; l++)
        tgtPair.push_back(tmp);
    
    makePairedLookup(pair_mrate, tgtPair, lookupPair);
    cerr<<"\nCreated paired lookup table"<<endl;
    cerr <<" First target string: "<< tgtPair[0][0] <<" last: "<< tgtPair[9][9] <<endl;
    cerr <<" First prior "<< lookupPair.priors(1,1) <<" last: "<< lookupPair.priors(10,10) <<endl;
    return 0;
}

int writeVCFHeader(std::ofstream& fo_vcf, string op_vcf_f, string bcf_file, string ped_file, string sample)
{
	fo_vcf.open(op_vcf_f.c_str());
	if(fo_vcf == NULL) { 
  		cerr<<"Unable to open vcf file for writing output. Exiting !"<<endl;
	}
	fo_vcf<<"##fileformat=VCFv4.1\n";
	fo_vcf<<"##source= "<<VERSION<<"\n";
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
	fo_vcf<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"<<sample<<"\n";
	return 0;
}

int callDenovoFromBCF(string ped_file, string bcf_file, double snp_mrate, 
               double indel_mrate, double poly_rate, double pair_mrate, 
               double mu_scale, string op_vcf_f, string model)
{
	

    // Create SNP lookup
    lookup_snp_t lookupSNP;
	vector<vector<string > > tgtSNP;
	callMakeSNPLookup(tgtSNP, lookupSNP, snp_mrate, poly_rate, model);	

    // Create INDEL lookup
    lookup_indel_t lookupIndel;
	vector<vector<string > > tgtIndel;
	callMakeINDELLookup(tgtIndel, lookupIndel, indel_mrate, poly_rate, model);    
    
    // Create paired lookup
    lookup_pair_t lookupPair;
    vector<vector<string > > tgtPair;
    callMakePairedLookup(tgtPair, lookupPair, pair_mrate, poly_rate); 

	// Iterate each position of BCF file
	Trio* trios;
    Pair* pairs;
    int trio_count = 0, pair_count = 0;
    parse_ped (ped_file, &trios, &pairs, trio_count, pair_count);

	//create output vcf -- assumes theres only one trio/pair
	string sample;
	if(trio_count > 0)
  		sample = trios[0].cID;
  	else
  		sample = pairs[0].tumorID;
    ofstream fo_vcf;
    if(op_vcf_f != "EMPTY") {
    	writeVCFHeader(fo_vcf, op_vcf_f, bcf_file, ped_file, sample);
  	}
  	
    qcall_t mom_snp, dad_snp, child_snp;
	indel_t mom_indel, dad_indel, child_indel;
	pair_t tumor, normal;	
	bcf_hdr_t *hout, *hin;
	bcf_t *bp = NULL;
	bp = vcf_open(bcf_file.c_str(), "rb");
	hin = hout = vcf_hdr_read(bp);
	bcf1_t *b;
	b = static_cast<bcf1_t *> (calloc(1, sizeof(bcf1_t))); 

	while ( vcf_read(bp, hin, b) > 0 ) {
		int j = 0, flag =0;
		//printf("Position Number %d", pos++);
		// PROCESS TRIOS
		for ( j=0; j<trio_count; j++) {
			int is_indel = bcf_2qcall(hout, b, trios[j],  &mom_snp, &dad_snp, 
									  &child_snp, &mom_indel, &dad_indel, 
									  &child_indel, flag);
			if ( is_indel == 0 ) {				
				trio_like_snp(child_snp, mom_snp, dad_snp, flag, 
								tgtSNP, lookupSNP, op_vcf_f, fo_vcf);   			
			}   
			else if ( is_indel == 1 ) {
				trio_like_indel(&child_indel, &mom_indel, &dad_indel, flag, 
								tgtIndel, lookupIndel, mu_scale, op_vcf_f, fo_vcf);   
			}
			else if ( is_indel < 0 ) {
				printf("\n BCF PARSING ERROR - Paired Sample!  %d\n Exiting !\n", is_indel);
				exit(1);
			}
		}
		  
		// PROCESS  PAIRS
		if(model == "auto") { // paired sample model not developed for XS, XD
			for ( j=0; j<pair_count; j++) {
				int is_indel = bcf2Paired (hout, b, pairs[j], &tumor, &normal, flag);
				if ( is_indel == 0 ) {
					pair_like (tumor, normal, tgtPair, lookupPair, flag, op_vcf_f, fo_vcf);
				}	
				else if ( is_indel < 0 ) {
					printf("\n BCF PARSING ERROR - Paired Sample!  %d\n Exiting !\n", is_indel);
					exit(1);
				}	
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
	string model;
	if(argc > 1)
		model = argv[1];
	else 
		usage();
	string ped_f = "EMPTY"; // ped file
	string bcf_f = "EMPTY"; // bcf file
	string op_vcf_f = "EMPTY"; // vcf output -- optional

	double indel_mrate = 1e-9; // indel mutation prior
	double snp_mrate = 1e-8; // snp mutation prior
	double poly_rate = 1e-3; // polymorphism rate - used in prior calculations
	double pair_mrate = 1e-9; // mutation prior for paired samples
	double mu_scale = 1.0; // scaling factor for indel priors

	float posterior_prob_cutoff = 0.0001;
	int min_RD = 10; // read depth cutoff
    
    // Read in Command Line arguments
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {{"ped", 1, 0, 0}, 
			{"bcf", 1, 0, 1}, {"snp_mrate", 1, 0, 2}, {"indel_mrate", 1, 0, 3}, 
			{"poly_rate", 1, 0, 4}, {"pair_mrates", 1, 0, 5}, {"indel_mu_scale", 1, 0, 6}, {"output_vcf", 1, 0, 7},
			{"posterior_prob_cutoff", 1, 0, 8}, {"read_depth_cutoff", 1, 0, 9}, {"h", 1, 0, 10}};
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
				cerr<<"\nSNP mutation rate: "<<snp_mrate;
				break;
			case 3:
				indel_mrate = atof(optarg);
				cerr<<"\nindel mutation rate: "<<indel_mrate;
				break;
			case 4:
				poly_rate = atof(optarg);
				cerr<<"\npolymorphism rate: "<<poly_rate;
				break;
			case 5:
				pair_mrate = atof(optarg);
				cerr<<"\npaired-sample mutation rate: "<<pair_mrate;
				break;
			case 6:
				mu_scale = atof(optarg);
				cerr<<"\nindel mutation rate scaling factor: "<<mu_scale;
				break;
			case 7:
				op_vcf_f = optarg;
				cerr<<"\noutput vcf file: "<<op_vcf_f;
				break;
			case 8:
				posterior_prob_cutoff = atof(optarg);
				cerr<<"\nposterior probability cutoff: "<<posterior_prob_cutoff;
				break;
			case 9:
				min_RD = atoi(optarg);
				cerr<<"\nread depth filter: "<<min_RD;
				break;
			default:
				usage();
    			exit(1);
		}
	}
  
	if ((bcf_f == "EMPTY") || (ped_f == "EMPTY")) {
		cerr<<"ERROR ! Please specify both the PED file and BCF file !"
        "Exiting!\n";
        usage();
	}  
  	
    // Create lookup table and read ped, BCF files. Separate for autosomes, X in males and X in females.
    if (model ==  "auto" || model == "XS" || model == "XD"){
    	callDenovoFromBCF(ped_f, bcf_f, snp_mrate, indel_mrate, poly_rate, pair_mrate, mu_scale, op_vcf_f, model);
	}
	else {
    	usage();
  	}

	cerr<<"\nDone !"<<endl;
	exit(0);
}


int main(int argc, char* argv[])
{
  cerr<<VERSION;  
  if(argc >= 2)
  	if (strcmp(argv[1], "dnm") == 0) return mainDNG(argc-1, argv+1);
  	else if (strcmp(argv[1], "phaser") == 0) return mainPhaser(argc-1, argv+1);
  usage();
}


