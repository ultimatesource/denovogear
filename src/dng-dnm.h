#ifndef DNG_DNM_H_
#define DNG_DNM_H_

#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <getopt.h>
#include <memory> // Needs to be put in dng/hts/hts.h
#include <assert.h> // Needs to be put in dng/hts/hts.h

#include <dng/hts/bcf.h>

#include "version.h"
#include "lookup.h"
#include "makeLookup.h"


int RD_cutoff = 10; // cutoff for the read depth filter
double PP_cutoff = 0.0001; // posterior probability cutoff

double indel_mrate = 0; // indel mutation prior is based on a log linear model unless specified by the user.
double snp_mrate = 1e-8; // snp mutation prior
double poly_rate = 1e-3; // polymorphism rate - used in prior calculations
double pair_mrate = 1e-9; // mutation prior for paired samples
double mu_scale = 1.0; // scaling factor for indel priors

typedef std::vector<std::vector<std::string> > lookup_table_t;

int callDenovoFromBCF(std::string ped_file, std::string bcf_file, std::string op_vcf_f, 
		      std::string model, bool is_vcf, std::string region);

/*
void usage() {
  std::cerr << "\nUsage:\n";
    std::cerr << "Autosomes:\n";
    std::cerr << "\tdng dnm auto --bcf bcf_f --ped ped_f [OR] dng dnm auto --vcf vcf_f --ped ped_f\n";
    std::cerr << "X chromosome in male offspring:\n";
    std::cerr << "\tdng dnm XS --bcf bcf_f --ped ped_f [OR] dng dnm XS --vcf vcf_f --ped ped_f\n";
    std::cerr << "X chromosome in female offspring:\n";
    std::cerr << "\tdng dnm XD --bcf bcf_f --ped ped_f [OR] dng dnm XD --vcf vcf_f --ped ped_f\n";
    //cerr<<"Phaser:\n";
    //cerr<<"\tdng phaser --dnm dnm_f --pgt pgt_f --bam bam_f --window INT[1000]\n";
    std::cerr << "\nInput:\n";
    std::cerr << "DNM:\n";
    std::cerr << "--ped:\t Ped file to describe relationship between the samples.\n";
    std::cerr << "--bcf:\t BCF file, contains per-sample read depths and genotype likelihoods.\n";
    std::cerr << "--vcf:\t VCF file, contains per-sample read depths and genotype likelihoods.\n";
    std::cerr << "Phaser:\n";
    std::cerr << "--dnm: Tab delimited list of denovo mutations to be phased, format: chr pos inherited_base denovo_base.[example: 1 2000 A C]\n";
    std::cerr << "--pgt: Tab delimited genotypes of child and parents at SNP sites near denovo sites, format: chr pos GT_child GT_parent1 GT_parent2.[example: 1 2000 AC AC AA]\n";
    std::cerr << "--bam: alignment file (.bam) of the child.\n";
    std::cerr << "--window: optional argument which is the maximum distance between the DNM and a phasing site. The default value is 1000.\n";
    std::cerr << "\nOutput:\n";
    std::cerr << "--output_vcf:\t vcf file to write the output to.\n";
    std::cerr << "\nParameters:\n";
    std::cerr << "--snp_mrate:\t Mutation rate prior for SNPs. [1e-8]\n";
    std::cerr << "--indel_mrate:\t Mutation rate prior for INDELs. [1e-9]\n";
    std::cerr << "--pair_mrate:\t Mutation rate prior for paired sample analysis. [1e-9]\n";
    std::cerr << "--indel_mu_scale:\t Scaling factor for indel mutation rate. [1]\n";
    std::cerr << "--pp_cutoff:\t Posterior probability threshold. [0.0001]\n";
    std::cerr << "--rd_cutoff:\t Read depth filter, sites where either one of the sample have read depth less than this threshold are filtered out. [10]\n";
    std::cerr << "--region:\t Region of the BCF file to perform denovo calling. [string of the form \"chr:start-end\"]\n";
    std::cerr << std::endl;
    exit(1);
}
*/

/*
int mainDNG(int argc, char *argv[]) {
  std::string model;
  if(argc > 1) {
    model = argv[1];
  } else {
    usage();
  }
  std::string ped_f = "EMPTY"; // ped file
  std::string bcf_f = "EMPTY"; // bcf/vcf file
  std::string op_vcf_f = "EMPTY"; // vcf output -- optional
  bool is_vcf = false;
  std::string region = "NULL";
  // Read in Command Line arguments
    while(1) {
        int option_index = 0;
        static struct option long_options[] = {
            {"ped", 1, 0, 0},
            {"bcf", 1, 0, 1},
            {"snp_mrate", 1, 0, 2},
            {"indel_mrate", 1, 0, 3},
            {"poly_rate", 1, 0, 4},
            {"pair_mrate", 1, 0, 5},
            {"indel_mu_scale", 1, 0, 6},
            {"output_vcf", 1, 0, 7},
            {"pp_cutoff", 1, 0, 8},
            {"rd_cutoff", 1, 0, 9},
            {"h", 1, 0, 10},
            {"vcf", 1, 0, 11},
            {"region", 1, 0, 12},
            {0, 0, 0, 0}
        };
        int c = getopt_long(argc - 1, argv + 1, "", long_options, &option_index);
        if(c == -1) {
            break;
        }
        switch(c) {
        case 0:
            ped_f = optarg;
            break;
        case 1:
            bcf_f = optarg;
            break;
        case 2:
            snp_mrate = atof(optarg);
	    std::cerr << "\nSNP mutation rate: " << snp_mrate;
            break;
        case 3:
            indel_mrate = atof(optarg);
	    std::cerr << "\nindel mutation rate: " << indel_mrate;
            break;
        case 4:
            poly_rate = atof(optarg);
	    std::cerr << "\npolymorphism rate: " << poly_rate;
            break;
        case 5:
            pair_mrate = atof(optarg);
	    std::cerr << "\npaired-sample mutation rate: " << pair_mrate;
            break;
        case 6:
            mu_scale = atof(optarg);
	    std::cerr << "\nindel mutation rate scaling factor: " << mu_scale;
            break;
        case 7:
            op_vcf_f = optarg;
	    std::cerr << "\noutput vcf file: " << op_vcf_f;
            break;
        case 8:
            PP_cutoff = atof(optarg);
	    std::cerr << "\nposterior probability cutoff: " << PP_cutoff;
            break;
        case 9:
            RD_cutoff = atoi(optarg);
	    std::cerr << "\nread depth filter: " << RD_cutoff;
            break;
        case 10:
            usage();
            break;
        case 11:
            bcf_f = optarg;
            is_vcf = true;
            break;
        case 12:
            region = optarg;
            break;
        default:
            usage();
        }
    }

    if((bcf_f == "EMPTY") || (ped_f == "EMPTY")) {
      std::cerr << "ERROR ! Please specify both the PED file and BCF file !"
             "Exiting!\n";
      //usage();
    }

    // Create lookup table and read ped, BCF files. Separate for autosomes, X in males and X in females.
    if(model ==  "auto" || model == "XS" || model == "XD") {
        callDenovoFromBCF(ped_f, bcf_f, op_vcf_f, model, is_vcf, region);
    } else {
        usage();
    }

    std::cerr << "\nDone !" << std::endl;
    exit(0);
}
*/

int callMakeSNPLookup(lookup_table_t &tgtSNP, lookup_snp_t &lookupSNP, std::string model) {
  std::vector<string> tmp;
    for(int l = 0; l < 10; l++) {
        tmp.push_back("NA");
    }
    for(int l = 0; l < 100; l++) {
        tgtSNP.push_back(tmp);
    }


    if(model == "auto") { makeSNPLookup(snp_mrate, poly_rate, tgtSNP, lookupSNP); }
    else if(model == "XS") { makeXSSNPLookup(snp_mrate, poly_rate, tgtSNP, lookupSNP); }
    else if(model == "XD") { makeXDSNPLookup(snp_mrate, poly_rate, tgtSNP, lookupSNP); }
    else {
        cerr << endl << "Invalid model for SNP lookup, exiting.";
        exit(1);
    }

    cerr << "\nCreated SNP lookup table\n";
    cerr << " First mrate: " << lookupSNP.mrate(1, 1) << " last: " << lookupSNP.mrate(100, 10) << endl;
    cerr << " First code: " << lookupSNP.code(1, 1) << " last: " << lookupSNP.code(100, 10) << endl;
    cerr << " First target string: " << tgtSNP[0][0] << " last: " << tgtSNP[99][9] << endl;
    cerr << " First tref: " << lookupSNP.tref(1, 1) << " last: " << lookupSNP.tref(100, 10) << endl;
    return 0;
}

int callMakeINDELLookup(lookup_table_t &tgtIndel, lookup_indel_t &lookupIndel, std::string model) {
  std::vector<std::string> tmp;
    for(int l = 0; l < 3; l++) {
        tmp.push_back("NA");
    }
    for(int l = 0; l < 9; l++) {
        tgtIndel.push_back(tmp);
    }


    if(model == "auto") { makeIndelLookup(poly_rate, tgtIndel, lookupIndel); }
    else if(model == "XS") { makeXSIndelLookup(poly_rate, tgtIndel, lookupIndel); }
    else if(model == "XD") { makeXDIndelLookup(poly_rate, tgtIndel, lookupIndel); }
    else {
      std::cerr << std::endl << "Invalid model for INDEL lookup, exiting.";
      exit(1);
    }

    std::cerr << "\nCreated indel lookup table";
    std::cerr << " First code: " << lookupIndel.code(1, 1) << " last: " << lookupIndel.code(9, 3) << std::endl;
    std::cerr << " First target string: " << tgtIndel[0][0] << " last: " << tgtIndel[8][2] << std::endl;
    std::cerr << " First prior: " << lookupIndel.priors(1, 1) << " last: " << lookupIndel.priors(9, 3) << std::endl;
    return 0;
}

int callMakePairedLookup(lookup_table_t &tgtPair, lookup_pair_t &lookupPair) {
  std::vector<std::string> tmp;
    for(int l = 0; l < 10; l++) {
        tmp.push_back("NA");
    }
    for(int l = 0; l < 10; l++) {
        tgtPair.push_back(tmp);
    }
    makePairedLookup(pair_mrate, tgtPair, lookupPair);
    std::cerr << "\nCreated paired lookup table" << std::endl;
    std::cerr << " First target string: " << tgtPair[0][0] << " last: " << tgtPair[9][9] << std::endl;
    std::cerr << " First prior " << lookupPair.priors(1, 1) << " last: " << lookupPair.priors(10, 10) << std::endl;
    return 0;
}

void writeVCFHeader(hts::bcf::File &vcfout, std::string &bcf_file, std::string &ped_file, std::vector<std::string> &samples) {

  vcfout.AddHeaderMetadata("##fileformat=VCFv4.1");
  vcfout.AddHeaderMetadata("##source=", PACKAGE_STRING);
  vcfout.AddHeaderMetadata("##input_bcf=", bcf_file);
  vcfout.AddHeaderMetadata("##input_ped=", ped_file);
  vcfout.AddHeaderMetadata("##cutoff_read_depth=", RD_cutoff);
  vcfout.AddHeaderMetadata("##cutoff_posterior_probability=", PP_cutoff);
  vcfout.AddHeaderMetadata("##mutation_rate_snp=", snp_mrate);
  vcfout.AddHeaderMetadata("##mutation_rate_indel=", indel_mrate);
  vcfout.AddHeaderMetadata("##mutation_rate_paired=", pair_mrate);
  vcfout.AddHeaderMetadata("##mutation_rate_polymorphism=", poly_rate);
  vcfout.AddHeaderMetadata("##mutation_rate_indel_scaling_constant=", mu_scale);

#ifndef NEWVCFOUT
  // Old format, doesn't work well with htslib
  vcfout.AddHeaderMetadata("##INFO=<ID=RD_MOM,Number=1,Type=Integer,Description=\"Read depth for MOM\">");
  vcfout.AddHeaderMetadata("##INFO=<ID=RD_DAD,Number=1,Type=Integer,Description=\"Read depth for DAD\">");
  vcfout.AddHeaderMetadata("##INFO=<ID=RD_NORMAL,Number=1,Type=Integer,Description=\"Read depth for Normal sample in Paired Sample Analysis\">");
  vcfout.AddHeaderMetadata("##INFO=<ID=MQ_MOM,Number=1,Type=Integer,Description=\"Mapping quality for MOM\">");
  vcfout.AddHeaderMetadata("##INFO=<ID=MQ_DAD,Number=1,Type=Integer,Description=\"Mapping quality for DAD\">");
  vcfout.AddHeaderMetadata("##INFO=<ID=MQ_NORMAL,Number=1,Type=Integer,Description=\"Mapping quality for Normal sample in Paired Sample Analysis\">");
  vcfout.AddHeaderMetadata("##FORMAT=<ID=NULL_CONFIG(child/mom/dad),Number=1,Type=String,Description=\"NULL trio configuration\">");
  vcfout.AddHeaderMetadata("##FORMAT=<ID=PP_NULL,Number=1,Type=Float,Description=\"Posterior probability for the NULL configuration\">");
  vcfout.AddHeaderMetadata("##FORMAT=<ID=ML_NULL,Number=1,Type=Float,Description=\"Maximum Likelihood for the NULL configuration\">");
  vcfout.AddHeaderMetadata("##FORMAT=<ID=ML_DNM,Number=1,Type=Float,Description=\"Maximum Likelihood for the DNM configuration\">");
  vcfout.AddHeaderMetadata("##INFO=<ID=SNPcode,Number=1,Type=Integer,Description=\"Code that shows whether DNM configuration shown is monomorphic or contains variation\">");
  vcfout.AddHeaderMetadata("##INFO=<ID=INDELcode,Number=1,Type=Integer,Description=\" \">"); //// Needs Detail
  vcfout.AddHeaderMetadata("##INFO=<ID=PairSNPcode,Number=1,Type=Integer,Description=\" \">"); //// Needs Detail
  vcfout.AddHeaderMetadata("##INFO=<ID=code,Number=1,Type=Integer,Description=\" \">"); //// Needs Detail
  vcfout.AddHeaderMetadata("##FORMAT=<ID=DNM_CONFIG(child/mom/dad),Number=1,Type=String,Description=\"DNM trio configuration\">");
  vcfout.AddHeaderMetadata("##FORMAT=<ID=PP_DNM,Number=1,Type=Float,Description=\"Posterior probability for the DNM configuration\">");
  vcfout.AddHeaderMetadata("##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Read Depth of the child\">");
  vcfout.AddHeaderMetadata("##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality of the child\">");
  vcfout.AddHeaderMetadata("##FORMAT=<ID=RD_T,Number=1,Type=Integer,Description=\"Read Depth of the tumor sample\">");
  vcfout.AddHeaderMetadata("##FORMAT=<ID=MQ_T,Number=1,Type=Integer,Description=\"Mapping quality of the normal sample\">");
  vcfout.AddSample(samples[0].c_str());
#else
  // Newer header, produces more standarized VCF output
  vcfout.AddHeaderMetadata("##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Read Depth\">");
  vcfout.AddHeaderMetadata("##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality\">");
  for(std::string sample : samples)
    vcfout.AddSample(sample.c_str());
#endif
  vcfout.WriteHeader();
}


#endif
