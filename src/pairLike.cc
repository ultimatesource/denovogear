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
#include <iostream>
#include <fstream>
#include <string.h>
#include "parser.h"
#include "lookup.h"

#include "newmatap.h"
#include "newmatio.h"


using namespace std;

// Calculate Pair PP
void pair_like(pair_t tumor, pair_t normal, vector<vector<string> > &tgtPair,
               lookup_pair_t &lookupPair, int flag, vector<hts::bcf::File> &vcfout,
               double pp_cutoff, int RD_cutoff, int &n_site_pass) {
    // Filter low read depths
    if(tumor.depth < RD_cutoff || normal.depth < RD_cutoff) {
        return;
    }
    n_site_pass += 1;
    Real a[10];
    Real maxlike_null, maxlike_denovo, pp_null, pp_denovo, denom;
    Matrix N(1, 10);
    Matrix T(10, 1);
    Matrix P(10, 10);
    Matrix DN(10, 10);
    Matrix PP(100, 10);
    int i, j, k, l;
    int coor = tumor.pos;
    char ref_name[50];
    strcpy(ref_name, tumor.chr);  // Name of the reference sequence

    //Load Likelihood matrices L(D|Gt) and L(D|Gn)
    for(j = 0; j != 10; ++j) {
        a[j] = pow(10, -normal.lk[j] / 10.);
    }
    N << a;

    for(j = 0; j != 10; ++j) {
        a[j] = pow(10, -tumor.lk[j] / 10.);
    }
    T << a;

    P = KP(N, T); // 10 * 10
    // Combine transmission probs L(Gc | Gm, Gf)
    DN = SP(P, lookupPair.priors); // 10 * 10


    // Find max likelihood of null configuration
    PP = SP(DN, lookupPair.norm);   //zeroes out configurations with mendelian error
    maxlike_null = PP.maximum2(i, j);

    // Find max likelihood of de novo trio configuration
    PP = SP(DN,
            lookupPair.denovo);   //zeroes out configurations with mendelian inheritance
    maxlike_denovo = PP.maximum2(k, l);

    denom = DN.sum();

    pp_denovo = maxlike_denovo / denom; // denovo posterior probability
    pp_null = 1 - pp_denovo; // null posterior probability

    // Check for PP cutoff
    if(pp_denovo > pp_cutoff) {

        //remove ",X" from alt, helps with VCF op.
        string alt = tumor.alt;
	formatAltAlleles(alt);
        //size_t start = alt.find(",X");
        //if(start != std::string::npos) {
        //    alt.replace(start, 2, "");
        //}

        cout << "DENOVO-PAIR-SNP TUMOR_ID: " << tumor.id << " NORMAL_ID: " << normal.id;
        cout << " chr: " << ref_name << " pos: " << coor << " ref: " << tumor.ref_base
             << " alt: " << alt;
        cout << " maxlike_null: " << maxlike_null << " pp_null: " << pp_null <<
             " tgt_null(normal/tumor): " << tgtPair[i - 1][j - 1];
        cout << " maxlike_dnm: " << maxlike_denovo << " pp_dnm: " << pp_denovo;
        cout << " tgt_dnm(normal/tumor): " << tgtPair[k - 1][l -
                1]; //<<" flag: "<<flag; // flag is a variable that could be set in denovogear.cc(site specific info)
        cout << " READ_DEPTH tumor: " << tumor.depth << " normal: " << normal.depth;
        cout << " MAPPING_QUALITY tumor: " << tumor.rms_mapQ << " normal: " <<
             normal.rms_mapQ;
        cout << " null_snpcode: " << lookupPair.snpcode(i, j);
        cout << " dnm_snpcode: " << lookupPair.snpcode(k, l);
        cout << endl;

	if(!vcfout.empty()) {
	  auto rec = vcfout[0].InitVariant();
	  rec.target(ref_name);
	  rec.position(coor);
	  rec.alleles(std::string(1, tumor.ref_base) + "," + alt); 
	  rec.quality(0);
	  rec.filter("PASS");
#ifndef NEWVCFOUT
	  // Old vcf output format
	  rec.info("RD_NORMAL", normal.depth);
	  rec.info("MQ_NORMAL", normal.rms_mapQ);
	  rec.samples("NULL_CONFIG(normal/tumor)", std::vector<std::string>{tgtPair[i - 1][j - 1]});
	  rec.samples("pair_null_code", std::vector<float>{static_cast<float>(lookupPair.snpcode(i,j))});
	  rec.samples("PP_NULL", std::vector<float>{static_cast<float>(pp_null)});
	  rec.samples("DNM_CONFIG(tumor/normal)", std::vector<std::string>{tgtPair[k - 1][l - 1]});
	  rec.samples("pair_denovo_code", std::vector<float>{static_cast<float>(lookupPair.snpcode(k,l))});
	  rec.samples("PP_DNM", std::vector<float>{static_cast<float>(pp_denovo)});
	  rec.samples("RD_T", std::vector<int32_t>{tumor.depth});
	  rec.samples("MQ_T", std::vector<int32_t>{tumor.rms_mapQ});
#else
	  
	  // Newer, more accurate VCF output format
	  rec.info("pair_null_code", static_cast<float>(lookupPair.snpcode(i,j)));
	  rec.info("pair_denovo_code", static_cast<float>(lookupPair.snpcode(k,l)));
	  rec.info("PP_NULL", static_cast<float>(pp_null));
	  rec.info("PP_DNM", static_cast<float>(pp_denovo));
	  rec.samples("RD", std::vector<int32_t>{tumor.depth, normal.depth});
	  rec.samples("MQ", std::vector<int32_t>{tumor.rms_mapQ, normal.rms_mapQ});
	  std::vector<std::string> configs;
	  boost::split(configs, tgt[i-1][j-1], boost::is_any_of("/"));
	  rec.samples("NULL_CONFIG", configs);
	  boost::split(configs, tgt[k-1][l-1], boost::is_any_of("/"));
	  rec.samples("DNM_CONFIG", configs);  
#endif
	  vcfout[0].WriteRecord(rec);
	}

	/*
        if(op_vcf_f != "EMPTY") {
            fo_vcf << ref_name << "\t";
            fo_vcf << coor << "\t";
            fo_vcf << ".\t"; // Don't know the rsID
            fo_vcf << tumor.ref_base << "\t";
            fo_vcf << alt << "\t";
            fo_vcf << "0\t"; // Quality of the Call
            fo_vcf << "PASS\t"; // passed the read depth filter
            fo_vcf << "RD_NORMAL=" << normal.depth;
            fo_vcf << ";MQ_NORMAL=" << normal.rms_mapQ << ";\t";
            fo_vcf << "NULL_CONFIG(normal/tumor):pair_null_code:PP_NULL:DNM_CONFIG(normal/tumor):pair_denovo_code:PP_DNM:RD_T:MQ_T\t";
            fo_vcf << tgtPair[i - 1][j - 1] << ":" << lookupPair.snpcode(i,
                    j) << ":" << pp_null << ":";
            fo_vcf << tgtPair[k - 1][l - 1] << ":" << lookupPair.snpcode(k,
                    l) << ":" << pp_denovo << ":" << tumor.depth << ":" << tumor.rms_mapQ;
            fo_vcf << "\n";
        }
	*/
    }
}
