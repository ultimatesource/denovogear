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
#include "denovogear.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include "parser.h"
#include "lookup.h"
//#include "newmatap.h"
//#include "newmatio.h"
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>
#include <boost/algorithm/string.hpp>

using namespace std;

//typedef double Real;
//typedef Eigen::MatrixXd Matrix;

// Calculate Pair PP
void pair_like(pair_t &tumor, pair_t &normal, vector<vector<string> > &tgtPair,
               lookup_pair_t &lookupPair, int flag, /*vector<hts::bcf::File> &vcfout*/ hts::bcf::File *vcfout,
               parameters &params, int &n_site_pass, Pair &pair) {

	int RD_cutoff = params.RD_cutoff;
	double pp_cutoff = params.PP_cutoff;


    // Filter low read depths
    if(tumor.depth < RD_cutoff || normal.depth < RD_cutoff) {
        return;
    }
    n_site_pass += 1;
    //Real a[10];
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
        N(0, j) = pow(10, -normal.lk[j] / 10.);
        //a[j] = pow(10, -normal.lk[j] / 10.);
    }
    //N << a;

    for(j = 0; j != 10; ++j) {
        T(j, 0) = pow(10, -tumor.lk[j] / 10.);
        //a[j] = pow(10, -tumor.lk[j] / 10.);
    }
    //T << a;

    P = kroneckerProduct(N, T);
    //P = KP(N, T); // 10 * 10
    // Combine transmission probs L(Gc | Gm, Gf)
    DN = P.cwiseProduct(lookupPair.priors);
    //DN = SP(P, lookupPair.priors); // 10 * 10


    // Find max likelihood of null configuration
    PP = DN.cwiseProduct(lookupPair.norm);
    maxlike_null = PP.maxCoeff(&i, &j);
    //PP = SP(DN, lookupPair.norm);   //zeroes out configurations with mendelian error
    //maxlike_null = PP.maximum2(i, j);

    // Find max likelihood of de novo trio configuration
    PP = DN.cwiseProduct(lookupPair.denovo);
    maxlike_denovo = PP.maxCoeff(&k, &l);
    //PP = SP(DN, lookupPair.denovo);   //zeroes out configurations with mendelian inheritance
    //maxlike_denovo = PP.maximum2(k, l);

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

        cout << "DENOVO-PAIR-SNP TUMOR_ID: " << tumor.id;
	cout << " NORMAL_ID: " << normal.id;
        cout << " chr: " << ref_name;
	cout << " pos: " << coor;
	cout << " ref: " << tumor.ref_base;
	cout << " alt: " << alt;
        cout << " maxlike_null: " << maxlike_null;
	cout << " pp_null: " << pp_null;
	cout << " tgt_null(normal/tumor): " << tgtPair[i][j];
        cout << " maxlike_dnm: " << maxlike_denovo;
	cout << " pp_dnm: " << pp_denovo;
        cout << " tgt_dnm(normal/tumor): " << tgtPair[k][l]; 
        cout << " READ_DEPTH tumor: " << tumor.depth;
	cout << " normal: " << normal.depth;
        cout << " MAPPING_QUALITY tumor: " << tumor.rms_mapQ;
	cout << " normal: " << normal.rms_mapQ;
        cout << " null_snpcode: " << lookupPair.snpcode(i, j);
        cout << " dnm_snpcode: " << lookupPair.snpcode(k, l);
        cout << endl;

        if(vcfout != nullptr) {
            auto rec = vcfout->InitVariant();
            unsigned int nsamples = vcfout->samples().second;
            rec.target(ref_name);
            rec.position(coor);
            rec.alleles(std::string(1, tumor.ref_base) + "," + alt);
            rec.quality(0);
            rec.filter("PASS");
#ifndef NEWVCFOUT
            // Old vcf output format
            rec.info("RD_NORMAL", normal.depth);
            rec.info("MQ_NORMAL", normal.rms_mapQ);

            std::vector<std::string> null_configs(nsamples, hts::bcf::str_missing);
            null_configs[pair.tpos] = tgtPair[i][j];
            rec.samples("NULL_CONFIG(normal/tumor)", null_configs);

            std::vector<float> pair_null_codes(nsamples, hts::bcf::float_missing);
            pair_null_codes[pair.tpos] = lookupPair.snpcode(i, j);
            rec.samples("pair_null_code", pair_null_codes);

            std::vector<float> pair_denovo_codes(nsamples, hts::bcf::float_missing);
            pair_denovo_codes[pair.tpos] = lookupPair.snpcode(k, l);
            rec.samples("pair_denovo_code", std::vector<float> {static_cast<float>(lookupPair.snpcode(k, l))});

            std::vector<float> pp_nulls(nsamples, hts::bcf::float_missing);
            pp_nulls[pair.tpos] = pp_null;
            rec.samples("PP_NULL", pp_nulls);

            std::vector<std::string> dnm_configs(nsamples, hts::bcf::str_missing);
            dnm_configs[pair.tpos] = tgtPair[k][l];
            rec.samples("DNM_CONFIG(tumor/normal)", dnm_configs);

            std::vector<float> pp_denovos(nsamples, hts::bcf::float_missing);
            pp_denovos[pair.tpos] = pp_denovo;
            rec.samples("PP_DNM", pp_denovos);

            std::vector<int32_t> rds(nsamples, hts::bcf::int32_missing);
            rds[pair.tpos] = tumor.depth;
            rec.samples("RD_T", rds);

            std::vector<int32_t> mqs(nsamples, hts::bcf::int32_missing);
            mqs[pair.tpos] = tumor.rms_mapQ;
            rec.samples("MQ_T", mqs);


#else

            // Newer, more accurate VCF output format
            rec.info("pair_null_code", static_cast<float>(lookupPair.snpcode(i, j)));
            rec.info("pair_denovo_code", static_cast<float>(lookupPair.snpcode(k, l)));
            rec.info("PP_NULL", static_cast<float>(pp_null));
            rec.info("PP_DNM", static_cast<float>(pp_denovo));

            std::vector<int32_t> rds(nsamples, hts::bcf::int32_missing);
            rds[pair.npos] = normal.depth;
            rds[pair.tpos] = tumor.depth;
            rec.samples("RD", rds);

            std::vector<int32_t> mqs(nsamples, hts::bcf::int32_missing);
            mqs[pair.npos] = normal.rms_mapQ;
            mqs[pair.tpos] = tumor.rms_mapQ;
            rec.samples("MQ", mqs);

            std::vector<std::string> configs;
            boost::split(configs, tgtPair[i][j], boost::is_any_of("/"));
            std::vector<std::string> null_configs(nsamples, hts::bcf::str_missing);
            null_configs[pair.npos] = configs[0];
            null_configs[pair.tpos] = configs[1];
            rec.samples("NULL_CONFIG", null_configs);

            boost::split(configs, tgtPair[k][l], boost::is_any_of("/"));
            std::vector<std::string> dnm_configs(nsamples, hts::bcf::str_missing);
            dnm_configs[pair.npos] = configs[0];
            dnm_configs[pair.tpos] = configs[1];
            rec.samples("DNM_CONFIG", dnm_configs);

#endif
            vcfout->WriteRecord(rec);
            rec.Clear();
        }
    }
}
