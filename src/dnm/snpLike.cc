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
#include <string>
#include <fstream>
#include "parser.h"
#include "lookup.h"
#include <string.h>
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>
#include <iomanip>
#include <boost/algorithm/string.hpp>


//using namespace dng::task;
using namespace std;


// Added for testing purposes. Remove if someone think's it's ugly.
inline void printMatrix(Matrix &m) {
    std::cout << std::scientific;
    for(int a = 0; a < m.rows(); a++) {
        for(int b = 0; b < m.cols(); b++) {
            std::cout << std::setprecision(6) << m(a, b) << "  ";
        }
        std::cout << std::endl;
    }
}

// Calculate DNM and Null PP
int trio_like_snp(qcall_t &child, qcall_t &mom, qcall_t &dad, int flag,
				   lookup_table_t &tgt, lookup_snp_t &lookup,
                   hts::bcf::File *vcfout, parameters &params, Trio &trio) {

	int RD_cutoff = params.RD_cutoff;
	double pp_cutoff = params.PP_cutoff;

    // Filter low read depths ( < 10 )
    if(child.depth < RD_cutoff || mom.depth < RD_cutoff || dad.depth < RD_cutoff) {
        return 0;
    }

    Real maxlike_null, maxlike_denovo, pp_null, pp_denovo, denom;

    Matrix M(1, 10);
    Matrix C(10, 1);
    Matrix D(10, 1);
    Matrix P(10, 10);
    Matrix F(100, 10);
    Matrix L(100, 10);
    Matrix T(100, 10);
    Matrix DN(100, 10);
    Matrix PP(100, 10);
    int i, j, k, l;
    int coor = child.pos;
    char ref_name[50];
    strcpy(ref_name, child.chr); // Name of the reference sequence

    //Load Likelihood matrices L(D|Gm), L(D|Gd) and L(D|Gc)
    for(size_t a = 0; a < 10; a++) {
        M(0, a) = pow(10, -mom.lk[a] / 10.0);
    }

    for(size_t a = 0; a < 10; a++) {
        D(a, 0) = pow(10, -dad.lk[a] / 10.0);
    }

    for(size_t a = 0; a < 10; a++) {
        C(a, 0) = pow(10, -child.lk[a] / 10.0);
    }

    // Take the Kronecker products
    P = kroneckerProduct(M, D);
    F = kroneckerProduct(P, C);

    // Combine transmission probs L(Gc | Gm, Gf)
    // Take the entry wise (aka Schur) product of the two matricies
    T = F.cwiseProduct(lookup.tp);

    // Combine prior L(Gm, Gf)
    switch(mom.ref_base) {
    case 'A':
        L = T.cwiseProduct(lookup.aref);
        break;
    case 'C':
        L = T.cwiseProduct(lookup.cref);
        break;
    case 'G':
        L = T.cwiseProduct(lookup.gref);
        break;
    case 'T':
        L = T.cwiseProduct(lookup.tref);
        break;

    default: L = T; break;
    }

    DN = L.cwiseProduct(lookup.mrate);

    // Find max likelihood of null configuration
    PP = DN.cwiseProduct(lookup.norm);
    maxlike_null = PP.maxCoeff(&i, &j);

    // Find max likelihood of de novo trio configuration
    PP = DN.cwiseProduct(lookup.denovo);
    maxlike_denovo = PP.maxCoeff(&k, &l);

    denom = DN.sum();


#ifdef DEBUG_ENABLED
    cout << "\n\nMOM\n\n";
    cout << setw(10) << setprecision(10) << M;
    cout << "\n\nDAD\n\n";
    cout << setw(10) << setprecision(10) << D;
    cout << "\n\nCHILD\n\n";
    cout << setw(10) << setprecision(10) << C;
    cout << "\n\nP = L(Gm, Gd)\n\n";
    cout << setw(10) << setprecision(10) << P;
    cout << "\n\nF = L(Gc, Gm, Gd)\n\n";
    cout << setw(10) << setprecision(10) << F;
    cout << "\n\nT = L(Gc, Gm, Gd) * L(Gc | Gm, Gf)\n\n";
    cout << setw(10) << setprecision(10) << T;
    cout << "\n\nlookup.mrate\n\n";
    cout << setw(10) << setprecision(10) << lookup.mrate;
    cout << "\n\nlookup.norm\n\n";
    cout << setw(10) << setprecision(10) << lookup.norm;
    cout << "\n\nlookup.denovo\n\n";
    cout << setw(10) << setprecision(10) << lookup.denovo;
    cout << "\n\nL = L(Gc, Gm, Gd) * L(Gc | Gm, Gf) * L(Gm, Gf)\n\n";
    cout << setw(10) << setprecision(10) << L;
    cout << "\n\nDN\n\n";
    cout << setw(10) << setprecision(20) << DN;
    cout << "\n\nPP normal\n\n";
    cout << setw(10) << setprecision(20) << PP;
    cout << "\n\nPP denovo\n\n";
    cout << setw(10) << setprecision(20) << PP;
    cout << "\nmax_normal i " << i << " j " << j << " GT " << tgt[i][j];
    cout << "\nmax_denovo k " << k << " l " << l << " GT " << tgt[k][l];;
    cout << "\nDenom is " << denom;
#endif

    pp_denovo = maxlike_denovo / denom; //denovo posterior probability
    pp_null = 1 - pp_denovo; //null posterior probability

    // Check for PP cutoff
    if(pp_denovo > pp_cutoff) {

        //remove ",X" from alt, helps with VCF op.
        string alt = mom.alt;
        size_t start = alt.find(",X");
        if(start != std::string::npos) {
            alt.replace(start, 2, "");
        }

        cout << "DENOVO-SNP CHILD_ID: " << child.id;
        cout << " chr: " << ref_name;
        cout << " pos: " << coor;
        cout << " ref: " << mom.ref_base;
        cout << " alt: " << alt;
        cout << " maxlike_null: " << maxlike_null;
        cout << " pp_null: " << pp_null;
        cout << " tgt_null(child/mom/dad): " << tgt[i][j];
        cout << " snpcode: " << lookup.snpcode(i, j);
        cout << " code: " << lookup.code(i, j);
        cout << " maxlike_dnm: " << maxlike_denovo;
        cout << " pp_dnm: " << pp_denovo;
        cout << " tgt_dnm(child/mom/dad): " << tgt[k][l];
        cout << " lookup: " << lookup.code(k, l);
        cout << " flag: " << flag;
        cout << " READ_DEPTH child: " << child.depth;
        cout << " dad: " << dad.depth;
        cout << " mom: " << mom.depth;
        cout << " MAPPING_QUALITY child: " << child.rms_mapQ;
        cout << " dad: " << dad.rms_mapQ;
        cout << " mom: " << mom.rms_mapQ;
        cout << endl;

        if(vcfout != nullptr) {
        	unsigned int nsamples = vcfout->samples().second;
            auto rec = vcfout->InitVariant();
            rec.target(ref_name);
            rec.position(coor);
            rec.alleles(std::string(1, mom.ref_base) + "," + alt);
            rec.quality(0);
            rec.filter("PASS");

#ifdef NEWVCFOUT
            // Newer, more accurate VCF output format
            rec.info("SNPcode", static_cast<float>(lookup.snpcode(i, j)));
            rec.info("code", static_cast<float>(lookup.code(i, j)));
            rec.info("PP_NULL", static_cast<float>(pp_null));
            rec.info("ML_NULL", static_cast<float>(maxlike_null));
            rec.info("PP_DNM", static_cast<float>(static_cast<float>(pp_denovo)));
            rec.info("ML_DNM", static_cast<float>(static_cast<float>(maxlike_denovo)));

            std::vector<int32_t> rds(nsamples, hts::bcf::int32_missing);
            rds[trio.cpos] = child.depth;
            rds[trio.mpos] = mom.depth;
            rds[trio.dpos] = dad.depth;
            rec.samples("RD", rds);

            std::vector<int32_t> mqs(nsamples, hts::bcf::int32_missing);
            mqs[trio.cpos] = child.rms_mapQ;
            mqs[trio.mpos] = mom.rms_mapQ;
            mqs[trio.dpos] = dad.rms_mapQ;
            rec.samples("MQ", mqs);

            std::vector<std::string> configs;
            boost::split(configs, tgt[i][j], boost::is_any_of("/"));
            std::vector<std::string> null_configs(nsamples, hts::bcf::str_missing);
            null_configs[trio.cpos] = configs[0];
            null_configs[trio.mpos] = configs[1];
            null_configs[trio.dpos] = configs[2];
            rec.samples("NULL_CONFIG", null_configs);

            boost::split(configs, tgt[k][l], boost::is_any_of("/"));
            std::vector<std::string> dnm_configs(nsamples, hts::bcf::str_missing);
            dnm_configs[trio.cpos] = configs[0];
            dnm_configs[trio.mpos] = configs[1];
            dnm_configs[trio.dpos] = configs[2];
            rec.samples("DNM_CONFIG", dnm_configs);


#else
            // Old vcf output format
            rec.info("RD_MOM", mom.depth);
            rec.info("RD_DAD", dad.depth);
            rec.info("MQ_MOM", mom.rms_mapQ);
            rec.info("MQ_DAD", dad.rms_mapQ);
            rec.info("SNPcode", static_cast<float>(lookup.snpcode(i, j)));
            rec.info("code", static_cast<float>(lookup.code(i, j)));

            std::vector<std::string> null_configs(nsamples, hts::bcf::str_missing);
            null_configs[trio.cpos] = tgt[i][j];
            rec.samples("NULL_CONFIG(child/mom/dad)", null_configs);

            std::vector<float> pp_nulls(nsamples, hts::bcf::float_missing);
            pp_nulls[trio.cpos] = pp_null;
            rec.samples("PP_NULL", pp_nulls);

            std::vector<float> maxlike_nulls(nsamples, hts::bcf::float_missing);
            maxlike_nulls[trio.cpos] = maxlike_null;
            rec.samples("ML_NULL", maxlike_nulls);

            std::vector<std::string> dnm_configs(nsamples, hts::bcf::str_missing);
            dnm_configs[trio.cpos] = tgt[k][l];
            rec.samples("DNM_CONFIG(child/mom/dad)", dnm_configs);

            std::vector<float> pp_denovos(nsamples, hts::bcf::float_missing);
            pp_denovos[trio.cpos] = pp_denovo;
            rec.samples("PP_DNM", pp_denovos);

            std::vector<float> ml_dnms(nsamples, hts::bcf::float_missing);
            ml_dnms[trio.cpos] = maxlike_denovo;
            rec.samples("ML_DNM", ml_dnms);

            std::vector<int32_t> rds(nsamples, hts::bcf::int32_missing);
            rds[trio.cpos] = child.depth;
            rec.samples("RD", rds);

            std::vector<int32_t> mqs(nsamples, hts::bcf::int32_missing);
            mqs[trio.cpos] = child.rms_mapQ;
            rec.samples("MQ", mqs);
#endif
            vcfout->WriteRecord(rec);
            rec.Clear();
        }

    }

    return 1;
}
