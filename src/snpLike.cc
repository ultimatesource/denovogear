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
#include <string>
#include <fstream>
#include "parser.h"
#include "lookup.h"
#include <string.h>
#include <Eigen/KroneckerProduct>
#include <iomanip>

using namespace std;

// TODO: Move all the definitions of Real into a single header
typedef double Real;
typedef Eigen::MatrixXd Matrix;

void printMatrix(Matrix &m) {
  std::cout << std::scientific;
  for(int a = 0; a < m.rows(); a++) {
    for(int b = 0; b < m.cols(); b++) {
      std::cout << std::setprecision(6) << m(a,b) << "  ";
    }
    std::cout << std::endl;
  }
}

// Calculate DNM and Null PP
// TODO: tgt object is already typedefed
void trio_like_snp(qcall_t child, qcall_t mom, qcall_t dad, int flag, vector<vector<string > > &tgt, lookup_snp_t &lookup, 
		   std::vector<hts::bcf::File> &vcfout, double pp_cutoff, int RD_cutoff, int &n_site_pass) {

    // Filter low read depths ( < 10 )
    if(child.depth < RD_cutoff || mom.depth < RD_cutoff || dad.depth < RD_cutoff) {
        return;
    }
    n_site_pass += 1;
    //Real a[10];
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
    //std::cout << "M =";
    for(size_t a = 0; a < 10; a++) {
      M(0,a) = pow(10, -mom.lk[a] / 10.0);
      //a[j] = pow(10, -mom.lk[a] / 10.0);
      //std::cout << " " << a[j]
    }
    //M << a;

    
    for(size_t a = 0; a < 10; a++) {
      D(a,0) = pow(10, -dad.lk[a] / 10.0);
      //a[j] = pow(10, -dad.lk[j] / 10.0);
    }
    //D << a;

    for(size_t a = 0; a < 10; a++) {
      C(a,0) = pow(10, -child.lk[a] / 10.0);
      //a[j] = pow(10, -child.lk[j] / 10.0);
    }
    //C << a;

    // Take the Kronecker products
    P = kroneckerProduct(M, D);
    F = kroneckerProduct(P, C);
    //P = KP(M, D); // 10 * 10
    //F = KP(P, C); // 100 * 10

    // Combine transmission probs L(Gc | Gm, Gf)
    // Take the entry wise (aka Schur) product of the two matricies
    T = F.cwiseProduct(lookup.tp);
    //T = SP(F, lookup.tp); //100 * 10

    // Combine prior L(Gm, Gf)
    switch(mom.ref_base) {
    case 'A':
      L = T.cwiseProduct(lookup.aref);
      //L = SP(T, lookup.aref); 
      break;
    case 'C':
      L = T.cwiseProduct(lookup.cref);
      //L = SP(T, lookup.cref); 
      break;
    case 'G':
      L = T.cwiseProduct(lookup.gref);
      //L = SP(T, lookup.gref); 
      break;
    case 'T':
      L = T.cwiseProduct(lookup.tref);
      //L = SP(T, lookup.tref); 
      break;
      
    default: L = T; break;
    }

    //std::cout << "L = " << std::endl << L << std::endl;

    DN = L.cwiseProduct(lookup.mrate);
    //DN = SP(L, lookup.mrate);

    // Find max likelihood of null configuration
    
    PP = DN.cwiseProduct(lookup.norm);
    //PP = SP(DN, lookup.norm); //zeroes out configurations with mendelian error
    maxlike_null = PP.maxCoeff(&i, &j);
    
    //maxlike_null = PP.maximum2(i, j);


    // Find max likelihood of de novo trio configuration
    PP = DN.cwiseProduct(lookup.denovo);
    //PP = SP(DN, lookup.denovo); //zeroes out configurations with mendelian inheritance
    maxlike_denovo = PP.maxCoeff(&k, &l);
    //maxlike_denovo = PP.maximum2(k, l);

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
    
    /*
    std::cout << "ref base = " << mom.ref_base << std::endl;
    //std::cout << "M = " << M << std::endl;
    //std::cout << "P ----- " << std::endl; printMatrix(P);  std::cout << "------- " << std::endl;
    //std::cout << "F ----- " << std::endl; printMatrix(F); std::cout << "------- " << std::endl; //
    //std::cout << "lookup.tp ----- " << std::endl; printMatrix(lookup.tp); std::cout << "------- " << std::endl;
    //std::cout << "lookup.norm ----- " << std::endl; printMatrix(lookup.norm); std::cout << "------- " << std::endl; //    
    //std::cout << "T ----- " << std::endl; printMatrix(T); std::cout << "------- " << std::endl; //
    //std::cout << "lookup.aref ----- " << std::endl; printMatrix(lookup.aref); std::cout << "------- " << std::endl; //
    std::cout << "PP ----- " << std::endl; printMatrix(PP); std::cout << "------- " << std::endl; //
   
    std::cout << "F(" << k << "," << l << ") = " << F(k,l) << std::endl;
    std::cout << "L(" << k << "," << l << ") = " << L(k,l) << std::endl;
    std::cout << "DN(" << k << "," << l << ") = " << DN(k,l) << std::endl;
    std::cout << "maxlike_denovo " << maxlike_denovo << "(" << k << "," << l <<")" << std::endl;
    std::cout << "PP_[" << PP.rows() << "x" << PP.cols() << "]" << "\t(" << k << "," << l << ") = " << PP(k,l) << std::endl;
    std::cout << "pp_denovo = " << pp_denovo << " (" << maxlike_denovo << "/" << denom << ") " << std::endl;
    std::cout << "tgt[" << i << "," << j << "] = " << tgt[i][j] << std::endl;
    std::cout << "tgt[" << k << "," << l << "] = " << tgt[k][l] << std::endl;
    std::cout << std::endl;
    */
    if(pp_denovo > pp_cutoff) {
      
        //remove ",X" from alt, helps with VCF op.
        string alt = mom.alt;
        size_t start = alt.find(",X");
        if(start != std::string::npos) {
            alt.replace(start, 2, "");
        }

        cout << "DENOVO-SNP CHILD_ID: " << child.id;
        cout << " chr: " << ref_name << " pos: " << coor << " ref: " << mom.ref_base <<
             " alt: " << alt;
        cout << " maxlike_null: " << maxlike_null << " pp_null: " << pp_null <<
             " tgt_null(child/mom/dad): " << tgt[i][j];
        cout << " snpcode: " << lookup.snpcode(i, j) << " code: " << lookup.code(i, j);
        cout << " maxlike_dnm: " << maxlike_denovo << " pp_dnm: " << pp_denovo;
        cout << " tgt_dnm(child/mom/dad): " << tgt[k][l] << " lookup: " <<
             lookup.code(k, l) << " flag: " << flag;
        cout << " READ_DEPTH child: " << child.depth << " dad: " << dad.depth <<
             " mom: " << mom.depth;
        cout << " MAPPING_QUALITY child: " << child.rms_mapQ << " dad: " << dad.rms_mapQ
             << " mom: " << mom.rms_mapQ;
        cout << endl;

	if(!vcfout.empty()) {
	  auto rec = vcfout[0].InitVariant();
	  rec.target(ref_name);
	  rec.position(coor);
	  rec.alleles(std::string(1, mom.ref_base) + "," + alt); 
	  rec.quality(0);
	  rec.filter("PASS");
#ifndef NEWVCFOUT
	  // Old vcf output format
	  rec.info("RD_MOM", mom.depth);
	  rec.info("RD_DAD", dad.depth);
	  rec.info("MQ_MOM", mom.rms_mapQ);
	  rec.info("MQ_DAD", dad.rms_mapQ);
	  rec.info("SNPcode", static_cast<float>(lookup.snpcode(i,j)));
	  rec.info("code", static_cast<float>(lookup.code(i,j)));
	  rec.samples("NULL_CONFIG(child/mom/dad)", std::vector<std::string>{tgt[i][j]});
	  rec.samples("PP_NULL", std::vector<float>{static_cast<float>(pp_null)});
	  rec.samples("ML_NULL", std::vector<float>{static_cast<float>(maxlike_null)});
	  rec.samples("DNM_CONFIG(child/mom/dad)", std::vector<std::string>{tgt[k][l]});
	  rec.samples("PP_DNM", std::vector<float>{static_cast<float>(pp_denovo)});
	  rec.samples("ML_DNM", std::vector<float>{static_cast<float>(maxlike_denovo)});
	  rec.samples("RD", std::vector<int32_t>{child.depth});
	  rec.samples("MQ", std::vector<int32_t>{child.rms_mapQ});
#else
	  
	  // Newer, more accurate VCF output format
	  rec.info("SNPcode", static_cast<float>(lookup.snpcode(i,j)));
	  rec.info("code", static_cast<float>(lookup.code(i,j)));
	  rec.info("PP_NULL", static_cast<float>(pp_null));
	  rec.info("ML_NULL", static_cast<float>(maxlike_null));
	  rec.info("PP_DNM", static_cast<float>(static_cast<float>(pp_denovo)));
	  rec.info("ML_DNM", static_cast<float>(static_cast<float>(maxlike_denovo)));
	  rec.samples("RD", std::vector<int32_t>{child.depth, mom.depth, dad.depth});
	  rec.samples("MQ", std::vector<int32_t>{child.rms_mapQ, mom.rms_mapQ, dad.rms_mapQ});
	  std::vector<std::string> configs;
	  boost::split(configs, tgt[i][j], boost::is_any_of("/"));
	  rec.samples("NULL_CONFIG", configs);
	  boost::split(configs, tgt[k][l], boost::is_any_of("/"));
	  rec.samples("DNM_CONFIG", configs);  
#endif
	  vcfout[0].WriteRecord(rec);
	  }
    
     }

}
