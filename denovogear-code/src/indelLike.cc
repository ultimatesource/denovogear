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

#define MIN_READ_DEPTH_INDEL 10

// parameters for the indel mutation rate model
const float INSERTION_SLOPE = -0.2994;
const float INSERTION_INTERCEPT = -22.8689;
const float DELETION_SLOPE = -0.2856;
const float DELETION_INTERCEPT = -21.9313;



using namespace std;

// Calculate DNM and Null PP
void trio_like_indel(indel_t *child,indel_t *mom, indel_t *dad, int flag, 
                     vector<vector<string > > & tgtIndel, 
                     lookup_indel_t & lookupIndel, double mu_scale, 
                     string op_vcf_f, ofstream& fo_vcf, double pp_cutoff, 
                     int RD_cutoff, int& n_site_pass)
{  
  // Read depth filter
  if (child->depth < RD_cutoff ||
    mom->depth < RD_cutoff || dad->depth < RD_cutoff) {
    return;
  }

  n_site_pass += 1;
  Real a[3];   
  Real maxlike_null,maxlike_denovo,pp_null,pp_denovo,denom;       
  Matrix M(1,3);
  Matrix C(3,1);
  Matrix D(3,1);
  Matrix P(3,3);
  Matrix F(9,3);
  Matrix L(9,3);
  Matrix T(9,3);
  Matrix DN(9,3);
  Matrix PP(9,3);
  int i,j,k,l;
  int coor = child->pos;
  char ref_name[50];
  strcpy( ref_name, child->chr); // Name of the reference sequence

  // this is the case where child has more than 1 indel alt allele, ignore for now
  if (strstr(child->alt, ",") != 0) {
  	return;
  }

	bool is_insertion = false; // insertion/deletion event
	int len_diff = strlen(mom->ref_base) - strlen(mom->alt); // diff b/w alt, ref
	if (len_diff < 0) {
		is_insertion = true;
		len_diff = -len_diff;
	} else if (len_diff > 0) {
		is_insertion = false;
	}

  // mu_scale is the variable used to scale the mutation rate prior
  double log_indel_mrate, indel_mrate, new_indel_mrate;
  if(is_insertion) 
    log_indel_mrate = mu_scale * (INSERTION_SLOPE*len_diff + INSERTION_INTERCEPT); 
  else
    log_indel_mrate = mu_scale * (DELETION_SLOPE*len_diff + DELETION_INTERCEPT); 
	indel_mrate = exp(log_indel_mrate); // antilog of log ratio
	
	for (int j = 1; j <= 9; j++) {
		for (int l = 1; l <= 3; l++) {
			new_indel_mrate = pow(indel_mrate, lookupIndel.hit(j,l)); // hit is 0,1 or 2 (number of indels in the trio config)
			lookupIndel.mrate(j,l) = new_indel_mrate; 
		}
  } 

  //Load likelihood vectors    
  for (j = 0; j != 3; ++j) a[j]=pow(10,-mom->lk[j]/10.);
    M<<a;
  for (j = 0; j != 3; ++j) a[j]=pow(10,-dad->lk[j]/10.);
    D<<a;
  for (j = 0; j != 3; ++j) a[j]=pow(10,-child->lk[j]/10.);
    C<<a;

  P = KP(M,D);    
  F = KP(P,C);    
  // combine with transmission probs 
  T = SP(F,lookupIndel.tp); 
  // combine with priors
  L = SP(T, lookupIndel.priors);
  // combine with mutation rate
  DN = SP(L,lookupIndel.mrate);    
  
  // Find max likelihood of null configuration
  PP = SP(DN,lookupIndel.norm);   //zeroes out configurations with mendelian error
  maxlike_null = PP.maximum2(i,j);       

  //Find max likelihood of de novo trio configuration
  PP = SP(DN,lookupIndel.denovo);   //zeroes out configurations with mendelian inheritance
  maxlike_denovo = PP.maximum2(k,l); 

  //make proper posterior probs
  denom = DN.sum();
  pp_denovo = maxlike_denovo/denom; // denovo posterior probability
  pp_null = 1 - pp_denovo; // null posterior probability

  // Check for PP cutoff
  if ( pp_denovo > pp_cutoff ) {

    //remove ",X" from alt, helps with VCF op.
    string alt = mom->alt;  
    size_t start = alt.find(",X");
    if(start != std::string::npos)
        alt.replace(start, 2, "");


  	cout<<"DENOVO-INDEL child id: "<<child->id;
  	cout<<" chr: "<<ref_name<<" pos: "<<coor<<" ref: "<<mom->ref_base<<" alt: "<<alt;
  	cout<<" maxlike_null: "<<maxlike_null<<" pp_null: "<<pp_null<<" tgt: "<<tgtIndel[i-1][j-1];
  	cout<<" snpcode: "<<lookupIndel.snpcode(i,j)<<" code: "<<lookupIndel.code(i,j);
  	cout<<" maxlike_dnm: "<<maxlike_denovo<<" pp_dnm: "<<pp_denovo;
  	cout<<" tgt: "<<tgtIndel[k-1][l-1]<<" lookup: "<<lookupIndel.code(k,l)<<" flag: "<<flag;
  	cout<<" READ_DEPTH child: "<<child->depth<<" dad: "<<dad->depth<<" mom: "<<mom->depth;
    cout<<" MAPPING_QUALITY child: "<<child->rms_mapQ<<" dad: "<<dad->rms_mapQ<<" mom: "<<mom->rms_mapQ;
    cout<<endl;

    if(op_vcf_f != "EMPTY") {
      fo_vcf<<ref_name<<"\t";
      fo_vcf<<coor<<"\t";
      fo_vcf<<".\t";// Don't know the rsID
      fo_vcf<<mom->ref_base<<"\t";
      fo_vcf<<alt<<"\t";
      fo_vcf<<"0\t";// Quality of the Call
      fo_vcf<<"PASS\t";// passed the read depth filter
      fo_vcf<<"RD_MOM="<<mom->depth<<";RD_DAD="<<dad->depth;
      fo_vcf<<";MQ_MOM="<<mom->rms_mapQ<<";MQ_DAD="<<dad->rms_mapQ;  
      fo_vcf<<";INDELcode="<<lookupIndel.snpcode(i,j)<<";\t";
      fo_vcf<<"NULL_CONFIG(child/mom/dad):PP_NULL:DNM_CONFIG(child/mom/dad):PP_DNM:RD:MQ\t";
      fo_vcf<<tgtIndel[i-1][j-1]<<":"<<pp_null<<":"; 
      fo_vcf<<tgtIndel[k-1][l-1]<<":"<<pp_denovo<<":"<<child->depth<<":"<<child->rms_mapQ; 
      fo_vcf<<"\n";
    }
  }    
}
