#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include "parser.h"
#include "lookup.h"

#ifndef NEWMAT_H
#define NEWMAT_H
#include "newmatap.h"
#include "newmatio.h"
#endif

using namespace std;

// Calculate Pair PP
void pair_like(pair_t tumor, pair_t normal, vector<vector<string> > &tgtPair, 
	       lookup_pair_t & lookupPair, int flag, string op_vcf_f, ofstream& fo_vcf, 
         double pp_cutoff, int RD_cutoff, int& n_site_pass)
{
  // Filter low read depths
  if (tumor.depth < RD_cutoff || normal.depth < RD_cutoff) {
    return;
  }
  n_site_pass += 1;
  Real a[10];   
  Real maxlike_null, maxlike_denovo, pp_null, pp_denovo, denom;   
  Matrix N(1,10);
  Matrix T(10,1);
  Matrix P(10,10);
  Matrix DN(10,10);
  Matrix PP(100,10);
  int i,j,k,l;  
  int coor = tumor.pos;
  char ref_name[50];
  strcpy( ref_name, tumor.chr); // Name of the reference sequence
    
  //Load Likelihood matrices L(D|Gt) and L(D|Gn) 
  for (j = 0; j != 10; ++j)	
  	a[j]=pow(10, -normal.lk[j]/10.); 
  N<<a;
  
  for (j = 0; j != 10; ++j) 
  	a[j] = pow(10, -tumor.lk[j]/10.);
  T<<a;

  P = KP(N, T); // 10 * 10
  // Combine transmission probs L(Gc | Gm, Gf)
  DN = SP(P, lookupPair.priors); // 10 * 10


  // Find max likelihood of null configuration
  PP = SP(DN, lookupPair.norm);   //zeroes out configurations with mendelian error
  maxlike_null = PP.maximum2(i,j);   
  
  // Find max likelihood of de novo trio configuration
  PP = SP(DN, lookupPair.denovo);   //zeroes out configurations with mendelian inheritance
  maxlike_denovo = PP.maximum2(k,l); 

  denom = DN.sum();

  pp_denovo = maxlike_denovo / denom; // denovo posterior probability
  pp_null = 1 - pp_denovo; // null posterior probability

  // Check for PP cutoff 
  if ( pp_denovo > pp_cutoff ) {
    cout<<"DENOVO-PAIR TUMOR ID: "<<tumor.id<<" NORMAL ID: "<<normal.id;
    cout<<" ref_name: "<<ref_name<<" coor: "<<coor<<" ref_base: "<<tumor.ref_base<<" ALT: "<<tumor.alt;
    cout<<" maxlike_null: "<<maxlike_null<<" pp_null: "<<pp_null<<" tgt_null(normal/tumor): "<<tgtPair[i-1][j-1];
    cout<<" maxlike_dnm: "<<maxlike_denovo<<" pp_dnm: "<<pp_denovo;
    cout<<" tgt_denovo(normal/tumor): "<<tgtPair[k-1][l-1];//<<" flag: "<<flag; // flag is a variable that could be set in denovogear.cc(site specific info)
    cout<<" READ_DEPTH tumor: "<<tumor.depth<<" normal: "<<normal.depth;
    cout<<" MAPPING_QUALITY tumor: "<<tumor.rms_mapQ<<" normal: "<<normal.rms_mapQ;
    cout<<" null_snpcode: "<<lookupPair.snpcode(i, j);
    cout<<" dnm_snpcode: "<<lookupPair.snpcode(k, l);
    cout<<endl;

    if(op_vcf_f != "EMPTY") {
      fo_vcf<<ref_name<<"\t";
      fo_vcf<<coor<<"\t";
      fo_vcf<<".\t";// Don't know the rsID
      fo_vcf<<tumor.ref_base<<"\t";
      fo_vcf<<tumor.alt<<"\t";
      fo_vcf<<"0\t";// Quality of the Call
      fo_vcf<<"PASS\t";// passed the read depth filter
      fo_vcf<<"RD_NORMAL="<<normal.depth<<";";
      fo_vcf<<";MQ_NORMAL="<<normal.rms_mapQ<<";";  
      fo_vcf<<";NULL_CONFIG(normal/tumor)="<<tgtPair[i-1][j-1]<<";PP_NULL="<<pp_null;  
      fo_vcf<<";PairSNPcode="<<lookupPair.snpcode(k, l)<<";\t";
      fo_vcf<<"DNM_CONFIG(normal/tumor):PP_DNM:RD_T:MQ_T\t";
      fo_vcf<<tgtPair[k-1][l-1]<<":"<<pp_denovo<<":"<<tumor.depth<<":"<<tumor.rms_mapQ; 
      fo_vcf<<"\n";
    }    
  }
}
