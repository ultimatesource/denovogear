#include <vector>
#include <iostream>
#include <string.h>
#include "parser.h"
#include "lookup.h"

#ifndef NEWMAT_H
#define NEWMAT_H
#include "newmatap.h"
#include "newmatio.h"
#endif

#define MIN_READ_DEPTH_INDEL 10

using namespace std;

// Calculate Pair PP
void pair_like(pair_t tumor, pair_t normal, vector<vector<string> > &tgtPair, 
	       lookup_pair_t & lookupPair, int flag)
{
  Real a[10];   
  Real maxlike_null, maxlike_denovo, pp_null, pp_denovo, denom;   
  Matrix N(1,10);
  Matrix T(10,1);
  Matrix P(10,10);
  Matrix DN(10,10);
  Matrix PP(100,10);
  int i,j,k,l;  
  int coor = tumor.pos;
  char ref_name[3];
  strcpy( ref_name, tumor.chr); // Name of the reference sequence
  
  // Filter low read depths
  if (tumor.depth < MIN_READ_DEPTH_SNP || normal.depth < MIN_READ_DEPTH_SNP) {
    return;
  }
  
  //Load Likelihood matrices L(D|Gt) and L(D|Gn) 
  for (j = 0; j != 10; ++j)	{ 
  	a[j]=pow(10,-normal.lk[j]/10.); 
  }
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
  if ( pp_denovo > 0.1 ) {
    cout<<"\nDENOVO-PAIR TUMOR ID: "<<tumor.id;
    cout<<" ref_name: "<<ref_name<<" coor: "<<coor<<" ref_base: "<<tumor.ref_base<<" ALT: "<<tumor.alt;
    cout<<" maxlike_null: "<<maxlike_null<<" pp_null: "<<pp_null<<" tgt: "<<tgtPair[i-1][j-1];
    cout<<" snpcode: "<<lookupPair.snpcode(i,j);
    cout<<" maxlike_dnm: "<<maxlike_denovo<<" pp_dnm: "<<pp_denovo;
    cout<<" tgt: "<<tgtPair[k-1][l-1]<<" flag: "<<flag;
    printf(" READ_DEPTH tumor: %d normal: %d", tumor.depth, normal.depth);
    printf(" MAPPING_QUALITY tumor: %d normal: %d\n", tumor.rms_mapQ, normal.rms_mapQ);
  }
	

}
