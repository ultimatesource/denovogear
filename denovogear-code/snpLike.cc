#include <vector>
#include <iostream>
#include <string.h>
#include "parser.h"
#include "lookup.h"
#include "newmatap.h"
#include "newmatio.h"

#define MIN_READ_DEPTH_SNP 10
 
using namespace std;


// Calculate DNM and Null PP
void trio_like_snp(qcall_t child, qcall_t mom, qcall_t dad, int flag, 
  vector<vector<string > > & tgt, lookup_snp_t & lookup)
{
  Real a[10];   
  Real maxlike_null, maxlike_denovo, pp_null, pp_denovo, denom;   
  Matrix M(1,10);
  Matrix C(10,1);
  Matrix D(10,1);
  Matrix P(10,10);
  Matrix F(100,10);
  Matrix L(100,10);
  Matrix T(100,10);
  Matrix DN(100,10);
  Matrix PP(100,10);
  int i,j,k,l;  
  int coor = child.pos;
  char ref_name[3];
  strcpy( ref_name, child.chr); // Name of the reference sequence
  
  // Filter low read depths ( < 10 )
  if (child.depth < MIN_READ_DEPTH_SNP ||
      mom.depth < MIN_READ_DEPTH_SNP || dad.depth < MIN_READ_DEPTH_SNP) {
    return;
  }
  
  //Load Likelihood matrices L(D|Gm), L(D|Gd) and L(D|Gc)
  for (j = 0; j != 10; ++j)	{ 
  	//cerr<<"\n position "<<coor<<" mom lik "<<mom.lk[j]<<" dad lik "<<dad.lk[j]<<" child lik "<<child.lk[j]; 
  	a[j]=pow(10,-mom.lk[j]/10.); 
  }
  M<<a;
  
  for (j = 0; j != 10; ++j) 
  	a[j] = pow(10, -dad.lk[j]/10.);
  D<<a;
  
  for (j = 0; j != 10; ++j) 
  	a[j] = pow(10, -child.lk[j]/10.);
  C<<a;
        
  P = KP(M, D); // 10 * 10
  F = KP(P, C); // 100 * 10
  // Combine transmission probs L(Gc | Gm, Gf)
  T = SP(F, lookup.tp); //100 * 10

  // Combine prior L(Gm, Gf) 
  switch(mom.ref_base) {
    case 'A':
	  L = SP(T,lookup.aref); break;

    case 'C':
	  L = SP(T,lookup.cref); break;

    case 'G':
	  L = SP(T,lookup.gref); break;

    case 'T':
	  L = SP(T,lookup.tref); break;

    default: L=T; break;
  }

  DN=SP(L,lookup.mrate);

  // Find max likelihood of null configuration
  PP=SP(DN,lookup.norm);   //zeroes out configurations with mendelian error
  maxlike_null = PP.maximum2(i,j);   
  
  // Find max likelihood of de novo trio configuration
  PP=SP(DN,lookup.denovo);   //zeroes out configurations with mendelian inheritance
  maxlike_denovo=PP.maximum2(k,l); 

  denom=DN.sum();

  pp_denovo=maxlike_denovo/denom; // denovo posterior probability
  pp_null=1-pp_denovo; // null posterior probability

  // Check for PP cutoff 
  if ( pp_denovo > 0.01 ) {
    cout<<"\nDENOVO-SNP CHILD ID: "<<child.id;
    cout<<" ref_name: "<<ref_name<<" coor: "<<coor<<" ref_base: "<<mom.ref_base<<" ALT: "<<mom.alt;
    cout<<" maxlike_null: "<<maxlike_null<<" pp_null: "<<pp_null<<" tgt: "<<tgt[i-1][j-1];
    cout<<" snpcode: "<<lookup.snpcode(i,j)<<" code: "<<lookup.code(i,j);
    cout<<" maxlike_dnm: "<<maxlike_denovo<<" pp_dnm: "<<pp_denovo;
    cout<<" tgt: "<<tgt[k-1][l-1]<<" lookup: "<<lookup.code(k,l)<<" flag: "<<flag;
    printf(" READ_DEPTH child: %d dad: %d mom: %d", child.depth, dad.depth, mom.depth);
    printf(" MAPPING_QUALITY child: %d dad: %d mom: %d\n", child.rms_mapQ, dad.rms_mapQ, mom.rms_mapQ);
  }
	
}
