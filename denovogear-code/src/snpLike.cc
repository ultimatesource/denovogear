#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include "parser.h"
#include "lookup.h"
#include <string.h>
#include "newmatap.h"
#include "newmatio.h"

using namespace std;


// Calculate DNM and Null PP
void trio_like_snp( qcall_t child, qcall_t mom, qcall_t dad, int flag, 
  vector<vector<string > > & tgt, lookup_snp_t & lookup, 
  string op_vcf_f, ofstream& fo_vcf, double pp_cutoff, int RD_cutoff, int& n_site_pass)
{
  // Filter low read depths ( < 10 )
  if (child.depth < RD_cutoff || mom.depth < RD_cutoff || dad.depth < RD_cutoff) {
    return;
  }
  n_site_pass += 1;
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
  char ref_name[50];
  strcpy(ref_name, child.chr); // Name of the reference sequence 
    
  //Load Likelihood matrices L(D|Gm), L(D|Gd) and L(D|Gc)
  for (j = 0; j != 10; ++j)	{ 
    #ifdef DEBUG_ENABLED
    cout<<"\n"<<" mom lik "<<mom.lk[j]<<" "<<pow(10, -mom.lk[j]/10); 
    cout<<" dad lik "<<dad.lk[j]<<" "<<pow(10, -dad.lk[j]/10);
    cout<<" child lik "<<child.lk[j]<<" "<<pow(10, -child.lk[j]/10); 
  	#endif    
    a[j]=pow(10,-mom.lk[j]/10.0); 
  }  
  M<<a;  

  for (j = 0; j != 10; ++j) 
    a[j] = pow(10, -dad.lk[j]/10.0);
  D<<a;

  for (j = 0; j != 10; ++j) 
    a[j] = pow(10, -child.lk[j]/10.0);
  C<<a;

  P = KP(M, D); // 10 * 10
  F = KP(P, C); // 100 * 10

  // Combine transmission probs L(Gc | Gm, Gf)
  T = SP(F, lookup.tp); //100 * 10
  #ifdef DEBUG_ENABLED  
  cout<<"\n\nMOM\n\n";
  cout << setw(10) << setprecision(10) << M;
  cout<<"\n\nDAD\n\n";
  cout << setw(10) << setprecision(10) << D;
  cout<<"\n\nCHILD\n\n";
  cout << setw(10) << setprecision(10) << C;
  cout<<"\n\nP = L(Gm, Gd)\n\n";
  cout << setw(10) << setprecision(10) << P;
  cout<<"\n\nF = L(Gc, Gm, Gd)\n\n";
  cout << setw(10) << setprecision(10) << F;
  cout<<"\n\nT = L(Gc, Gm, Gd) * L(Gc | Gm, Gf)\n\n";
  cout << setw(10) << setprecision(10) << T;
  #endif

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

  #ifdef DEBUG_ENABLED
  cout<<"\n\nlookup.mrate\n\n";
  cout << setw(10) << setprecision(10) << lookup.mrate;
  cout<<"\n\nlookup.norm\n\n";
  cout << setw(10) << setprecision(10) << lookup.norm;
  cout<<"\n\nlookup.denovo\n\n";
  cout << setw(10) << setprecision(10) << lookup.denovo;
  cout<<"\n\nL = L(Gc, Gm, Gd) * L(Gc | Gm, Gf) * L(Gm, Gf)\n\n";
  cout << setw(10) << setprecision(10) << L;
  #endif

  DN=SP(L,lookup.mrate); 

  // Find max likelihood of null configuration
  PP=SP(DN,lookup.norm);   //zeroes out configurations with mendelian error
  maxlike_null = PP.maximum2(i,j);   

  #ifdef DEBUG_ENABLED
  cout<<"\n\nDN\n\n";
  cout << setw(10) << setprecision(20) << DN;
  cout<<"\n\nPP normal\n\n";
  cout << setw(10) << setprecision(20) << PP;
  #endif

  // Find max likelihood of de novo trio configuration
  PP=SP(DN,lookup.denovo);   //zeroes out configurations with mendelian inheritance
  maxlike_denovo=PP.maximum2(k,l);

  #ifdef DEBUG_ENABLED
  cout<<"\n\nPP denovo\n\n";
  cout << setw(10) << setprecision(20) << PP;
  #endif

  denom=DN.sum();

  #ifdef DEBUG_ENABLED
  cout<<"\nmax_normal i "<<i<<" j "<<j<<" GT "<<tgt[i-1][j-1];
  cout<<"\nmax_denovo k "<<k<<" l "<<l<<" GT "<<tgt[k-1][l-1];;
  cout<<"\nDenom is "<<denom;    
  #endif  

  pp_denovo=maxlike_denovo/denom; //denovo posterior probability
  pp_null=1-pp_denovo; //null posterior probability

  // Check for PP cutoff 
  if (pp_denovo > pp_cutoff) 
  {
    cout<<"DENOVO-SNP CHILD_ID: "<<child.id;
    cout<<" chr: "<<ref_name<<" pos: "<<coor<<" ref: "<<mom.ref_base<<" alt: "<<mom.alt;
    cout<<" maxlike_null: "<<maxlike_null<<" pp_null: "<<pp_null<<" tgt_null(child/mom/dad): "<<tgt[i-1][j-1];
    cout<<" snpcode: "<<lookup.snpcode(i,j)<<" code: "<<lookup.code(i,j);
    cout<<" maxlike_dnm: "<<maxlike_denovo<<" pp_dnm: "<<pp_denovo;
    cout<<" tgt_dnm(child/mom/dad): "<<tgt[k-1][l-1]<<" lookup: "<<lookup.code(k,l)<<" flag: "<<flag;
    cout<<" READ_DEPTH child: "<<child.depth<<" dad: "<<dad.depth<<" mom: "<<mom.depth;
    cout<<" MAPPING_QUALITY child: "<<child.rms_mapQ<<" dad: "<<dad.rms_mapQ<<" mom: "<<mom.rms_mapQ;
    cout<<endl;

    if(op_vcf_f != "EMPTY") {
      fo_vcf<<ref_name<<"\t";
      fo_vcf<<coor<<"\t";
      fo_vcf<<".\t";// Don't know the rsID
      fo_vcf<<mom.ref_base<<"\t";
      fo_vcf<<mom.alt<<"\t";
      fo_vcf<<"0\t";// Quality of the Call
      fo_vcf<<"PASS\t";// passed the read depth filter
      fo_vcf<<"RD_MOM="<<mom.depth<<";RD_DAD="<<dad.depth;
      fo_vcf<<";MQ_MOM="<<mom.rms_mapQ<<";MQ_DAD="<<dad.rms_mapQ;  
      fo_vcf<<";NULL_CONFIG="<<tgt[i-1][j-1]<<";PP_NULL="<<pp_null;  
      fo_vcf<<";ML_NULL="<<maxlike_null<<";ML_DNM="<<maxlike_denovo;
      fo_vcf<<";SNPcode="<<lookup.snpcode(i,j)<<";code="<<lookup.code(i,j)<<";\t";
      fo_vcf<<"DNM_CONFIG:PP_DNM:RD:MQ\t";
      fo_vcf<<tgt[k-1][l-1]<<":"<<pp_denovo<<":"<<child.depth<<":"<<child.rms_mapQ; 
      fo_vcf<<"\n";
    }
  }

}
