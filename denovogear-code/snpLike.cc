#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#define WANT_STREAM       // include iostream and iomanipulators
#include "bcf.h"
#include "parser.h"
#include "lookup.h"


#include "bgzf.h"
#include "newmatap.h"
#include "newmatio.h"

#ifdef use_namespace
using namespace RBD_LIBRARIES;
#endif

#ifndef VERSION
#define VERSION "dummy"
#endif


using namespace std;
/* commands */
int read_lookup(vector<vector<string > > & tgt, lookup_t & lookup)
{
  string line;
  int tmp[1000];
  float tmp2[1000];
  float tmp3[1000];
  float tmp4[1000];
  float tmp5[1000];
  float tmp6[1000];
  float tmp7[1000];
  float tmp8[1000];
  float tmp9[1000];
  float tmp10[1000];

  //vector<string> tmp4; 
  float blah;
  int i=0,j=0,k=0;
  string blah2;
  
  vector<string> tmpV;
  for (int l=0;l<10;l++)
      tmpV.push_back("NA");
  for (int l=0;l<100;l++)
      tgt.push_back(tmpV);
  
  lookup.aref.resize(100,10);
  lookup.cref.resize(100,10);
  lookup.gref.resize(100,10);
  lookup.tref.resize(100,10);  
  lookup.snpcode.resize(100,10);  
  lookup.code.resize(100,10);  
  lookup.tp.resize(100,10);  
  lookup.mrate.resize(100,10);
  lookup.denovo.resize(100,10);
  lookup.norm.resize(100,10);
  

  ifstream lookup_file("lookup.txt");
  if (!lookup_file) {
    cerr <<"cannot open lookup table"<<endl; 
    exit(1);
  }
  
  cerr<<"Entered read lookup"<<endl;
  
  while (getline(lookup_file,line)) {
    istringstream iss(line);
    iss >> blah;
    tmp[k]= (int) blah;
    iss >> blah;
    tmp2[k]=blah;
    iss >> blah;
    tmp3[k]=blah;
    iss >> blah2;
    tgt[i][j++]=blah2;

    if (j==10) {
      i++;
      j=0;
    }
      
    iss >> blah;
    tmp4[k]=blah;
    iss >> blah;
    tmp5[k]=blah;
    iss >> blah;
    tmp6[k]=blah;
    iss >> blah;
    tmp7[k]=blah;
    iss >> blah;
    tmp8[k]=blah;
    iss >> blah;
    tmp9[k]=blah;
    iss >> blah;
    tmp10[k]=blah;
    k++;
  }
  lookup_file.close();
  lookup.code << tmp2;
  lookup.tp << tmp3;
  lookup.snpcode << tmp;
  lookup.mrate << tmp4;
  lookup.denovo << tmp5;
  lookup.norm << tmp6;
  lookup.aref << tmp7;
  lookup.cref << tmp8;
  lookup.gref << tmp9;
  lookup.tref << tmp10;

  cerr<<" Read "<<k <<" elements from lookup"<<endl;
  cerr <<" First mrate"<< lookup.mrate(1,1) <<" Last "<< lookup.mrate(100,10)<<endl;
  cerr <<" First "<< lookup.code(1,1) <<" Last "<< lookup.code(100,10)<<endl;
  cerr <<" First "<< tgt[0][0] <<" Last "<< tgt[99][9] <<endl;
  cerr <<" First "<< lookup.tref(1,1) <<" Last "<< lookup.tref(100,10) <<endl;
  return 0; 
}


void trio_like_snp(qcall_t child, qcall_t mom, qcall_t dad, int flag, vector<vector<string > > & tgt, lookup_t & lookup)
{
  Real a[10];   
  Real maxlike_null,maxlike_denovo,pp_null,pp_denovo,denom,numer;   

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

  
  //Load vectors

 
  for (j = 0; j != 10; ++j)	{ 
  	//cerr<<"\n position "<<coor<<" mom lik "<<mom.lk[j]<<" dad lik "<<dad.lk[j]<<" child lik "<<child.lk[j]; 
  	a[j]=pow(10,-mom.lk[j]/10.); 
  }
  M<<a;
  
  for (j = 0; j != 10; ++j) 
  	a[j]=pow(10,-dad.lk[j]/10.);
  D<<a;
  
  for (j = 0; j != 10; ++j) 
  	a[j]=pow(10,-child.lk[j]/10.);
  C<<a;
        

  P=KP(M,D);
  F=KP(P,C);
  T=SP(F,lookup.tp); //combine with transmission probs 


  switch(mom.ref_base){
    case 'A':
	  L=SP(T,lookup.aref); break;

    case 'C':
	  L=SP(T,lookup.cref); break;

    case 'G':
	  L=SP(T,lookup.gref); break;

    case 'T':
	  L=SP(T,lookup.tref); break;

    default: L=T; break;
  }

  DN=SP(L,lookup.mrate);
  PP=SP(DN,lookup.norm);   //zeroes out configurations with mendelian error
  maxlike_null = PP.maximum2(i,j);   
  
  //Find max likelihood of de novo trio configuration
  PP=SP(DN,lookup.denovo);   //zeroes out configurations with mendelian inheritance
  maxlike_denovo=PP.maximum2(k,l); 

  //make proper posterior probs
  denom=DN.sum();
  numer=PP.sum();
  pp_denovo=maxlike_denovo/denom;
  pp_null=1-pp_denovo;

  if ( pp_denovo > 0.001 ) {
    cout<<"\nDENOVO-SNP CHILD ID: "<<child.id;
    cout<<" ref_name: "<<ref_name<<" coor: "<<coor<<" ref_base: "<<mom.ref_base;
    cout<<" maxlike_null: "<<maxlike_null<<" pp_null: "<<pp_null<<" tgt: "<<tgt[i-1][j-1];
    cout<<" snpcode: "<<lookup.snpcode(i,j)<<" code: "<<lookup.code(i,j);
    cout<<" maxlike_dnm: "<<maxlike_denovo<<" pp_dnm: "<<pp_denovo;
    cout<<" tgt: "<<tgt[k-1][l-1]<<" lookup: "<<lookup.code(k,l)<<" flag: "<<flag;
    printf(" READ_DEPTH child: %d dad: %d mom: %d", child.depth, dad.depth, mom.depth);
    printf(" MAPPING_QUALITY child: %d dad: %d mom: %d\n", child.rms_mapQ, dad.rms_mapQ, mom.rms_mapQ);
  }
     

}


