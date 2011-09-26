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

#define MRATE 5e-7
#define MIN_READ_DEPTH 0
#define MIN_READ_DEPTH_INDEL 10
#define MIN_MAPQ 40 

int read_indelLookup (vector<vector<string> > &tgtIndel, lookup_t & lookupIndel);
void trio_like_indel(indel_t *child, indel_t *mom, indel_t *dad, int flag, vector<vector<string> > &tgtIndel, lookup_t & lookupIndel);


using namespace std;

/* commands */
int read_indelLookup( vector<vector<string > >  & tgtIndel, lookup_t & lookupIndel)
{
            string line;
            int tmp[27];
            float tmp2[27];
            float tmp3[27];
            float tmp4[27];
            float tmp5[27];
            float tmp6[27];
            float tmp7[27];


  vector<string> tmpIndel;
    lookupIndel.priors.resize(9, 3);
    lookupIndel.snpcode.resize(9,3);  
    lookupIndel.code.resize(9,3);  
    lookupIndel.tp.resize(9,3);  
    lookupIndel.mrate.resize(9,3);
    lookupIndel.denovo.resize(9,3);
    lookupIndel.norm.resize(9,3);
  
  for (int l=0;l<3;l++)
     tmpIndel.push_back("NA");
  for (int l=0;l<9;l++)
     tgtIndel.push_back(tmpIndel);

                //vector<string> tmp4; 
            float blah;
            int i=0,j=0,k=0;
            string blah2;
            
            ifstream lookup_file("lookup_indel.txt");
            if (!lookup_file) {
                cerr <<"cannot open lookup table"<<endl; 
                exit(1);
            }
            
            while (getline(lookup_file,line)) {
                istringstream iss(line);
                iss >> blah;
                tmp[k]= (int) blah;
                iss >> blah;
                tmp2[k]=blah;
                iss >> blah;
                tmp3[k]=blah;
                iss >> blah2;
                tgtIndel[i][j++]=blah2;
                
                if (j==3) {
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
            
                k++;
            }
            lookup_file.close();
            lookupIndel.code << tmp2;
            lookupIndel.tp << tmp3;
            lookupIndel.snpcode << tmp;
            lookupIndel.mrate << tmp4;
            lookupIndel.denovo << tmp5;
            lookupIndel.norm << tmp6;
            lookupIndel.priors << tmp7;
            

            cerr<<"\n Read "<<k <<" elements from INDEL lookup"<<endl;
            cerr <<" First mrate"<< lookupIndel.mrate(1,1) <<" Last "<< lookupIndel.mrate(9,3)<<endl;
            cerr <<" First "<< lookupIndel.code(1,1) <<" Last "<< lookupIndel.code(9,3)<<endl;
            cerr <<" First "<< tgtIndel[0][0] <<" Last "<< tgtIndel[8][2] <<endl;
            cerr <<" First "<< lookupIndel.priors(1,1) <<" Last "<< lookupIndel.priors(9,3) <<endl;
            return 0; 
        }



void trio_like_indel(indel_t *child,indel_t *mom, indel_t *dad, int flag, vector<vector<string > > & tgtIndel, lookup_t & lookupIndel)
{
    
    
    
    
    Real a[3];   
    Real maxlike_null,maxlike_denovo,pp_null,pp_denovo,denom,numer;       
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
    char ref_name[3];
    strcpy( ref_name, child->chr); // Name of the reference sequence
    

    
    // this is the complex case where child has more than 1 alt, ignore for now
    if (strstr(child->alt, ",")!=0) {
        // printf("\nMultiple alt, position %ld\n", child->pos); // testing - Avinash 
    	return;
    }


    //Currently, samtools creates indel records are created for every trio in the BCF regardless of whether they have data
    //skip sites with both parents missing or child data missing (could still see de novo with data in child and one parent e.g. DD/NN/RR)
    if (child->depth < MIN_READ_DEPTH_INDEL){return;}
    if (mom->depth < MIN_READ_DEPTH_INDEL && dad->depth < MIN_READ_DEPTH_INDEL ) {return;}    


    
    //Load vectors    
    for (j = 0; j != 3; ++j) a[j]=pow(10,-mom->lk[j]/10.);
    M<<a;
    for (j = 0; j != 3; ++j) a[j]=pow(10,-dad->lk[j]/10.);
    D<<a;
    for (j = 0; j != 3; ++j) a[j]=pow(10,-child->lk[j]/10.);
    C<<a;
    
    
	P=KP(M,D);    
	F=KP(P,C);    
	T=SP(F,lookupIndel.tp); //combine with transmission probs 
       
    L  = SP(T, lookupIndel.priors);
	DN = SP(L,lookupIndel.mrate);    
    PP = SP(DN,lookupIndel.norm);   //zeroes out configurations with mendelian error
    maxlike_null = PP.maximum2(i,j);       
    
	//Find max likelihood of de novo trio configuration
	PP=SP(DN,lookupIndel.denovo);   //zeroes out configurations with mendelian inheritance
    maxlike_denovo=PP.maximum2(k,l); 
    
    
    
    //make proper posterior probs
    
    denom=DN.sum();
    numer=PP.sum();
    pp_denovo=maxlike_denovo/denom;
    pp_null=1-pp_denovo;
    
    
    if ( pp_denovo > 0.001 ) {
    	cout<<"\nDENOVO-INDEL CHILD ID: "<<child->id;
    	cout<<" ref_name: "<<ref_name<<" coor: "<<coor<<" ref_base: "<<mom->ref_base<<" ALT: "<<mom->alt;
    	cout<<" maxlike_null: "<<maxlike_null<<" pp_null: "<<pp_null<<" tgt: "<<tgtIndel[i-1][j-1];
    	cout<<" snpcode: "<<lookupIndel.snpcode(i,j)<<" code: "<<lookupIndel.code(i,j);
    	cout<<" maxlike_dnm: "<<maxlike_denovo<<" pp_dnm: "<<pp_denovo;
    	cout<<" tgt: "<<tgtIndel[k-1][l-1]<<" lookup: "<<lookupIndel.code(k,l)<<" flag: "<<flag;
    	printf(" READ_DEPTH child: %d dad: %d mom: %d", child->depth, dad->depth, mom->depth);
    	printf(" MAPPING_QUALITY child: %d dad: %d mom: %d\n", child->rms_mapQ, dad->rms_mapQ, mom->rms_mapQ);
  	}
    
}
