#include <vector>
#include <iostream>
#include <string.h>
#include "parser.h"
#include "lookup.h"
#include "newmatap.h"
#include "newmatio.h"

#define MIN_READ_DEPTH_INDEL 10

using namespace std;

// Calculate DNM and Null PP
void trio_like_indel(indel_t *child,indel_t *mom, indel_t *dad, int flag, 
					 vector<vector<string > > & tgtIndel, 
					 lookup_indel_t & lookupIndel, double mu_scale)
{    
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
    char ref_name[3];
    strcpy( ref_name, child->chr); // Name of the reference sequence
        
    // this is the case where child has more than 1 alt, ignore for now
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

    // mu_scale is the variable used to scale the mu function
    double indel_mrate, new_indel_mrate;
	if(is_insertion) 
		indel_mrate = mu_scale * (-22.8689 - (0.2994 * len_diff)); // insertion
	else
		indel_mrate = mu_scale * (-21.9313 - (0.2856 * len_diff)); // deletion
	indel_mrate = exp(indel_mrate); // antilog of log ratio
	
	for (int j = 1; j <= 9; j++)
		for (int l = 1; l <= 3; l++) {
			new_indel_mrate = pow(indel_mrate, lookupIndel.hit(j,l)); // hit is 0,1 or 2
			lookupIndel.mrate(j,l) = new_indel_mrate; 
		}
	
	//Currently, samtools creates indel records are created for every trio in the BCF regardless of whether they have data
    //skip sites with both parents missing or child data missing (could still see de novo with data in child and one parent e.g. DD/NN/RR)
    if (child->depth < MIN_READ_DEPTH_INDEL ||
		mom->depth < MIN_READ_DEPTH_INDEL || dad->depth < MIN_READ_DEPTH_INDEL) {
		return;
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
    if ( pp_denovo > 0.0001 ) {
    	cout<<"\nDENOVO-INDEL child id: "<<child->id;
    	cout<<" ref_name: "<<ref_name<<" coor: "<<coor<<" ref_base: "<<mom->ref_base<<" ALT: "<<mom->alt;
    	cout<<" maxlike_null: "<<maxlike_null<<" pp_null: "<<pp_null<<" tgt: "<<tgtIndel[i-1][j-1];
    	cout<<" snpcode: "<<lookupIndel.snpcode(i,j)<<" code: "<<lookupIndel.code(i,j);
    	cout<<" maxlike_dnm: "<<maxlike_denovo<<" pp_dnm: "<<pp_denovo;
    	cout<<" tgt: "<<tgtIndel[k-1][l-1]<<" lookup: "<<lookupIndel.code(k,l)<<" flag: "<<flag;
    	printf(" READ_DEPTH child: %d dad: %d mom: %d", child->depth, dad->depth, mom->depth);
    	printf(" MAPPING_QUALITY child: %d dad: %d mom: %d\n", child->rms_mapQ, dad->rms_mapQ, mom->rms_mapQ);
  	}
    
}
