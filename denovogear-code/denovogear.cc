#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "bcf.h"
#include "bgzf.h"
#include "denovogear.h"

#define WANT_STREAM       // include iostream and iomanipulators
#include "newmatap.h"
#include "newmatio.h"

using namespace std;

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

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


int parse(char* ped_file, char* bcf_file) 
{
   	printf("\nPED file : %s, BCF file : %s", ped_file, bcf_file);   
   	int trio_count = parse_ped(ped_file, trios); // Parse PED files and read in trios
   	printf ("\nThe number of trios in the ped file : %d\n", trio_count); 
   
	
   
   	/*
   	for(i=0; i<trio_count; i++) {
		printf("\n%s %s %s %s\n", trios[i].fID, trios[i].cID, trios[i].mID, trios[i].dID);
   	}
   	*/
     	
	vector<vector<string > > tgt;
	lookup_t lookup;
    read_lookup(tgt, lookup);
    
	vector<vector<string > > tgtIndel;
	lookup_t lookupIndel;   	
    read_indelLookup(tgtIndel, lookupIndel);
    qcall_t mom_snp, dad_snp, child_snp;
   	indel_t mom_indel, dad_indel, child_indel;
   	
   	
   	
   	//printf("\nfID 0 %s", trios[0].fID);
	bcf_hdr_t *hout, *hin;
	bcf_t *bp = vcf_open(bcf_file, "rb");
	hin = hout = vcf_hdr_read(bp);
	bcf1_t *b;
	b = static_cast<bcf1_t *> (calloc(1, sizeof(bcf1_t)));  
	while ( vcf_read(bp, hin, b) > 0 ) {
   		int j = 0, flag =0;
   		//printf("Position Number %d", pos++);
   		for ( j=0; j<trio_count; j++) {
   			int is_indel = bcf_2qcall(hout, b, trios[j],  &mom_snp, &dad_snp, &child_snp, &mom_indel, &dad_indel, &child_indel, flag);
   			//cout<<"\nis_indel is "<<is_indel;
   			
   			if ( is_indel == 0 ) {
   				trio_like_snp(child_snp, mom_snp, dad_snp, flag, tgt, lookup);   			
   			}   
   			else if ( is_indel == 1 ) {
			    trio_like_indel(&child_indel, &mom_indel, &dad_indel, flag, tgtIndel, lookupIndel);  
   			}
   			else if ( is_indel < 0 ) {
   				printf("\n BCF PARSING ERROR !  %d", is_indel);
   				printf("\n Exiting !");
   				exit(1);
   			}
   		}   		  
	}
	
	bcf_hdr_destroy(hin);
	bcf_destroy(b);
    bcf_close(bp);
   	return 0;
   
}


int main(int argc, char *argv[])
{

  parse(argv[1], argv[2]);
  printf("\nDone !\n");
  return 0;
}
