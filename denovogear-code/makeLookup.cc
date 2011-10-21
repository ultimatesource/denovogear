#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <set>
#include <vector>

#include "makeLookup.h"
#include "lookup.h"
using namespace std;

const double g_kMinLike = 1e-26;
double g_Mrate = 1e-8;
double g_IndelMrate = 1e-9;
const double g_kIndelFDR = 0.05;
double g_PolyRate = 1e-3;
const double g_kIndelTheta = 1e-4; // NOT USED CURRENTLY
double g_t_prob, g_priors[4], g_kDenovorate; 
int g_inf, g_n_u_alleles;
string g_gts;
bool g_dflag, g_nflag;

void readSNPLookup(ofstream& fout, vector<vector<string > > & tgt, 
					lookup_t & lookup, float lines[][1000])
{
	static int snp_index = 0, k = 0, l=0;
	fout<<g_n_u_alleles<<" "<<g_inf<<" "<<g_t_prob<<" "<<g_gts<<" ";
	fout<<g_kDenovorate<<" "<<g_dflag<<" "<<g_nflag;
	for(int i=0; i<4; i++) {
		fout<<" "<<g_priors[i];
	}
	lines[0][l] = g_n_u_alleles;
	lines[1][l] = g_inf;
	lines[2][l] = g_t_prob;
	lines[3][l] = g_kDenovorate;
	lines[4][l] = g_dflag;
	lines[5][l] = g_nflag;
	lines[6][l] = g_priors[0];
	lines[7][l] = g_priors[1];
	lines[8][l] = g_priors[2];
	lines[9][l++] = g_priors[3];
	//cout<<endl<<snp_index<<" "<<k<<" "<<l;
	tgt[k][snp_index++] = g_gts;
	if (snp_index == 10) {
		k++;
		snp_index = 0;
	}		
	fout<<"\n";
}

void readIndelLookup(ofstream& fout, vector<vector<string > > & tgt, 
					  lookup_t & lookup, float lines[][27])
{
	static int indel_index = 0, k = 0, l = 0;
	fout<<g_n_u_alleles<<" "<<g_inf<<" "<<g_t_prob<<" "<<g_gts<<" ";
	fout<<g_kDenovorate<<" "<<g_dflag<<" "<<g_nflag;
	fout<<" "<<g_priors[0]; // ONLY ONE INDEL PRIOR
	lines[0][l] = g_n_u_alleles;
	lines[1][l] = g_inf;
	lines[2][l] = g_t_prob;
	lines[3][l] = g_kDenovorate;
	lines[4][l] = g_dflag;
	lines[5][l] = g_nflag;
	lines[6][l++] = g_priors[0];
	//cout<<endl<<indel_index<<" "<<k;
	tgt[k][indel_index++] = g_gts;
	if (indel_index == 3) {
		k++;
		indel_index = 0;
	}
	fout<<"\n";
}

void getIndelPriors(string g_gts, int n_uniqa)
{
    string m_d_alleles = g_gts.substr(2, 4);
    set<char> m_d_uniq_a;

    m_d_uniq_a.insert(m_d_alleles[0]);
    m_d_uniq_a.insert(m_d_alleles[1]);
    m_d_uniq_a.insert(m_d_alleles[2]);
    m_d_uniq_a.insert(m_d_alleles[3]);
    int m_d_nuniq = m_d_uniq_a.size();

    int nhit = 0;
	for( int j =0; j<4; j++) {
		if (m_d_alleles[j] == 'R')
			nhit++;
	}

	//ref not in the genotypes
	if ((nhit == 0) && (m_d_nuniq == 1)) { //0 copies of ref in parents
	   		g_priors[0] = (3.0/5.0) * (1.0/5.0);
	} 

	else if (nhit > 0) { //ref in the genotypes
	   	if (nhit == 4) { // 4 copies of the ref, site look monomorphic in parents
	   		g_priors[0] = 1;
	   	} 
	    else if (nhit == 3) { // 3 copy of ref in parents
	    	g_priors[0] = (3.0/5.0) * (4.0/5.0) * 0.5;
	    } 
	    else if (nhit == 2) { // two 
			if ((m_d_alleles[0] == m_d_alleles[1]) && (m_d_alleles[2] == m_d_alleles[3])) {
				g_priors[0] = (2.0/5.0) * (1.0/5.0) * 0.5;
			} 
			else if ((m_d_alleles[0] != m_d_alleles[1]) &&  (m_d_alleles[2] != m_d_alleles[3])) { 
				g_priors[0] = (2.0/5.0) * (2.0/5.0);
			} 
			else {
				cout<<"\nref[i] all  all[0]  all[1]  all[2]  all[3]\n"; 
				cout<<"\nException 1 in g_priors\n";
				exit(1);
			}
		}
	    else if (nhit == 1) { //one 
			g_priors[0] = (2.0/5.0) * (2.0/5.0) * 0.5;
		} 
	    else {
		   	cout<<"\nExiting! Exception 2 in g_priors\n";
		   	exit(1);
		}
	}

	if (nhit == 4) {
		g_priors[0] *= g_kIndelFDR;
	} else {
	   	g_priors[0] *= (1 - g_kIndelFDR);
	}

}

void getSNPPriors(string g_gts, int n_uniqa)
{
    string m_d_alleles = g_gts.substr(2, 4);
    char ref[] = { 'A', 'C', 'G', 'T' };
    set<char> m_d_uniq_a;

    m_d_uniq_a.insert(m_d_alleles[0]);
    m_d_uniq_a.insert(m_d_alleles[1]);
    m_d_uniq_a.insert(m_d_alleles[2]);
    m_d_uniq_a.insert(m_d_alleles[3]);
    int m_d_nuniq = m_d_uniq_a.size();
    for (int i=0; i<4; i++) { // LOOP OVER ALL 4 POSSIBLE BASES IN REFERENCE
		if (n_uniqa > 3) {
			g_priors[i] = g_kMinLike; 
			continue;
		}
		if (n_uniqa == 3) {
			g_priors[i] = g_PolyRate * g_PolyRate;
			continue;
		}
	 	
	 	int nhit = 0;
	 	for( int j =0; j<4; j++) {
	 		//cout<<"\nm_d ref"<<m_d_alleles[j]<<" "<<ref[i];
	 		if (m_d_alleles[j] == ref[i])
	 			nhit++;
		}
	
		//ref not in the genotypes
		if (nhit == 0) {
	    	if (m_d_nuniq == 1) { //0 copies of ref in parents
	    		g_priors[i] = 0.002 * (3.0/5.0) * (1.0/5.0);
	    	} 
	    	else { //triallelic  
	    		g_priors[i] = g_PolyRate * g_PolyRate;
	    	} 
		}
		
		//ref in the genotypes
		else if (nhit > 0) {
	    	if (nhit == 4) { // 4 copies of the ref
	    		g_priors[i] = 0.998;
	    	} 
	    	else if (nhit == 3) { // 3 copy of ref in parents
	    		g_priors[i] = 0.002 * (3.0/5.0) * (4.0/5.0) * 0.5;
	    	} 
	    	else if (nhit == 2) { // two 
				if ((m_d_alleles[0] == m_d_alleles[1]) && (m_d_alleles[2] == m_d_alleles[3])) {
					g_priors[i] = 0.002 * (2.0/5.0) * (1.0/5.0) * 0.5;
				} 
				else if ((m_d_alleles[0] != m_d_alleles[1]) &&  (m_d_alleles[2] != m_d_alleles[3])) { 
					g_priors[i] = 0.002 * (2.0/5.0) * (2.0/5.0);
				} 
				else {
					cout<<"\nref[i] all  all[0]  all[1]  all[2]  all[3]\n"; 
					cout<<"\nException 1 in g_priors\n";
					exit(1);
				}
			}
	    	else if (nhit == 1) { //one 
		 		g_priors[i] = 0.002 * (2.0/5.0) * (2.0/5.0) * 0.5;
			} 
	    	else {
		    	cout<<"\nExiting! Exception 2 in g_priors\n";
		    	exit(1);
			}
		}
		else {
			cout<<"\nExiting! Exception 3 in g_priors\n";
			exit(1);
		}
		//cout<<"\n"<<g_priors[i];
    }
}	

void makeSNPLookup(vector<vector<string > > & tgt, lookup_t & lookup)
{
	float lines[10][1000]; 
	string seq1[] = { "A","A","A","A","C","C","C","G","G","T" };
	string seq2[] = { "A","C","G","T","C","G","T","G","T","T" };
	ofstream fout("snp_lookup.txt");
	fout.precision(10);
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
	
	for( int did=0; did < 10; did++) {
		for( int cid=0; cid < 10; cid++) {
			for( int mid=0; mid < 10; mid++) {	
				g_kDenovorate = 1.0 - g_Mrate;
                g_dflag = false; 
                g_nflag = false;
				set<string> u_alleles;
				g_gts = seq1[cid];
				u_alleles.insert(seq1[cid]);
				g_gts += seq2[cid];
				u_alleles.insert(seq2[cid]);
				g_gts += "/";
				g_gts += seq1[mid];
				u_alleles.insert(seq1[mid]);
				g_gts += seq2[mid];
				u_alleles.insert(seq2[mid]);
				g_gts += "/";
				g_gts += seq1[did];
				u_alleles.insert(seq1[did]);
				g_gts += seq2[did];
				u_alleles.insert(seq2[did]);
				g_n_u_alleles = u_alleles.size();
				string alleles = seq1[cid] + seq2[cid] + seq1[mid] + seq2[mid] + seq1[did] + seq2[did];
				//cout<<"\n"<<i<<" GT string: "<<g_gts;
				//cout<<" The number of unique u_alleles is: "<<g_n_u_alleles;
				
				getSNPPriors(alleles, g_n_u_alleles);

				// 4 alleles => abort
				if (g_n_u_alleles == 4) {
					g_t_prob = 0.0;
					g_inf = 9;
                    g_dflag = false; 
                    g_nflag = true;
					readSNPLookup(fout, tgt, lookup, lines);
					continue;
				}

				if (g_n_u_alleles == 3) {
					g_inf = 9;
					if ((seq1[cid] == seq1[did] || seq1[cid] == seq2[did]) 
						&& (seq2[cid] == seq1[mid] || seq2[cid] == seq2[mid])) {
						if (seq1[cid] == seq2[cid])
							g_t_prob = 0.25;
						else if ((seq1[mid] == seq2[mid]) || (seq1[did] == seq2[did]))
							g_t_prob = 0.5;	
						else
							g_t_prob = 0.25;
					}
					else
						g_t_prob = 0;
                    g_dflag = false; 
                    g_nflag = true;
					readSNPLookup(fout, tgt, lookup, lines);
					continue;
				}
				

				// 1. Child is missing data or homozygous
				if (seq1[cid] == seq2[cid]) {
					// 1b. mom and dad homozygous for same allele; class D
		    		if (seq1[did] == seq2[did] &&  seq1[mid] ==  seq2[mid]) {
			    		 g_inf = 6; 
			    		 g_t_prob = 1.0;
			    	}		    
		    		
		    		// 1c. mom homozygous or missing data, dad heterozygous; class E
		    		else if ((seq1[did] != seq2[did]) &&  (seq1[mid] == seq2[mid])) {
		    			g_inf = 7;
		    			g_t_prob = 0.5;
		    		}
		    
		    		// 1d. dad homozygous or missing data, mom heterozygous; class F
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] == seq2[did])) { 
			    		g_inf = 8;
			    		g_t_prob = 0.5;
			    	}
		    
		    		// 1e. both parents het;   class G
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] != seq2[did])) { 
			    		g_inf = 9;
			    		g_t_prob = 0.25;
			    	}
				}
		
		
				// 2. Child is Het 
				else {
		     		g_inf = 9;

		    		// 2a.mom and dad homozygous for same allele;
		    		if ((seq1[mid] == seq2[mid]) && (seq1[did] == seq2[did])) {
		    			g_t_prob = 1.0;
		    		}

		    		//2b. mom homozygous or missing data, dad heterozygous; 
		    		else if ((seq1[did] != seq2[did]) && (seq1[mid] == seq2[mid])) {
		    			g_t_prob = 0.5;
		    		}
		    
		    		// 2c. dad homozygous or missing data, mom heterozygous;
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] == seq2[did])) { 
			    		g_t_prob = 0.5;
			    	}

		    		// 2d. dad het, mom heterozygous;
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] != seq2[did])) {
		    			g_t_prob = 0.5;
		    		}

				}
		
				if (g_inf == 0) {
					cout<<"\nExiting! Problems with assigning SNP trio config.\n";
					exit(1);					
				}

				// Check for MIs, This following routine is always run
				int p1 = 0, p2 = 0, p12 = 0, p22 = 0, sum1 = 0, sum2 = 0;
				bool test2, test3;
				// Check if child allele 1 is in mom
				if ((seq1[cid] == seq2[mid]) || (seq1[cid] == seq1[mid])) {
					p1++;
				}
				// Check if child allele 1 is in dad
				if ((seq1[cid] == seq2[did]) || (seq1[cid] == seq1[did])) {
					p2++;
				}
				// If child allele 1 is not found in either parent test3 = 1
				if ((p1==0) && (p2==0)) {
					test2 = true;
				} else {
					test2 = false;
				}    
				// Check if child allele 2 is in mom
				if (seq2[cid] == seq1[mid] || seq2[cid] == seq2[mid])
					p12++;
				// Check if child allele 2 is in dad
				if (seq2[cid] == seq1[did] || seq2[cid] == seq2[did])
					p22++;
				// If child allele 2 is not found in either parent test3 = 1 
				if ((p12==0) && (p22==0)) {
					test3 = true;
				} else {
					test3 = false;
				}
				
				sum1 = p1 + p12;
				sum2 = p2 + p22;      
				
				if ((test2 == false) || (test3 == false)) { //both child alleles seen in at least one parent
					if (sum1 == 0) { // both child alleles are seen in dad
						g_inf = 1; 
					} else if (sum2 == 0) { // both child alleles are seen in mum
						g_inf = 2;
					}				    
				}
				
				if (test2 == true || test3 == true) {
					g_inf = 3;
				}
				//For configuration inf[k]=1 and inf[k]=2, make parsimony assumption that only one mutation occurred
				if(g_inf == 1) { // both alleles in dad, denovo from mom
					g_t_prob = 1.0;
					g_kDenovorate = g_Mrate;
					if((alleles[0] == 'A' || alleles[1] == 'A' || alleles[0] == 'G' || alleles[1] == 'G') && (alleles[2] == 'C' || alleles[3] == 'T'))
							g_inf = 5; // transversion
					else if((alleles[0] == 'C' || alleles[1] == 'C' || alleles[0] == 'T' || alleles[1] == 'T') && (alleles[2] == 'A' || alleles[3] == 'G'))
							g_inf = 5; // transversion
					else 
							g_inf = 4; // transition	
					g_dflag = true; // normal and denovo flags for DNG
				    g_nflag = false; 
				} else if(g_inf == 2) { // both alleles in mum, denovo from dad
					g_t_prob = 1.0;
					g_kDenovorate = g_Mrate;
					if((alleles[0] == 'A' || alleles[1] == 'A' || alleles[0] == 'G' || alleles[1] == 'G') && (alleles[4] == 'C' || alleles[5] == 'T'))
					   g_inf = 5; // transversion
					else if((alleles[0] == 'C' || alleles[1] == 'C' || alleles[0] == 'T' || alleles[1] == 'T') && (alleles[4] == 'A' || alleles[5] == 'G'))
					   g_inf = 5; // transversion
					else 
					   g_inf = 4; // transition	
					g_dflag = true; // normal and denovo flags for DNG
				    g_nflag = false; 
				} else if (g_inf == 3) {
					g_t_prob = 1.0;
					bool flagA = false, flagC = false, flagG = false, flagT = false;
					for(unsigned int k=0; k < alleles.length(); k++) {
						if(alleles[k] == 'C')
							flagC = true;
						else if (alleles[k] == 'T') 
							flagT = true;
						else if (alleles[k] == 'A') 
							flagA = true;
						else if (alleles[k] == 'G') 
							flagG = true;
					}

					if (seq1[cid] == seq2[cid]) { // #homozygous child
						g_kDenovorate = g_Mrate * g_Mrate;
					} 
    				else if (flagC && flagT) { // transition
	    				g_inf = 4;
	    				g_kDenovorate = g_Mrate;
	    			} 
    				else if (flagA && flagG) { // transition
	    				g_inf = 4;
	    				g_kDenovorate = g_Mrate;
	    			} 
    				else {  // transversion
	    				g_inf = 5;
	    				g_kDenovorate = g_Mrate;
	    			}

				    g_dflag = true; 
				    g_nflag = false; // normal and denovo flags for DNG
				} else {
					g_dflag = false; 
					g_nflag = true;
				}
				readSNPLookup(fout, tgt, lookup, lines);										
			}
		}
	}
	fout.close();
	lookup.snpcode << lines[0];
	lookup.code << lines[1];
	lookup.tp << lines[2];
	lookup.mrate << lines[3];
	lookup.denovo << lines[4];
	lookup.norm << lines[5];
	lookup.aref << lines[6];
	lookup.cref << lines[7];
	lookup.gref << lines[8];
	lookup.tref << lines[9];
}

void makeIndelLookup(vector<vector<string > > & tgt, lookup_t & lookupIndel)
{
	//string gts[] = { "RR", "RD", "DD" };
	float lines[7][27];
	string seq1[] = { "R", "R", "D" };
	string seq2[] = { "R", "D", "D" };
	ofstream fout("indel_lookup.txt");
	fout.precision(10);
	lookupIndel.priors.resize(9, 3);
    lookupIndel.snpcode.resize(9,3);  
    lookupIndel.code.resize(9,3);  
    lookupIndel.tp.resize(9,3);  
    lookupIndel.mrate.resize(9,3);
	lookupIndel.hit.resize(9,3);
    lookupIndel.denovo.resize(9,3);
    lookupIndel.norm.resize(9,3);
	
	for( int did=0; did < 3; did++) {
		for( int cid=0; cid < 3; cid++) {
			for( int mid=0; mid < 3; mid++) {				
				g_kDenovorate = 0;
                g_dflag = false; 
                g_nflag = false;
				set<string> u_alleles;
				g_gts = seq1[cid];
				u_alleles.insert(seq1[cid]);
				g_gts += seq2[cid];
				u_alleles.insert(seq2[cid]);
				g_gts += "/";
				g_gts += seq1[mid];
				u_alleles.insert(seq1[mid]);
				g_gts += seq2[mid];
				u_alleles.insert(seq2[mid]);
				g_gts += "/";
				g_gts += seq1[did];
				u_alleles.insert(seq1[did]);
				g_gts += seq2[did];
				u_alleles.insert(seq2[did]);
				g_n_u_alleles = u_alleles.size();
				string alleles = seq1[cid] + seq2[cid] + seq1[mid] + seq2[mid] + seq1[did] + seq2[did];
				//cout<<"\n"<<" GT string: "<<g_gts;
				//cout<<" The number of unique u_alleles is: "<<g_n_u_alleles;				
				getIndelPriors(alleles, g_n_u_alleles);

				// 1. Child is missing data or homozygous
				if (seq1[cid] == seq2[cid]) {
					// 1b. mom and dad homozygous for same allele; class D
		    		if (seq1[did] == seq2[did] && seq1[mid] ==  seq2[mid]) {
			    		 g_inf = 6; 
			    		 g_t_prob = 1.0;
			    	}		    
		    		
		    		// 1c. mom homozygous or missing data, dad heterozygous; class E
		    		else if ((seq1[did] != seq2[did]) && (seq1[mid] == seq2[mid])) {
		    			g_inf = 7;
		    			g_t_prob = 0.5;
		    		}
		    
		    		// 1d. dad homozygous or missing data, mom heterozygous; class F
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] == seq2[did])) { 
			    		g_inf = 8;
			    		g_t_prob = 0.5;
			    	}
		    
		    		// 1e. both parents het;   class G
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] != seq2[did])) { 
			    		g_inf = 9;
			    		g_t_prob = 0.25;
			    	}
				}
		
		
				// 2. Child is Het 
				else {
		     		g_inf = 9;

		    		// 2a.mom and dad homozygous for same allele;
		    		if ((seq1[mid] == seq2[mid]) && (seq1[did] == seq2[did])) {
		    			g_t_prob = 1.0;
		    		}

		    		//2b. mom homozygous or missing data, dad heterozygous; 
		    		else if ((seq1[did] != seq2[did]) && (seq1[mid] == seq2[mid])) {
		    			g_t_prob = 0.5;
		    		}
		    
		    		// 2c. dad homozygous or missing data, mom heterozygous;
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] == seq2[did])) { 
			    		g_t_prob = 0.5;
			    	}

		    		// 2d. dad het, mom heterozygous;
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] != seq2[did])) {
		    			g_t_prob = 0.5;
		    		}
				}
		
				if (g_inf == 0) {
					cout<<"\nExiting! Problems with assigning SNP trio config.\n";
					exit(1);					
				}

				// Check for MIs
				int p1 = 0, p2 = 0, p12 = 0, p22 = 0;
				bool test2, test3;

				if ((seq1[cid] == seq2[mid]) || (seq1[cid] == seq1[mid])) {
					p1++;
				}
				if ((seq1[cid] == seq2[did]) || (seq1[cid] == seq1[did])) {
					p2++;
				}
				if ((p1==0) && (p2==0)) {
					test2 = true;
				} else {
					test2 = false;
				}    
		
				if ((seq2[cid] == seq1[mid]) || (seq2[cid] == seq2[mid]))
					p12++;
				if ((seq2[cid] == seq1[did]) || (seq2[cid] == seq2[did]))
					p22++;
				if ((p12==0) && (p22==0)) {
					test3 = true;
				} else { 
					test3 = false;
				}
				
				p1 = p1+p12;
				p2 = p2+p22;      				
				if ((test2 == false) || (test3 == false)) {
					if (p1==0) { 
						g_inf = 1;
					} else if (p2==0) {
						g_inf = 2;
					}				    
				}
				if (test2 == true || test3 == true) {
					g_inf = 3;
				}

				if (g_inf < 6) {
					g_t_prob = 0;
				}

				if (g_inf==3) {
					
					g_t_prob = 1;
		    		if (cid == 0) { // two mutation
			    		g_kDenovorate = 2; 
			    	}
		    		if (cid == 1) { // single mutation
						//cout<<"\ng inf == 3"<<g_gts;
			    		g_kDenovorate = 1;
			    	} 
		    		if (cid == 2) { // two mutation
		    			g_kDenovorate = 2; 
		    		}
		    		g_dflag = true; 
				    g_nflag = false;
		    	} else {
			    	g_dflag = false; 
				    g_nflag = true;
			    }   
			    
				readIndelLookup(fout, tgt, lookupIndel, lines);
			}
		}
	}

	fout.close();
	
	lookupIndel.snpcode << lines[0];
	lookupIndel.code << lines[1];
    lookupIndel.tp << lines[2];    
    lookupIndel.hit << lines[3];
    lookupIndel.denovo << lines[4];
    lookupIndel.norm << lines[5];
    lookupIndel.priors << lines[6];

}

void makeLookup(string table_type, double Mrate, double IndelMrate, double PolyRate, 
				vector<vector<string > > & tgt, lookup_t & lookup) 
{
	g_Mrate = Mrate;
	g_IndelMrate = IndelMrate;
	g_PolyRate = PolyRate;
	if(table_type == "point") {
		makeSNPLookup(tgt, lookup);
	} 
	else if(table_type == "indel") {
		makeIndelLookup(tgt, lookup);
	}
	else {
		cout<<"\nInvalid Table Type !";
		exit(1);
	}
}
