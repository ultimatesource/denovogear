/*
Make lookup table for paired samples, SNPs and Indels. Calculate prior probabilities, transmission
probabilities and denovo, normal flags for all genotype configurations. 
*/
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <set>

#include "makeLookup.h"

//#define LOOKUP_ENABLED 1 // enable this flag if you want to generate lookup tables

#ifdef LOOKUP_ENABLED
	ofstream fout_snp("snp_lookup.txt");
	ofstream fout_pair("pair_lookup.txt");
	ofstream fout_indel("indel_lookup.txt");
#endif



using namespace std;

const double g_kMinLike = 1e-26; // Minimum Likelihood value
double g_Mrate = 0; // SNP mutation rate
double g_IndelMrate = 0; // Indel mutation rate
const double g_kIndelFDR = 0.05;
double g_PolyRate = 0;
double g_t_prob, g_priors[4], g_kDenovorate; 
int g_inf, g_n_u_alleles; // Transition/Transversion, Number of unique alleles
string g_gts; // Genotype string of trio
bool g_dflag, g_nflag; // Denovo flag, Normal Flag
int ta_c = 0, ta_c2 = 0, ta_c3 = 0;
int g_khit = 0;

// Write to SNP Lookup table file and struct
//void readSNPLookup(ofstream& fout, vector<vector<string > > & tgt,
	 //float lines[][1000]) -- OLD
void readSNPLookup(vector<vector<string > > & tgt,
	 float lines[][1000])
{
	static int snp_index = 0, k = 0, l=0;
	#ifdef LOOKUP_ENABLED
		fout_snp<<g_n_u_alleles<<" "<<g_inf<<" "<<g_t_prob<<" "<<g_gts<<" ";
		fout_snp<<g_kDenovorate<<" "<<g_dflag<<" "<<g_nflag;
		for(int i=0; i<4; i++) {
			fout_snp<<" "<<g_priors[i];
		}
	#endif
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
	tgt[k][snp_index++] = g_gts;
	if (snp_index == 10) {
		k++;
		snp_index = 0;
	}	
	#ifdef LOOKUP_ENABLED	
		fout_snp<<"\n";
	#endif
}

// Write to Indel lookup table file and struct
//void readIndelLookup(ofstream& fout, vector<vector<string > > & tgt, 
	//float lines[][27]) - OLD
void readIndelLookup(vector<vector<string > > & tgt, 
	float lines[][27])     
{
	static int indel_index = 0, k = 0, l = 0;

	#ifdef LOOKUP_ENABLED	
		fout_indel<<g_n_u_alleles<<" "<<g_inf<<" "<<g_t_prob<<" "<<g_gts<<" ";
		fout_indel<<g_kDenovorate<<" "<<g_dflag<<" "<<g_nflag;
		fout_indel<<" "<<g_priors[0]; // ONLY ONE INDEL PRIOR
	#endif
	lines[0][l] = g_n_u_alleles;
	lines[1][l] = g_inf;
	lines[2][l] = g_t_prob;
	lines[3][l] = g_khit;
	lines[4][l] = g_dflag;
	lines[5][l] = g_nflag;
	lines[6][l++] = g_priors[0];
	tgt[k][indel_index++] = g_gts;
	if (indel_index == 3) {
		k++;
		indel_index = 0;
	}
	#ifdef LOOKUP_ENABLED
		fout_indel<<"\n";
	#endif
}

// Calculate Indel priors
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

// Calculate SNP priors
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
      
        if(m_d_nuniq == 3) {
          g_priors[i] = 0.002 * 0.002 / 414;  // prior split equally amongst all triallelic cases
          continue;
        }
          
        else if (n_uniqa == 3) { // 3rd allele in child case
            //cout<<"\nc1 "<<g_gts;
			g_priors[i] = g_PolyRate * g_PolyRate;
			continue;
		}
	 	
	 	int nhit = 0; 
	 	for( int j = 0; j<4; j++) {
	 		if (m_d_alleles[j] == ref[i])
	 			nhit++;
		}
	
		//ref not in the genotypes
		if (nhit == 0) { // 0 copies of ref in parents
	    	if (m_d_nuniq == 1) { // biallelic
	    		g_priors[i] = 0.995 * 0.002 * (3.0/5.0) * (1.0/5.0);
	    	} 
	    	else if (m_d_nuniq == 2) { // triallelic 
              g_priors[i] = 0.002 * 0.002 / 414; // prior split equally amongst all triallelic cases
            } 
            else { 
              g_priors[i] = g_PolyRate * g_PolyRate;
	    	}          
		}
		
		//ref in the genotypes
		else if (nhit > 0) {
	    	if (nhit == 4) { // 4 copies of the ref
	    		g_priors[i] = 0.995 * 0.998;
	    	} 
	    	else if (nhit == 3) { // 3 copies of ref in parents
	    		g_priors[i] = 0.995 * 0.002 * (3.0/5.0) * (4.0/5.0) * 0.5;
	    	} 
	    	else if (nhit == 2) { // 2 copies of the ref in parents
				if ((m_d_alleles[0] == m_d_alleles[1]) && (m_d_alleles[2] == m_d_alleles[3])) {
					g_priors[i] = 0.995 * 0.002 * (2.0/5.0) * (1.0/5.0) * 0.5;
				} 
				else if ((m_d_alleles[0] != m_d_alleles[1]) &&  (m_d_alleles[2] != m_d_alleles[3])) { 
					g_priors[i] = 0.995 * 0.002 * (2.0/5.0) * (2.0/5.0);
				} 
				else {
					cout<<"\n"<<ref[i]<<" "<<g_gts; 
					cout<<"\nException 1 in g_priors\n";
					exit(1);
				}
			}
	    	else if (nhit == 1) { //one 
		 		g_priors[i] = 0.995 * 0.002 * (2.0/5.0) * (2.0/5.0) * 0.5;
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
    }
}	

// Make Indel lookup table
void makeIndelLookup(double IndelMrate, double PolyRate, 
	vector<vector<string > > & tgt, lookup_indel_t & lookupIndel)
{
	#ifdef LOOKUP_ENABLED
		fout_indel.precision(10);
	#endif

	g_IndelMrate = IndelMrate;
	g_PolyRate = PolyRate;

	float lines[7][27];
	string seq1[] = { "R", "R", "D" };
	string seq2[] = { "R", "D", "D" };
	
	lookupIndel.priors.resize(9, 3);
    lookupIndel.snpcode.resize(9,3);  
    lookupIndel.code.resize(9,3);  
    lookupIndel.tp.resize(9,3);  
    lookupIndel.mrate.resize(9,3);// Now Indel mrate is calculated based on length of indel dynamically
    lookupIndel.hit.resize(9,3);
    lookupIndel.denovo.resize(9,3);
    lookupIndel.norm.resize(9,3);
	
	for( int did=0; did < 3; did++) {
		for( int cid=0; cid < 3; cid++) {
			for( int mid=0; mid < 3; mid++) {				
				g_khit = 0;
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
				getIndelPriors(alleles, g_n_u_alleles);

				// Child is missing data or homozygous
				if (seq1[cid] == seq2[cid]) {
					// Mom and Dad homozygous for same allele
		    		if (seq1[did] == seq2[did] && seq1[mid] ==  seq2[mid]) {
			    		 g_inf = 6; 
			    		 g_t_prob = 1.0;
			    	}		    
		    		
		    		// Mom homozygous or missing data, dad heterozygous
		    		else if ((seq1[did] != seq2[did]) && (seq1[mid] == seq2[mid])) {
		    			g_inf = 7;
		    			g_t_prob = 0.5;
		    		}
		    
		    		// Dad homozygous or missing data, mom heterozygous
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] == seq2[did])) { 
			    		g_inf = 8;
			    		g_t_prob = 0.5;
			    	}
		    
		    		// Both parents heterozygous
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] != seq2[did])) { 
			    		g_inf = 9;
			    		g_t_prob = 0.25;
			    	}
				}
		
            
                // Child is Het 
                else {
                      g_inf = 9;
                      // Mom and dad homozygous for same allele
                      if ((seq1[mid] == seq2[mid]) && (seq1[did] == seq2[did])) {
                        g_t_prob = 1.0;
                      }

                      // Mom homozygous or missing data, dad heterozygous
                      else if ((seq1[did] != seq2[did]) && (seq1[mid] == seq2[mid])) {
                          g_t_prob = 0.5;
                      }
              
                      // Dad homozygous or missing data, mom heterozygous
                      else if ((seq1[mid] != seq2[mid]) && (seq1[did] == seq2[did])) { 
                          g_t_prob = 0.5;
                      }

                      // Dad het, mom heterozygous
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
              
                // Check for C1 in M1 and M2
				if ((seq1[cid] == seq2[mid]) || (seq1[cid] == seq1[mid])) {
					p1++;
				}
                // Check for C1 in D1 and D2
				if ((seq1[cid] == seq2[did]) || (seq1[cid] == seq1[did])) {
					p2++;
				}
                // C1 not found in M & D => test2 is true
				if ((p1==0) && (p2==0)) {
					test2 = true;
				} else {
					test2 = false;
				}    
		
                // Check for C2 in M1 and M2
				if ((seq2[cid] == seq1[mid]) || (seq2[cid] == seq2[mid]))
					p12++;
                // Check for C2 in D1 and D2
				if ((seq2[cid] == seq1[did]) || (seq2[cid] == seq2[did]))
					p22++;
                // C2 not found in M & D => test3 is true
				if ((p12==0) && (p22==0)) {
					test3 = true;
				} else { 
					test3 = false;
				}
				
				p1 = p1+p12;
				p2 = p2+p22;      
                
                // Both child alleles seen in one parent
				if ((test2 == false) || (test3 == false)) {
					if (p1==0) { // Both alleles seen in dad
						g_inf = 1;
					} else if (p2==0) { // Both alleles seen in mom
						g_inf = 2;
					}				    
				}
                // C1 or C2 not found in M and D
				if (test2 == true || test3 == true) {
					g_inf = 3;
				}

				if (g_inf < 6) {
					g_t_prob = 0;
				}
                // DNM
				if (g_inf==3) {					
					g_t_prob = 1;
		    		if (cid == 0) { // two mutations
			    		g_khit = 2; 
			    	}
		    		if (cid == 1) { // single mutation
			    		g_khit = 1;
			    	} 
		    		if (cid == 2) { // two mutations
		    			g_khit = 2; 
		    		}
		    		g_dflag = true; 
				    g_nflag = false;
		    	} else {
			    	g_dflag = false; 
				    g_nflag = true;
			    }   
			    
				//readIndelLookup(fout, tgt, lines);
				readIndelLookup(tgt, lines);
			}
		}
	}

	#ifdef LOOKUP_ENABLED
		fout_indel.close();
	#endif
	
	lookupIndel.snpcode << lines[0];
	lookupIndel.code << lines[1];
    lookupIndel.tp << lines[2];    
    lookupIndel.hit << lines[3];
    lookupIndel.denovo << lines[4];
    lookupIndel.norm << lines[5];
    lookupIndel.priors << lines[6];

}
// Create SNP lookup
void makeSNPLookup(double SNPMrate, double PolyRate, 
	vector<vector<string > > & tgt, lookup_snp_t & lookup)
{
	#ifdef LOOKUP_ENABLED
		fout_snp.precision(10);
	#endif

	g_Mrate = SNPMrate;
	g_PolyRate = PolyRate;

	float lines[10][1000]; 
	string seq1[] = { "A","A","A","A","C","C","C","G","G","T" };
	string seq2[] = { "A","C","G","T","C","G","T","G","T","T" };
	
	
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

				// 4 unique alleles 
				if (g_n_u_alleles == 4) {
					g_t_prob = 0.0;
					g_inf = 9;
                    g_dflag = false; 
                    g_nflag = true;
					//readSNPLookup(fout, tgt, lines);
					readSNPLookup(tgt, lines);
					continue;
				}
                
                // triallelic case
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
                      g_dflag = false; 
                      g_nflag = true;
					}
					else {
					  if (seq1[cid] != seq1[did] && seq1[cid] != seq2[did] && seq1[cid] != seq1[mid] && seq1[cid] != seq2[mid] 
							&& seq2[cid] != seq1[did] && seq2[cid] != seq2[did] && seq2[cid] != seq1[mid] && seq2[cid] != seq2[mid]) // TWO MUTATIONS
						g_kDenovorate = g_Mrate * g_Mrate;
					  else
					  	g_kDenovorate = g_Mrate; // ONE MUTATION
                      g_t_prob = 1;
                      g_dflag = true; 
                      g_nflag = false;
                    }                   
					readSNPLookup(tgt, lines);
					continue;
				}
				

				// Child is missing data or homozygous
				if (seq1[cid] == seq2[cid]) {
					// Mom and dad homozygous for same allele
		    		if (seq1[did] == seq2[did] &&  seq1[mid] ==  seq2[mid]) {
			    		 g_inf = 6; 
			    		 g_t_prob = 1.0;
			    	}		    
		    		
		    		// Mom homozygous or missing data, dad heterozygous
		    		else if ((seq1[did] != seq2[did]) &&  (seq1[mid] == seq2[mid])) {
		    			g_inf = 7;
		    			g_t_prob = 0.5;
		    		}
		    
		    		// Dad homozygous or missing data, mom heterozygous
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] == seq2[did])) { 
			    		g_inf = 8;
			    		g_t_prob = 0.5;
			    	}
		    
		    		// Both parents het
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] != seq2[did])) { 
			    		g_inf = 9;
			    		g_t_prob = 0.25;
			    	}
				}
		
		
				// Child is Heterozygous 
				else {
		     		g_inf = 9;

		    		// Mom and Dad homozygous for same allele
		    		if ((seq1[mid] == seq2[mid]) && (seq1[did] == seq2[did])) {
		    			g_t_prob = 1.0;
		    		}

		    		// Mom homozygous or missing data, dad heterozygous
		    		else if ((seq1[did] != seq2[did]) && (seq1[mid] == seq2[mid])) {
		    			g_t_prob = 0.5;
		    		}
		    
		    		// Dad homozygous or missing data, mom heterozygous
		    		else if ((seq1[mid] != seq2[mid]) && (seq1[did] == seq2[did])) { 
			    		g_t_prob = 0.5;
			    	}

		    		// Dad heterozygous, mom heterozygous;
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
				// If child allele 1 is not found in either parent test3 is true
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
				// If child allele 2 is not found in either parent test3 is true 
				if ((p12==0) && (p22==0)) {
					test3 = true;
				} else {
					test3 = false;
				}
				
				sum1 = p1 + p12;
				sum2 = p2 + p22;      
              
				// both child alleles seen in one parent
				if ((test2 == false) || (test3 == false)) { 
					if (sum1 == 0) { // both child alleles seen in dad
						g_inf = 1; 
					} else if (sum2 == 0) { // both child alleles seen in mum
						g_inf = 2;
					}				    
				}
				
				if (test2 == true || test3 == true) {
					g_inf = 3;
				}
              
				// both alleles in dad, denovo from mom
                if(g_inf == 1) {  // AA/BB/AA
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
				} else if(g_inf == 2) { // both alleles in mum, denovo from dad AA/AA/BB
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
				} else if (g_inf == 3) { // DNM, the ones we are interested in 
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

					if (seq1[cid] == seq2[cid]) { // homozygous child
						g_kDenovorate = g_Mrate * g_Mrate;
						//cout<<"\nAlleles is: "<<alleles;
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
				    g_nflag = false; 
				} else {
					g_dflag = false; 
					g_nflag = true;
				}
				//readSNPLookup(fout, tgt, lines);	- OLD
				readSNPLookup(tgt, lines);									
			}
		}
	}

	#ifdef LOOKUP_ENABLED
		fout_snp.close();
	#endif
	
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

// Lookup table for paired samples
void makePairedLookup(double pairMrate, vector<vector<string > > & tgt, lookup_pair_t & lookup)
{
	#ifdef LOOKUP_ENABLED
		fout_pair.precision(10);
	#endif
	
	double d_flag[100], n_flag[100], n_alleles[100], priors[100]; 
	string seq1[] = { "A","A","A","A","C","C","C","G","G","T" };
	string seq2[] = { "A","C","G","T","C","G","T","G","T","T" };
	
	

	lookup.priors.resize(10, 10); 
	lookup.denovo.resize(10, 10);
    lookup.norm.resize(10, 10);
    lookup.snpcode.resize(10, 10);

	// Iterate through all genotypes
	for( int tum = 0; tum < 10; tum++) {
		for( int nor = 0; nor < 10; nor++) {	
			int index = tum*10 + nor;
           	d_flag[index] = false; // denovo flag
			n_flag[index] = true; // normal flag
            priors[index] = 1.0 - pairMrate;
			set<string> u_alleles;
			g_gts = seq1[nor];
			u_alleles.insert(seq1[nor]);
			g_gts += seq2[nor];
			u_alleles.insert(seq2[nor]);
			g_gts += "/";
			g_gts += seq1[tum];
			u_alleles.insert(seq1[tum]);
			g_gts += seq2[tum];
			u_alleles.insert(seq2[tum]);			
			n_alleles[index] = u_alleles.size(); // number of unique alleles
			string alleles = seq1[nor] + seq2[nor] + seq1[tum] + seq2[tum];
			tgt[tum][nor] = g_gts; // genotype string

			// set prior when tumor different from normal sample
			if (seq1[nor] != seq1[tum]) {
				if (seq1[nor] != seq2[tum]) {
					if (seq2[nor] != seq1[tum]) {
						d_flag[index] = true; // denovo flag
						n_flag[index] = false; // normal flag					
						priors[index] = pairMrate * pairMrate;
					}
					else {
						d_flag[index] = true; // denovo flag
						n_flag[index] = false; // normal flag					
						priors[index] = pairMrate;
					}
				}	
				else if (seq2[nor] != seq1[tum]) {
					if (seq2[nor] != seq2[tum]) {
						d_flag[index] = true; // denovo flag
						n_flag[index] = false; // normal flag						
						priors[index] = pairMrate * pairMrate;
					}
					else {
						d_flag[index] = true; // denovo flag
						n_flag[index] = false; // normal flag				
						priors[index] = pairMrate;
					}
				}
			}
			else if (seq2[nor] != seq2[tum]) {
						d_flag[index] = true; // denovo flag
						n_flag[index] = false; // normal flag						
						priors[index] = pairMrate;
			}

			#ifdef LOOKUP_ENABLED
				fout_pair<<n_alleles[index]<<" "<<g_gts<<" "<<d_flag[index];
				fout_pair<<" "<<n_flag[index]<<" "<<priors[index]<<"\n";
			#endif

		}
	}
	#ifdef LOOKUP_ENABLED
		fout_pair.close();
	#endif

	lookup.snpcode << n_alleles;
	lookup.denovo << d_flag;
	lookup.norm << n_flag;
	lookup.priors << priors;
}

