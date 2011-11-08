void makePairedLookup(double pairMrate, vector<vector<string > > & tgt, lookup_pair_t & lookupIndel)
{
	
	float d_flag[100], n_flag[100], n_alleles[100]; 
	string seq1[] = { "A","A","A","A","C","C","C","G","G","T" };
	string seq2[] = { "A","C","G","T","C","G","T","G","T","T" };
	ofstream fout("pair_lookup.txt");
	fout.precision(10); 
	
	lookupIndel.priors.resize(10, 10); 
	lookupIndel.denovo.resize(10, 10);///////////
    lookupIndel.norm.resize(10, 10);
    static int snp_index = 0, k = 0, l=0;

	// Iterate through all genotypes
	for( int tum = 0; tum < 10; tum++) {
		for( int nor = 0; nor < 10; nor++) {	
			int index = tum*10 + nor;
           	d_flag[index] = false; // denovo flag
			n_flag[index] = true; // normal flag
            g_kDenovorate = 1.0 - pairMrate;
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

			if (seq1[nor] != seq1[tum]) 
				if (seq1[nor] != seq2[tum]) 
					if (seq2[nor] != seq1[tum]) {
						d_flag[index] = true; // denovo flag
						n_flag[index] = false; // normal flag					
						g_kDenovorate = pairMrate * pairMrate;
					}
					else {
						d_flag[index] = true; // denovo flag
						n_flag[index] = false; // normal flag					
						g_kDenovorate = pairMrate;
					}	
				else if (seq2[nor] != seq1[tum]) 
					if (seq2[nor] != seq2[tum]) {
						d_flag[index] = true; // denovo flag
						n_flag[index] = false; // normal flag						
						g_kDenovorate = pairMrate * pairMrate;
					}
					else {
						d_flag[index] = true; // denovo flag
						n_flag[index] = false; // normal flag				
						g_kDenovorate = pairMrate;
					}

		}
	}
	lookupIndel.snpcode << n_alleles;
	lookup.denovo << d_flag;
	lookup.norm << n_flag;
}