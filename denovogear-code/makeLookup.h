using namespace std;

void printUsage();
void writeLookup(ofstream& fout, bool snp);
void getIndelPriors(string g_gts1, int n_uniqa);
void getSNPPriors(string g_gts1, int n_uniqa);
void makeSNPLookup();
void makeIndelLookup();
void makeLookup(string table_type, double Mrate, double IndelMrate, double PolyRate);