#define MAX_TRIOS 1000
#define ID_LENGTH 100
#define LINE_LENGTH 1000

typedef struct  {
	char fID[ID_LENGTH];	// family ID
	char cID[ID_LENGTH];	// child sampleID
	char dID[ID_LENGTH];	// dad sampleID
	char mID[ID_LENGTH];	// mom sampleID
} Trio;

int parse_ped(const char* ped_file, Trio* trios);
