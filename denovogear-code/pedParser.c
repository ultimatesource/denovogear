#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "pedParser.h"


int parse_ped(const char* ped_file, Trio* trios)
{
	
	FILE *fp;
	fp = fopen( ped_file, "r");
	if (fp == NULL) {
         printf("Unable to open PED file, Exiting !");
         exit(0);
	}
	
	char line[LINE_LENGTH];
	char fID1[ID_LENGTH];	// family ID
	char cID1[ID_LENGTH];	// child sampleID
	char dID1[ID_LENGTH];	// dad sampleID
	char mID1[ID_LENGTH];	// mom sampleID
	int trio_count = 0;
	
	while (fgets (line , 100 , fp) != NULL) {
		sscanf(line, "%s %s %s %s", fID1, cID1, dID1, mID1);
		if( strcmp(mID1, "0") && strcmp(dID1, "0") ) {
			strcpy( trios[trio_count].fID, fID1);
			strcpy( trios[trio_count].cID, cID1);
			strcpy( trios[trio_count].dID, dID1);
			strcpy( trios[trio_count].mID, mID1);
			trio_count++;
		}
	}
	
	fclose(fp);
	return trio_count;
}
