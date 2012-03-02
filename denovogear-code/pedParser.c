#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pedParser.h"

// Parse PED file and get trio information
void parse_ped(const char* ped_file, Trio* trios, Pair* pairs, int& trio_count, 
               int& pair_count)
{
	
	FILE *fp;
	fp = fopen( ped_file, "r");
	if (fp == NULL) {
         printf("\nUnable to open PED file, Exiting !\n");
         exit(0);
	}
	
	char line[LINE_LENGTH];
	char col[4][ID_LENGTH];	

	
	while (fgets (line , 10000 , fp) != NULL) {
		sscanf(line, "%s %s %s %s", col[0], col[1], col[2], col[3]);     
      
		if ((strcmp(col[2], "0") != 0) && (strcmp(col[3], "0") != 0)) {
			strcpy (trios[trio_count].fID, col[0]); // Get Family ID
			strcpy (trios[trio_count].cID, col[1]); // Get Child ID
			strcpy (trios[trio_count].dID, col[2]); // Get Dad ID
          	strcpy (trios[trio_count].mID, col[3]); // Get Mom ID
			trio_count++; // Increase trio count
			#ifdef DEBUG_ENABLED
				printf("\nLine\t%s", line);
				printf("\nfID\t%s", col[0]);
				printf("\ncID\t%s", col[1]);
				printf("\ndID\t%s", col[2]);
				printf("\nmID\t%s", col[3]);
			#endif
		}
      
        else if ((strcmp(col[2], "0") != 0) && (strcmp(col[3], "0") == 0)) {
          strcpy (pairs[pair_count].pairID, col[0]); // Get Pair ID
          strcpy (pairs[pair_count].tumorID, col[1]); // Get Tumor Sample ID
          strcpy (pairs[pair_count].normalID, col[2]); // Get Normal Sample ID
          pair_count++; // Increase pair count
        }
	}
	
	fclose(fp);
}
