#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pedParser.h"

#define DEBUG_ENABLED
// Parse PED file and get trio information
void parse_ped(const char* ped_file, Trio** t, Pair** p, int& trio_count, 
               int& pair_count)
{
	Trio* trios;
	Pair* pairs;

	int trios_allocated = 10;
	trios = (Trio*) calloc (trios_allocated, sizeof(Trio));
	if(trios == NULL) {
		printf("\nError allocating memory(1). Exiting!");
		exit(1);
	}

	int pairs_allocated = 10;
	pairs = (Pair*) calloc (pairs_allocated, sizeof(Pair));
	if(pairs == NULL) {
		printf("\nError allocating memory(1). Exiting!");
		exit(1);
	}

	FILE *fp;
	fp = fopen( ped_file, "r");
	if (fp == NULL) {
         printf("\nUnable to open PED file, Exiting !\n");
         exit(0);
	}
	
	char line[LINE_LENGTH];
	char col[4][ID_LENGTH];	
	
	while (fgets (line, 10000, fp) != NULL) {
		sscanf(line, "%s %s %s %s", col[0], col[1], col[2], col[3]);     
      
		if ((strcmp(col[2], "0") != 0) && (strcmp(col[3], "0") != 0)) {
			

			strcpy (trios[trio_count].fID, col[0]); // Get Family ID
			strcpy (trios[trio_count].cID, col[1]); // Get Child ID
			strcpy (trios[trio_count].dID, col[2]); // Get Dad ID
          	strcpy (trios[trio_count].mID, col[3]); // Get Mom ID
			trio_count++; // Increase trio count

			//check memory usage
			if(trio_count > trios_allocated) {
				trios_allocated += 10;
				Trio* temp_trios = (Trio*) realloc(trios, trios_allocated*sizeof(Trio));	
				if(temp_trios == NULL) {
					printf("\nError allocating memory(2). Exiting!");
					exit(1);
				}
				else 
					trios=temp_trios;							 
			}

			#ifdef DEBUG_ENABLED
				printf("\nLine\t%s", line);
				printf("\nfID\t%s", col[0]);
				printf("\ncID\t%s", col[1]);
				printf("\ndID\t%s", col[2]);
				printf("\nmID\t%s", col[3]);
				printf("\ntrios_allocated\t%d", trios_allocated);
				printf("\ntrio_count\t%d", trio_count);
			#endif
		}
      
        else if ((strcmp(col[2], "0") != 0) && (strcmp(col[3], "0") == 0)) {
			strcpy (pairs[pair_count].pairID, col[0]); // Get Pair ID
			strcpy (pairs[pair_count].tumorID, col[1]); // Get Tumor Sample ID
			strcpy (pairs[pair_count].normalID, col[2]); // Get Normal Sample ID
			pair_count++; // Increase pair count

            //check memory usage
			if(pair_count > pairs_allocated) {
				pairs_allocated += 10;
				Pair* temp_pairs = (Pair*) realloc(pairs, pairs_allocated*sizeof(Pair));	
				if(temp_pairs == NULL) {
					printf("\nError allocating memory(2). Exiting!");
					exit(1);
				}
				else 
					pairs=temp_pairs;							 
			}
        }
	}
	
	fclose(fp);
	*p=pairs;
	*t=trios;
}
