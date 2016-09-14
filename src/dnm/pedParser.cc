/*
 * Copyright (c) 2010, 2011 Genome Research Ltd.
 * Copyright (c) 2012, 2013 Donald Conrad and Washington University in St. Louis
 * Authors: Donald Conrad <dconrad@genetics.wustl.edu>,
 * Avinash Ramu <aramu@genetics.wustl.edu>
 * This file is part of DeNovoGear.
 *
 * DeNovoGear is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <string.h>
#include <iostream>
#include <vector>

#include "pedParser.h"
using namespace std;

void parse_ped(string ped_file, std::vector<Trio> &trios, std::vector<Pair> &pairs) {

	FILE *fp;
    fp = fopen(ped_file.c_str(), "r");
    if(fp == NULL) {
        printf("\nUnable to open PED file, Exiting !\n");
        exit(0);
    }

    char line[LINE_LENGTH];
    char col[4][ID_LENGTH];
    int nline = 0;

    unsigned int position = 0;
    while(fgets(line, 10000, fp) != NULL) {
        if(nline++ > 10000) {
            printf("Error allocating memory for number of lines in PED file.");
            exit(1);
        }
        sscanf(line, "%s %s %s %s", col[0], col[1], col[2], col[3]);

        if((strcmp(col[2], "0") != 0) && (strcmp(col[3], "0") != 0)) {
        	Trio tmp;
            strcpy(tmp.fID, col[0]);  // Get Family ID
            strcpy(tmp.cID, col[1]);  // Get Child ID
            strcpy(tmp.dID, col[2]);  // Get Dad ID
            strcpy(tmp.mID, col[3]);  // Get Mom ID
            trios.push_back(tmp);

        } else if((strcmp(col[2], "0") != 0) && (strcmp(col[3], "0") == 0)) {
        	Pair tmp;
        	strcpy(tmp.pairID, col[0]);  // Get Pair ID
        	strcpy(tmp.tumorID, col[1]);  // Get Tumor Sample ID
        	strcpy(tmp.normalID, col[2]);  // Get Normal Sample ID
        	pairs.push_back(tmp);

        }
    }
    fclose(fp);

}


/*
// Parse PED file and get trio information
void parse_ped1(string ped_file, Trio **t, Pair **p, int &trio_count,
               int &pair_count) {
    Trio *trios;

    int trios_allocated = 10000;
    trios = new(nothrow) Trio[10000];
    if(trios == NULL) {
        printf("\nError allocating memory(1). Exiting!");
        exit(1);
    }

    int pairs_allocated = 10000;
    Pair *pairs = new(nothrow) Pair[10000];
    if(pairs == NULL) {
        printf("\nError allocating memory(1). Exiting!");
        exit(1);
    }

    FILE *fp;
    fp = fopen(ped_file.c_str(), "r");
    if(fp == NULL) {
        printf("\nUnable to open PED file, Exiting !\n");
        exit(0);
    }

    char line[LINE_LENGTH];
    char col[4][ID_LENGTH];
    int nline = 0;

    while(fgets(line, 10000, fp) != NULL) {
        if(nline++ > 10000) {
            printf("Error allocating memory for number of lines in PED file.");
            exit(1);
        }
        sscanf(line, "%s %s %s %s", col[0], col[1], col[2], col[3]);

        if((strcmp(col[2], "0") != 0) && (strcmp(col[3], "0") != 0)) {
            strcpy(trios[trio_count].fID, col[0]);  // Get Family ID
            strcpy(trios[trio_count].cID, col[1]);  // Get Child ID
            strcpy(trios[trio_count].dID, col[2]);  // Get Dad ID
            strcpy(trios[trio_count].mID, col[3]);  // Get Mom ID
            trio_count++;
            //check memory usage
            if(trio_count > trios_allocated) {
                trios_allocated += 10;
                Trio *temp_trios = (Trio *) realloc(trios, trios_allocated * sizeof(Trio));
                if(temp_trios == NULL) {
                    printf("\nError allocating memory(2). Exiting!");
                    exit(1);
                } else {
                    trios = temp_trios;
                }
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
        } else if((strcmp(col[2], "0") != 0) && (strcmp(col[3], "0") == 0)) {
            strcpy(pairs[pair_count].pairID, col[0]);  // Get Pair ID
            strcpy(pairs[pair_count].tumorID, col[1]);  // Get Tumor Sample ID
            strcpy(pairs[pair_count].normalID, col[2]);  // Get Normal Sample ID
            pair_count++; // Increase pair count

            //check memory usage
            if(pair_count > pairs_allocated) {
                pairs_allocated += 10;
                Pair *temp_pairs = (Pair *) realloc(pairs, pairs_allocated * sizeof(Pair));
                if(temp_pairs == NULL) {
                    printf("\nError allocating memory(2). Exiting!");
                    exit(1);
                } else {
                    pairs = temp_pairs;
                }
            }
        }
    }
    fclose(fp);
    *p = pairs;
    *t = trios;
}
*/
