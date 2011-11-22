/* Implements parental phasing by looking at the genotypes of parents at phasing 
sites within a specified window, possible to infer parent as reads are from same
molecule as DNM, uses Samtools to pull the required reads */

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <iostream>


using namespace std;

extern "C" { // from utils/sam_view.c
	int main_samview(int argc, char *argv[]);// to get the individual reads
}

void getReads(char* bam_f, int chr1, long pos1, long pos2)
{
	long interval = pos2-pos1;
	if (interval < 0) 
		interval = -interval;
	string reads_f = "reads_temp.txt"; // store reads in this file
	char temp_string[50];
	sprintf(temp_string, "%d:%ld", chr1, interval);
	cout<<"\nTEMP STRING IS "<<temp_string;
	int argc1 = 3;
	char command[10] = "view";
	char file[10] = "check.txt";
	char region[10] = "1:1000";
	char* argv1[] = {command, file, region};
	cout<<main_samview(argc1, argv1);
	return;
}

int mainPhaser( int argc, char* argv[])
{
	char DNM_f[50] = "EMPTY", parentGT_f[50] = "EMPTY", bam_f[50] = "EMPTY";
	long window = 1000; // default window size is 1000
   
    // Read in Command Line arguments
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {{"dnm", 1, 0, 0}, 
			{"pgt", 1, 0, 1}, {"bam", 1, 0, 2}, {"window", 1, 0, 3},};
		int c = getopt_long (argc, argv, "", long_options, &option_index);
		if (c == -1)
			break;
		switch(c) 
		{      
			case 0:
				strcpy(DNM_f, optarg); // File with list of DNMs to be phased
				break;
			case 1:
				strcpy(parentGT_f, optarg); // File with parental genotypes for phasing sites
				break;	
			case 2:
				strcpy(bam_f, optarg); // BAM file 
				break;
			case 3:
				window = atoi(optarg); // size of window for phasing sites
				break;		
		}
	}
  
	if(!strcmp(DNM_f, "EMPTY") || !strcmp(parentGT_f, "EMPTY") || !strcmp(bam_f, "EMPTY")) {
		cout<<"INPUT ERROR! Please specify the list of DNMs, list of Parental GTs and the BAM file!"
        "\nFor example:\n\t./denovogear phaser --dnm dnm1.txt --pgt pgt1.txt --bam bam1.bam"
        "\nExiting!\n";
		exit(1);
	} 

	cout<<"\nList of DNMs: "<<DNM_f<<", List of parental GTs: "<<parentGT_f<<endl;
	ifstream fin1(DNM_f, ios::in);
	
	if(fin1.is_open()) {
		int chr1;
		long pos1;
		char token1[200]; // token for parsing
		int line_n1 = 0;
		while (fin1.good()) {
			line_n1++;	
			fin1.getline(token1, 20, '\t' );
			chr1 = atoi(token1);// chromosome of DNM
			fin1.getline(token1, 20);
			pos1 = atoi(token1);// posn of DNM
			cout<<"\nline_n1 "<<line_n1<<" chr1 "<<chr1<<" pos1 "<<pos1;
			ifstream fin2(parentGT_f, ios::in);
			if(fin2.is_open()) {
				int chr2;
				string gt1, gt2;
				long pos2;
				char token2[200]; // token for parsing
				int line_n2 = 0;
				while (fin2.good()) {
					line_n2++;	
					fin2.getline(token2, 20, '\t');
					chr2 = atoi(token2); // chr of phasing site
					fin2.getline(token2, 20, '\t');
					pos2 = atoi(token2); // position of phasing site
					fin2.getline(token2, 20, '\t');
					gt1 = token2; // genotype of first parent
					fin2.getline(token2, 20);
					gt2 = token2; // genotype of second parent
					int dist = pos2 - pos1; // distance b/w DNM and phasing site
					if(dist < 0)
						dist = -dist;
					if(chr2 > chr1 || 
					  	(chr2 == chr1 && ((pos2 - pos1) > window))) 
						break;
					else if ((chr2 == chr1) && (dist <= window)) {
						if(gt1[0] == 'N' || gt2[0] == 'N') // GT not available
							continue;
						if(gt1[0] != gt1[1] && gt2[0] != gt2[1]) // both parents het, not informative
							continue;
						else { 
							getReads(bam_f, chr1, pos1, pos2); // get reads corresponding to the positions
							cout<<"\n\tline_n2 "<<line_n2<<" chr2 ";
							cout<<chr2<<" pos2 "<<pos2<<" gt1 "<<gt1<<" gt2 "<<gt2;    
						}
					}

				}
				cout<<"\nThe number of lines read 2 is "<<line_n2;
			}
			else {
				cout<<"\nUnable to open parent GT file: "<<parentGT_f<<" ! Exiting!\n";
				exit(1);
			}
		}
		cout<<"\nThe number of lines read is "<<line_n1;
	}
	else {
		cout<<"\nUnable to open DNM file: "<<DNM_f<<" ! Exiting!\n";
		exit(1);
	}		
	fin1.close();
	cout<<"\n";
	return 0;
}



