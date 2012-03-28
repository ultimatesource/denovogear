/* 

   Implements parental phasing by looking at the genotypes of parents at phasing 
   sites within a specified window, possible to infer parent as reads are from same
   molecule as DNM, uses Samtools to pull the required reads 

   Usage - ./denovogear phaser --dnm dnm_f --pgt pgt_f --bam bam_f --window [1000]

*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>



using namespace std;
const int g_kFileNameLength = 500;

extern "C" { // from utils/sam_view.c
  int main_samview(int argc, char *argv[]);// to get the individual reads
}

// Format Seq string according to cigar operation
void formatSeq(string& seq, int& pos, char op, int num)
{
  //cout<<"\nOriginal seq is "<<seq;
  //cout<<"\nop "<<op<<" num "<<num<<" pos "<<pos;
  switch (op) {
  case 'S': { // soft clip, erase characters
    seq.erase(pos, num);
    break;
  }
  case 'M': { // match, retain as is
    pos+=num;
    break;
  }
  case 'D': { // deletion, insert a '-'
    seq.insert(pos, num, '-');
    pos+=num;
    break;
  }
  case 'I': { // insertion, remove insert
    seq.erase(pos, num);
    break;
  }
  default: { // unknown cigar operation
    cerr<<op;
    assert(false);
    break;
  } 
  }

  //cout<<"\nModified seq is "<<seq;
}

// Extract CIGAR from reads, expand reads and get the denovo and hap base
int processReads(char* reads_file, long dnm_pos, long hap_pos, string gt1, string gt2, char variant_base)
{
  ifstream fin1(reads_file, ios::in);	
  if (!fin1.is_open()) {
    return -1;
  }
  else {			
    string l;
    getline(fin1, l);
    map<string, char> read_dnm, read_hap;
    map<string, int> pair_count;
    while (fin1.good()) {
      //cout<<"\n\nNext Read";
      istringstream iss(l);
      vector<string> fields;
      string token;
      iss>>token;
      while(iss) {
	     fields.push_back(token);
	     iss>>token;
      }			

      string qname = fields[0];
      int flag = atoi(fields[1].c_str());
      
      // filter unmapped reads, 3rd bit is 1
      if((flag & 4) == 4 ) {
        //cout<<"\n\nUnmapped fragment ! Skipping";
        //cout<<" flag "<<flag;
        getline(fin1, l);
        continue;
      }

      // FILTER QC FAILED READS
      if((flag & 512) == 512 ) {
        cout<<"\n\nFailed Quality Control";
        cout<<" flag "<<flag;
        cout<<"\n\n";
        getline(fin1, l);
        continue;
      }

      // CONSIDER PCR DUPLICATES AS WE ARE NOT GENOTYPING
      /*if((flag & 1024) == 1024 ) {
        //cout<<"\n\nPCR DUPLICATE";
        //cout<<" flag "<<flag;
        //cout<<"\n\n";
        getline(fin1, l);
        continue;
      }*/
        
      string cigar = fields[5];
      string old_cigar = cigar;
      //cerr<<"\n NEW CIGAR "<<cigar<<"TT";
      string seq = fields[9];
      string old_seq = seq;
      long pos = atol(fields[3].c_str());
      int cig_index = 0;
      //cout<<"\nseq is "<<seq;
      while(cigar != "" && cigar != "*") {
      	//cout<<"\n CIGAR "<<cigar;
      	stringstream ss;
      	ss<<cigar;
      	char ch;
      	int num;
      	ss>>num>>ch; // get the single cigar operation
      	string operation;
      	stringstream ss2;
      	ss2<<num<<ch;
      	operation = ss2.str();
      	//cout<<"\nxtracted is "<<operation<<endl;
      	formatSeq(seq, cig_index, ch, num);
      	int prune = operation.length();
      	//				cerr<<"\n "<<prune;
      	cigar = cigar.substr(prune);
      	//				cerr<<" NEW CIGAR "<<cigar<<"**";
      }
      int len = seq.length();	
      string bases;
      long offset = 0;
      if((dnm_pos >= pos) && ((pos + len) > dnm_pos)) {
	      offset = dnm_pos - pos;
	      //cout<<"\nOffset "<<offset<<"BASE "<<seq[offset];
	      read_dnm[qname] = seq[offset];			
      }
      if((hap_pos >= pos) && ((pos + len) > hap_pos)) {
	      offset = hap_pos - pos;
	      read_hap[qname] = seq[offset];				
      }		
      if(read_dnm.count(qname) > 0 && read_hap.count(qname) > 0)	{
      	string bases;
      	//cout<<"\n PAIR BASES 1 is "<<bases;
        if ( read_dnm[qname] == 'G') {
          /*cout<<"\nInside G";
          cout<<"\nqname"<<qname;
          cout<<"\nold seq is "<<old_seq;
          cout<<"\nnew seq is "<<seq;
          cout<<"\ncigar is "<<old_cigar;
          cout<<"\nflag is "<<flag;
          cout<<"\n\nx is "<<x;
          cout<<"\nhap pos "<<hap_pos<<" pos "<<pos<<" offset "<<offset;*/
        }
      	bases += read_dnm[qname];
      	bases += read_hap[qname];
      	//bases[0] = read_dnm[qname];
      	//bases[1] = read_hap[qname];
      	//cout<<"BASES 2 is "<<bases;
      	if(pair_count.count(bases) == 0) {
      	  pair_count[bases] = 1;
      	}
      	else
      	  pair_count[bases]++;
      }
      getline(fin1, l);
    }
    fin1.close();
    //cout<<"\n";
    map<string, int>::iterator it;
    cout<<endl<<"\tHAP POS "<<hap_pos<<" p1: "<<gt1<<" p2: "<<gt2;
    cout<<" Number of denovo-phasing pairs found: "<<pair_count.size();
    for(it = pair_count.begin(); it != pair_count.end(); it++) {
      char dnm_b = (*it).first[0];
      char hap_b = (*it).first[1];
      int count = (*it).second;
      string parent_of_origin = "N/A";

      if ((hap_b == gt1[0]) || (hap_b == gt1[1])) {
        if ((hap_b != gt2[0]) && (hap_b != gt2[1]))
          if (variant_base == dnm_b)
            parent_of_origin = "p1";
          else 
            parent_of_origin = "p2";
      }

      else if ((hap_b == gt2[0]) || (hap_b == gt2[1])) {
        if ((hap_b != gt1[0]) && (hap_b != gt1[1]))
         if (variant_base == dnm_b)
            parent_of_origin = "p2";
          else 
            parent_of_origin = "p1";
      }

      if ((hap_b == gt1[0]) && (hap_b == gt1[1])) {
          if (variant_base == dnm_b)
            parent_of_origin = "p1";
          else 
            parent_of_origin = "p2";
      }

      else if ((hap_b == gt2[0]) && (hap_b == gt2[1])) {
         if (variant_base == dnm_b)
            parent_of_origin = "p2";
          else 
            parent_of_origin = "p1";
      }
     

      cout<<"\n\t\tBase at DNM position: "<<dnm_b<<" Base at phasing position: "<<hap_b<<"\t";
      cout<<" INFERRED PARENT OF ORIGIN for DNM: "<<parent_of_origin<<" SUPPORTING READ COUNT: "<<count;
      
    }
    return 0;
  }
}

// Call samtools to get the reads from the bam file
void getReadsFromBAM(char* bam_f, int chr1, long dnm_pos, long hap_pos, char* temp_file)
{
  long interval = hap_pos-dnm_pos;
  if (interval < 0) 
    interval = -interval;
  string reads_f = "reads_temp.txt"; // store reads in this file
  char program[] = "samtools\0";
  char command[] = "view\0";
  //char op1[] = "-S\0"; // SAM format option
  char op[] = "-o\0";
  strcpy(temp_file, "XXXXXX");
  int fd;
  fd = mkstemp(temp_file);
  //char file[100];
  //strcpy(file, bam_f);
  char region[30];
  long int start = hap_pos, end = dnm_pos;
  if (dnm_pos < hap_pos) {
    start = dnm_pos;
    end = hap_pos;
  }
  sprintf(region, "%d:%ld-%ld", chr1, start, end);
  char* argv1[] = {program, command, op, temp_file, bam_f, region};
  int argc1 = sizeof(argv1)/sizeof(char *);
  main_samview(argc1-1, argv1+1);
  return;
}

// Main
int mainPhaser( int argc, char* argv[])
{
  char DNM_f[g_kFileNameLength] = "EMPTY", parentGT_f[g_kFileNameLength] = "EMPTY", bam_f[g_kFileNameLength] = "EMPTY";
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
	window = atoi(optarg); // size of window for phasing sites( > insert size )
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
	// DNM FILE FORMAT - chr posn inherited_base variant_base
  if(fin1.is_open()) { // PARSE THROUGH DNMs
    int chr1;
    long dnm_pos;
    char inherited_base, variant_base;
    char token1[200]; // token for parsing
    int line_n1 = 0;
    fin1.getline(token1, 20, '\t' );
    while (fin1.good()) {
      line_n1++;	
      
      if (strstr(token1, "chr")) {
            char* token1_p = token1 + 3;
            chr1 = atoi(token1_p);
      } 
      else      
        chr1 = atoi(token1);// chromosome of DNM, PROBLEM WITH SEX CHR
      fin1.getline(token1, 20, '\t');
      dnm_pos = atol(token1);// posn of DNM
      fin1.getline(token1, 20, '\t');
      inherited_base = token1[0];// inherited base at DNM position
      fin1.getline(token1, 20);
      variant_base = token1[0];// variant base i.e DNM base
      cout<<"\n\nDNM_pos "<<chr1<<":"<<dnm_pos<<"\tINHERITED "<<inherited_base<<"\tVARIANT "<<variant_base;

      if(inherited_base==variant_base) {
        cout<<"\nInherited base same as variant base";
        cout<<"\nline_n1 "<<line_n1<<" chr1 "<<chr1<<" dnm_pos "<<dnm_pos;  
        cout<<"\nExiting!";
        exit(1);
      }
      
      fstream fin2(parentGT_f, ios::in);
      if(fin2.is_open()) { // PARSE Parental GTs
      	int chr2;
      	string gt1, gt2, gt3;
      	long hap_pos;
      	char token2[200]; // token for parsing
      	int line_n2 = 0;
        fin2.getline(token2, 20, '\t');
      	while (fin2.good()) {
      	  line_n2++;	      	  
          if (strstr(token2, "chr")) {
            char* token2_p = token2 + 3;
            //cout<<"\ntoken2_p "<<token2_p;
            chr2 = atoi(token2_p);
            //cout<<"\nchr2 "<<chr2;
          } 
          else
      	    chr2 = atoi(token2); // chr of phasing site, PROBLEM WITH SEX CHR
      	  fin2.getline(token2, 20, '\t');
      	  hap_pos = atol(token2); // position of phasing site
      	  fin2.getline(token2, 20, '\t');
      	  gt1 = token2; // genotype of first parent
      	  fin2.getline(token2, 20, '\t');
      	  gt2 = token2; // genotype of second parent
      	  fin2.getline(token2, 20);
          gt3 = token2; // genotype of child
          long dist = hap_pos - dnm_pos; // distance b/w DNM and phasing site
      	  if (dist < 0)
      	    dist = -dist;
      	  //cout<<"\n\tline_n2 "<<line_n2<<" dist "<<dist;
      	  //cout<<" hap_pos "<<chr2<<":"<<hap_pos<<" gt1 "<<gt1<<" gt2 "<<gt2<<" gt3 "<<gt3; 
      	  if ((chr2 > chr1) || 
              ((chr2 == chr1) && ((hap_pos - dnm_pos) > window))) 
      	    break;
      	  else if ((chr2 == chr1) && (dist <= window) && (hap_pos != dnm_pos)) {
      	    if (gt1[0] == 'N' || gt2[0] == 'N') // GT not available
      	      continue;
      	    if (gt1[0] != gt1[1] && gt2[0] != gt2[1]) // both parents het, not informative
      	      continue;
            // ignore triallelic case for now
            if (gt3[0] == gt3[1]) // child hom, not informative
              continue; 
      	    else { 
      	      //cout<<"\nInformative position "<<chr1<<" "<<dnm_pos<<" "<<hap_pos;
      	      char temp_file[20];
              getReadsFromBAM(bam_f, chr1, dnm_pos, hap_pos, temp_file); // get reads corresponding to the positions
              processReads(temp_file, dnm_pos, hap_pos, gt1, gt2, variant_base);
              remove(temp_file);
      	    }
      	  }
          fin2.getline(token2, 20, '\t');
      	}
  	    fin2.close();
  	    //cout<<"\nThe number of lines read GT file is "<<line_n2<<" DNM is "<<line_n1;
      }
      else {
  	    cout<<"\nUnable to open parent GT file: "<<parentGT_f<<" ! Exiting!\n";
  	    exit(1);
      }
      fin1.getline(token1, 20, '\t' );
    }
    fin1.close();
    cout<<"\nThe number of lines read DNM file is "<<line_n1;
  }
  else {
    cout<<"\nUnable to open DNM file: "<<DNM_f<<" ! Exiting!\n";
    exit(1);
  }	
  cout<<"\n";
  return 0;
}