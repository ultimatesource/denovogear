/*
 * Copyright (c) 2010, 2011 Genome Research Ltd.
 * Copyright (c) 2012, 2013 Donald Conrad and Washington University in St. Louis
 * Copyright (c) 2014 Reed A. Cartwright
 * Authors: Donald Conrad <dconrad@genetics.wustl.edu>,
 *          Avinash Ramu <aramu@genetics.wustl.edu>
 *          Reed A. Cartwright <reed@cartwrig.ht>
 *
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

/*

   Implements parental phasing by looking at the genotypes of parents at phasing
   sites within a specified window, possible to infer parent as reads are from same
   molecule as DNM, uses Samtools to pull the required reads

   Usage - ./denovogear phaser --dnm dnm_f --pgt pgt_f --bam bam_f --window [1000]
   Notes - skips hard clipped reads.
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

#include "version.h"

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <dng/hts/bam.h>
#include <dng/app.h>
#include <dng/task/phaser.h>

using namespace std;
using namespace dng::task;

const int g_kFileNameLength = 500;

extern "C" { // from utils/sam_view.c
    int main_samview(int argc, char *argv[]);// to get reads from BAM
}

// Format Seq string according to cigar operation
int formatSeq(string &seq, int &pos, char op, int num) {
    //cout<<"\nOriginal seq is "<<seq;
    //cout<<"\nop "<<op<<" num "<<num<<" pos "<<pos;
    switch(op) {
    case 'S': { // soft clip, erase characters
        seq.erase(pos, num);
        break;
    }
    case 'M': { // match, retain as is
        pos += num;
        break;
    }
    case 'D': { // deletion, insert a '-'
        seq.insert(pos, num, '-');
        pos += num;
        break;
    }
    case 'I': { // insertion, remove insert
        seq.erase(pos, num);
        break;
    }
    default: { // unknown cigar operation, skip hard clipped reads
        //cerr<<"Unknown CIGAR, skipping read";
        return 1;
    }
    }
    return 0;
    //cout<<"\nModified seq is "<<seq;
}

inline char indexed_char(std::size_t x) {
    static const char table[] = "NACNGNNNTN";
    return table[x];
}


// Extract CIGAR from reads, expand reads and get the denovo and hap base
int processReads(const char *bam_f, std::string &chr1, long dnm_pos,
                 long hap_pos, string &gt1, string &gt2, char variant_base) {

    map<string, char> read_dnm,
        read_hap;// the key is the query name, the char is the base. paired reads have the same query name.
    map<string, int> pair_count;

    // Open bam file for reading
    htsFile *fp = hts_open(bam_f, "r");
    bam_hdr_t *hdr = sam_hdr_read(fp);
    bam1_t *rec = bam_init1();

    // get specific region
    // TODO: C++'ize this
    char region[30];
    long int start = hap_pos, end = dnm_pos;
    if(dnm_pos < hap_pos) {
        start = dnm_pos;
        end = hap_pos;
    }
    sprintf(region, "%s:%ld-%ld", chr1.c_str(), start, end);

    // Find only sequences that fall between the dnm_pos and hap_pos
    // TODO: this is very inefficient, being called multiple times. Just preprocess all the potential positions
    //       and do the query once
    hts_idx_t *sam_indx = sam_index_load(fp, bam_f);
    hts_itr_t *sam_itr = sam_itr_querys(sam_indx, hdr, region);
    while(sam_itr_next(fp, sam_itr, rec) >= 0) {
        const char *qname = bam_get_qname(rec); // field 1
        uint32_t flag = rec->core.flag; // field 2

        // filter unmapped reads, 3rd bit is 1 or mate rescue XT:A:M
        uint8_t *aux;
        bool XT = ((aux = bam_aux_get(rec, "XT")) && (bam_aux2A(aux) == 'M'));
        if((flag & 4) == 4 || XT) {
            continue;
        }

        if((flag & 512) == 512) {
            std::cout << "\n\nFailed Quality Control";
            std::cout << " flag " << flag << std::endl;
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

        uint32_t cigar_len = rec->core.n_cigar;
        uint32_t *cigar = bam_get_cigar(rec); // field 6

        int32_t seq_len = rec->core.l_qseq;
        uint8_t *seq = bam_get_seq(rec); // field 10
        std::string
        seq_formatted; // A string copy of the sequence, will be processed to match the cigar
        for(int a = 0; a < seq_len; a++) {
            // A = 1, C = 2, G = 4, T = 8
            seq_formatted += indexed_char(bam_seqi(seq, a));
        }

        int32_t pos = rec->core.pos; // field 4

        // reformat the sequence based on the cigar string
        int cig_index = 0;
        for(int a = 0; a < cigar_len; a++) {
            int op = bam_cigar_op(cigar[a]);
            int op_len = bam_cigar_oplen(cigar[a]);

            switch(op) {
            case BAM_CSOFT_CLIP:
                // 'S' soft clip, erase characters
                seq_formatted.erase(cig_index, op_len);
                break;
            case BAM_CMATCH:
                // 'M' match, retain as is
                cig_index += op_len;
                break;
            case BAM_CDEL:
                // 'D' deletion, insert a '-'
                seq_formatted.insert(cig_index, op_len, '-');
                break;
            case BAM_CINS:
                // 'I' insertion, remove insert
                seq_formatted.erase(cig_index, op_len);
                break;
            default:
                std::cout << "Unable to handle cigar operation (htslib code " << op << ")." <<
                          std::endl;
                return 1;
            }
        }

        int formatted_len = seq_formatted.length();
        long offset = 0;
        if((dnm_pos >= pos) && ((pos + formatted_len) > dnm_pos)) {
            offset = dnm_pos - pos;
            read_dnm[qname] = seq_formatted[offset - 1];
        }
        if((hap_pos >= pos) && ((pos + formatted_len) > hap_pos)) {
            offset = hap_pos - pos;
            read_hap[qname] = seq_formatted[offset - 1];
        }
        if(read_dnm.count(qname) > 0 && read_hap.count(qname) > 0)  {
            std::string bases;
            bases += read_dnm[qname];
            bases += read_hap[qname];

            if(pair_count.count(bases) == 0) {
                pair_count[bases] = 1;
            } else {
                pair_count[bases]++;
            }
        }
    }

    map<string, int>::iterator it;
    if(pair_count.size() > 0) {
        cout << endl << "\tHAP POS " << hap_pos << " p1: " << gt1 << " p2: " << gt2;
    } else {
        return 0;
    }
    //cout<<" Number of denovo-phasing pairs found: "<<pair_count.size();
    for(it = pair_count.begin(); it != pair_count.end(); it++) {
        char dnm_b = (*it).first[0];
        char hap_b = (*it).first[1];
        int count = (*it).second;
        string parent_of_origin = "N/A";

        if((hap_b == gt1[0]) || (hap_b == gt1[1])) {
            if((hap_b != gt2[0]) && (hap_b != gt2[1])) {
                if(variant_base == dnm_b) {
                    parent_of_origin = "p1";
                } else {
                    parent_of_origin = "p2";
                }
            }
        }

        else if((hap_b == gt2[0]) || (hap_b == gt2[1])) {
            if((hap_b != gt1[0]) && (hap_b != gt1[1])) {
                if(variant_base == dnm_b) {
                    parent_of_origin = "p2";
                } else {
                    parent_of_origin = "p1";
                }
            }
        }

        if((hap_b == gt1[0]) && (hap_b == gt1[1])) {
            if(variant_base == dnm_b) {
                parent_of_origin = "p1";
            } else {
                parent_of_origin = "p2";
            }
        } else if((hap_b == gt2[0]) && (hap_b == gt2[1])) {
            if(variant_base == dnm_b) {
                parent_of_origin = "p2";
            } else {
                parent_of_origin = "p1";
            }
        }

        cout << "\n\t\tBase at DNM position: " << dnm_b << " Base at phasing position: "
             << hap_b << "\t";
        cout << " INFERRED PARENT OF ORIGIN for DNM: " << parent_of_origin <<
             " SUPPORTING READ COUNT: " << count;

    }
    return 10;

}



// Call samtools to get the reads from the bam file
// TODO: DELETE
void getReadsFromBAM(char *bam_f, string chr1, long dnm_pos, long hap_pos,
                     char *temp_file) {
    std::cout << std::endl << "---------------------------------" << std::endl;

    htsFile *fp_b = hts_open(bam_f, "r");
    bam_hdr_t *hdr_b = sam_hdr_read(fp_b);
    bam1_t *rec_b = bam_init1();

    char region[30];
    long int start = hap_pos, end = dnm_pos;
    if(dnm_pos < hap_pos) {
        start = dnm_pos;
        end = hap_pos;
    }
    sprintf(region, "%s:%ld-%ld", chr1.c_str(), start, end);
    std::cout << std::endl << "REGION = " << region << std::endl;

    hts_idx_t *sam_indx = sam_index_load(fp_b, bam_f);
    hts_itr_t *sam_itr = sam_itr_querys(sam_indx, hdr_b, region);
    bam1_t *rec_i = bam_init1();
    while(sam_itr_next(fp_b, sam_itr, rec_b) >= 0) {
        //    std::cout << "HERE" << std::endl;
        //}
        //  while(sam_read1(fp_b, hdr_b, rec_b) > 0) {
        const char *qname_b = bam_get_qname(rec_b);
        std::cout << "qname = " << qname_b << std::endl;

        uint32_t flag_b = rec_b->core.flag;
        std::cout << "flag = " << flag_b << std::endl;

        std::cout << "pos = " << rec_b->core.pos << std::endl;

        uint32_t cigar_b_len = rec_b->core.n_cigar;
        //std::cout << "cigar len = " << cigar_b_len << std::endl;
        uint32_t *cigar_b = bam_get_cigar(rec_b);
        std::cout << "cigar = ";
        for(int a = 0; a < cigar_b_len; a++) {
            std::cout << bam_cigar_oplen(cigar_b[a]) << "[" << bam_cigar_op(
                          cigar_b[a]) << "]";
            //std::cout << "cigar = " << bam_cigar_oplen(cigar_b[0]) << std::endl;
            //std::cout << "cigar op = " << bam_cigar_op(cigar_b[0]) << std::endl;
        }
        std::cout << std::endl;

        uint8_t *seq_b = bam_get_seq(rec_b);
        int32_t seq_b_len = rec_b->core.l_qseq;
        //std::cout << "seq_len = " << seq_b_len << std::endl;
        std::cout << "seq = ";
        for(int a = 0; a < seq_b_len; a++) {
            //std::cout << static_cast<int>(seq_b[a]);
            std::cout << indexed_char(bam_seqi(seq_b, a));
        }
        std::cout << std::endl;

        uint32_t aux_len = bam_get_l_aux(rec_b);
        std::cout << "size of auxiliary data = " << aux_len << std::endl;
        const uint8_t *aux_b = bam_aux_get(rec_b, "XS");//bam_get_aux(rec_b, "XS");
        std::cout << "aux = " << aux_b << std::endl;
        std::cout << "    = " << bam_aux2A(aux_b) << std::endl;
        //std::cout << "aux = ";
        //for(int a = 0; a < aux_len; a++)
        //  std::cout << (char)aux_b[a];
        //std::cout << std::endl;

        std::cout << std::endl;

    }

    //std::cout << "bam_f = " << bam_f << std::endl;
    hts::bam::File child_bam(bam_f, "r");
    //hts::bam::File child_bam(bam_f, "r", "1:1-999999999");
    const bam_hdr_t *hdr = child_bam.header();
    hts::bam::Alignment line;
    //while(child_bam.Read(line)) {
    //  std::cout << line.target_id() << std::endl;
    //}

    //std::cout << "HERE: " << child_bam.Read(line) << std::endl;
    //std::cout << "\tid = " << line.target_id() << std::endl;
    /*
      long interval = hap_pos - dnm_pos;
      if(interval < 0) {
      interval = -interval;
      }
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
      if(dnm_pos < hap_pos) {
      start = dnm_pos;
      end = hap_pos;
      }
      sprintf(region, "%s:%ld-%ld", chr1.c_str(), start, end);
      char *argv1[] = {program, command, op, temp_file, bam_f, region};
      int argc1 = sizeof(argv1) / sizeof(char *);
      main_samview(argc1 - 1, argv1 + 1);
      return;
    */
}

// Main
int Phaser::operator()(Phaser::argument_type &arg) {
    //int main(int argc, char *argv[]) {
    std::cerr << PACKAGE_STRING << " --- SNV Phaser" << std::endl;
    /*
    char DNM_f[g_kFileNameLength] = "EMPTY",
      parentGT_f[g_kFileNameLength] = "EMPTY", bam_f[g_kFileNameLength] = "EMPTY";
    long window = 1000; // default window size is 1000

      // Read in Command Line arguments
      while(1) {
        int option_index = 0;
        static struct option long_options[] = {
    {"dnm", 1, 0, 0},
    {"pgt", 1, 0, 1},
    {"bam", 1, 0, 2},
    {"window", 1, 0, 3},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
        };
        int c = getopt_long(argc, argv, "", long_options, &option_index);
        if(c == -1) {
    break;
        }
        switch(c) {
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
        case 'h':
    cerr << "Usage:\n"
         << "  dng phaser [options]\n\n"
         << "Options:\n"
         << "  --dnm filename: File with list of DNMs to be phased\n"
         << "  --pgt filename: File with parental genotypes for phasing sites\n"
         << "  --bam filename: bam file\n"
         << "  --window size: size of window for phasing sites (> insert size)\n"
         << endl;
    return EXIT_SUCCESS;
        default:
    cerr << "ERROR: Unknown flag.\n";
    return EXIT_FAILURE;
        }
      }
    */

    if(arg.dnm.empty() || arg.pgt.empty() || arg.bam.empty()) {
        throw std::runtime_error("INPUT ERROR! Please specify the list of DNMs, list of Parental GTs and the BAM file!"
                                 "\nFor example:\n\tdng phaser --dnm dnm1.txt --pgt pgt1.txt --bam bam1.bam"
                                 "\nExiting!");
    }

    const char *DNM_f = arg.dnm.c_str();
    const char *parentGT_f = arg.pgt.c_str();
    const char *bam_f = arg.bam.c_str();
    long window = arg.window;

    /*
      if(!strcmp(DNM_f, "EMPTY") || !strcmp(parentGT_f, "EMPTY")
         || !strcmp(bam_f, "EMPTY")) {
        cout << "INPUT ERROR! Please specify the list of DNMs, list of Parental GTs and the BAM file!"
    "\nFor example:\n\tdng phaser --dnm dnm1.txt --pgt pgt1.txt --bam bam1.bam"
    "\nExiting!\n";
        exit(1);
      }
    */

    cout << "\nList of DNMs: " << DNM_f << ", List of parental GTs: " << parentGT_f
         << endl;
    ifstream fin1(DNM_f, ios::in);
    // DNM FILE FORMAT - chr posn inherited_base variant_base
    if(fin1.is_open()) { // PARSE THROUGH DNMs
        string chr1;
        long dnm_pos;
        char inherited_base, variant_base;
        char token1[200]; // token for parsing
        int line_n1 = 0;
        fin1.getline(token1, 20, '\t');
        while(fin1.good()) {
            line_n1++;

            /*if (strstr(token1, "chr")) {
              char* token1_p = token1 + 3;
              chr1 = atoi(token1_p);
              }
              else
              chr1 = atoi(token1);// chromosome of DNM, PROBLEM WITH SEX CHR*/
            //cout <<endl<<"token1 is "<<token1;
            chr1 = token1;
            fin1.getline(token1, 20, '\t');
            dnm_pos = atol(token1);// posn of DNM
            fin1.getline(token1, 20, '\t');
            inherited_base = token1[0];// inherited base at DNM position
            fin1.getline(token1, 20);
            variant_base = token1[0];// variant base i.e DNM base
            cout << "\nDNM_pos " << chr1 << ":" << dnm_pos << "\tINHERITED " <<
                 inherited_base << "\tVARIANT " << variant_base;
            int returnsum = 0;
            if(inherited_base == variant_base) {
                cout << "\nInherited base same as variant base";
                cout << "\nline_n1 " << line_n1 << " chr1 " << chr1 << " dnm_pos " << dnm_pos;
                cout << "\nExiting!";
                exit(1);
            }

            fstream fin2(parentGT_f, ios::in);
            if(fin2.is_open()) { // PARSE Parental GTs
                string chr2;
                string gt1, gt2, gt_c;
                long hap_pos;
                //char token2[200]; // token for parsing
                char token2[200];
                int line_n2 = 0;
                fin2.getline(token2, 20, '\t');
                while(fin2.good()) {

                    /*if (strstr(token2, "chr")) {
                      char* token2_p = token2 + 3;
                      //cout<<"\ntoken2_p "<<token2_p;
                      chr2 = atoi(token2_p);
                      //cout<<"\nchr2 "<<chr2;
                      }
                      else
                      chr2 = atoi(token2); // chr of phasing site, PROBLEM WITH SEX CHR*/

                    line_n2++;
                    chr2 = token2;
                    fin2.getline(token2, 20, '\t');
                    hap_pos = atol(token2); // position of phasing site
                    fin2.getline(token2, 20, '\t');
                    gt_c = token2; // genotype of child
                    fin2.getline(token2, 20, '\t');
                    gt1 = token2; // genotype of first parent
                    fin2.getline(token2, 20, '\n');
                    gt2 = token2; // genotype of second parent
                    fin2.getline(token2, 20, '\t'); // for the next lines chr
                    //cout<<endl<<"c\t"<<gt_c<<"\t"<<"p"<<"\t"<<gt1<<"\t"<<gt2;

                    if(gt1[0] == 'N' || gt2[0] == 'N'  || gt_c[0] == 'N') { // GT not available
                        continue;
                    }
                    if(gt1 == gt2) { // both parents het or both hom, not informative, // ignore triallelic case for now
                        continue;
                    }
                    if(gt_c[0] == gt_c[1]) { // child hom, not informative
                        continue;
                    }


                    //cout<<"\n\tline_n2 "<<line_n2<<" dist "<<dist;
                    //cout<<" hap_pos "<<chr2<<":"<<hap_pos<<" gt1 "<<gt1<<" gt2 "<<gt2<<" gt3 "<<gt3;
                    /*if ((chr2 > chr1) ||
                      ((chr2 == chr1) && ((hap_pos - dnm_pos) > window))) // the GTs are sorted by position, helps avoid iterating through entire file
                      break;*/

                    long dist = hap_pos - dnm_pos; // distance b/w DNM and phasing site
                    if(dist < 0) {
                        dist = -dist;
                    }
                    //cout<<endl<<"dist "<<dist<<"window "<<window;
                    if((chr2 == chr1) && (dist <= window) && (hap_pos != dnm_pos)) {
                        //char temp_file[20];
                        //getReadsFromBAM(bam_f, chr1, dnm_pos, hap_pos,temp_file); // get reads corresponding to the positions
                        returnsum += processReads(bam_f, chr1, dnm_pos, hap_pos, gt1, gt2,
                                                  variant_base);
                        //remove(temp_file);
                    }

                }
                fin2.close();
                if(returnsum == 0) {
                    cout << " - Insufficient reads present to phase this site.";
                }
                //cout<<"\nThe number of lines read GT file is "<<line_n2<<" DNM is "<<line_n1;
            } else {
                cout << "\nUnable to open parent GT file: " << parentGT_f << " ! Exiting!\n";
                exit(1);
            }
            fin1.getline(token1, 20, '\t');
        }
        fin1.close();
        cout << "\nThe number of lines read DNM file is " << line_n1;
    } else {
        cout << "\nUnable to open DNM file: " << DNM_f << " ! Exiting!\n";
        exit(1);
    }
    cout << "\n";
    exit(0);
}


typedef dng::CommandLineApp<dng::task::Phaser> App;

class PhaserApp : App {
public:

    PhaserApp(int argc, char *argv[]) : App(argc, argv) {}

    int operator()() {
        using namespace std;
        if(arg.version) {
            return CmdVersion();
        }
        if(arg.help) {
            return CmdHelp();
        }
        return task_(arg);
    }

protected:
    int CmdHelp() const {
        string usage_name(arg.run_name);
        if(usage_name.substr(0, 4) == "dng-") {
            usage_name[3] = ' ';
        }
        cerr << "Usage:\n"
             << "  dng phaser [options]\n\n"
             << "Options:\n"
             << "  --dnm filename: File with list of DNMs to be phased\n"
             << "  --pgt filename: File with parental genotypes for phasing sites\n"
             << "  --bam filename: bam file\n"
             << "  --window size: size of window for phasing sites (> insert size)\n"
             << endl;
        return EXIT_SUCCESS;
    }
};


int main(int argc, char *argv[]) {
    try {
        return PhaserApp(argc, argv)();
    } catch(std::exception &e) {
        std::cerr << e.what() << std::endl;
    }

    return EXIT_FAILURE;
}
