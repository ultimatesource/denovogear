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

#pragma once
#ifndef OLD_PARSER_H_
#define OLD_PARSER_H_

#include <vector>
#include "pedParser.h"
#include <dng/hts/bcf.h>
#include <array>
//#include "htslib/vcf.h"

// TODO: qcall_t, pair_t, and indel_t can be merged. Use struct inheritance.
// TODO: Create namespace

#define MIN_READ_DEPTH 10
#define MIN_READ_DEPTH_INDEL 10
#define MIN_MAPQ 40
#define MAX_QCALL_LINE 2048

// maximum phred-scaled likelihood, as close to zero as possible (8 bit)
#define MAX_PL 255

// Stores the list of sample values (PL) for each sample
typedef std::vector<std::vector<uint32_t>> sample_vals_int;

// Store SNP info
typedef struct { //New struct for Q calls
    char chr[50]; /* Chromosome Number. */
    long pos; /* position, the first base in a chromosome has offset zero. */
    char ref_base; /* Either A, C, G or T */
    char alt[20]; /* ALT string */
    int depth; /* number of mapped reads */
    int rms_mapQ; /* RMS mapping quality */
    int min_lk; /* minimum lk capped at 255 */
    int lk[10];   /* log likelihood ratio, capped at 255 */
    char id[ID_LENGTH]; /* string for sample ID */
} snp_object_t;


// Store Paired Samples
//typedef struct { //New struct for Q calls
//    char chr[50]; /* Chromosome Number. */
//    long pos; /* position, the first base in a chromosome has offset zero. */
//    char ref_base; /* Either A, C, G or T */
//    char alt[20]; /* ALT string */
//    int depth; /* number of mapped reads */
//    int rms_mapQ; /* RMS mapping quality */
//    int min_lk; /* minimum lk capped at 255 */
//    int lk[10];   /* log likelihood ratio, capped at 255 */
//    char id[ID_LENGTH]; /* string for sample ID */
//} pair_t;

typedef snp_object_t qcall_t;
typedef snp_object_t pair_t;

// Store INDEL
typedef struct { //New struct for Indels
    char chr[50]; /* chromosome Number. */
    long pos; /* position, the first base in a chromosome has offset zero. */
    char ref_base[10000]; /* Either A, C, G or T */
    char alt[10000]; /* ALT string */
    int depth; /* number of mapped reads */
    int rms_mapQ; /* RMS mapping quality */
    int min_lk; /* minimum lk capped at 255 */
    int lk[3];   /* log likelihood ratio, capped at 255, ORDER = ref/ref, ref/alt & alt/alt */
    char id[ID_LENGTH]; /* string for sample ID */
    int length; /* to be filled in by Ruth */
} indel_t;


// Need both N and X to represent all alleles (X for the old VCF version, after v4.1 use N)
// TODO: turn into switch statement - a table lookup for only 5 values is boosting speed and make the code less managable
static int8_t nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, -1, 4,
    4, 4, 4, 4,  3, 4, 4, 4, -1, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, -1, 4,
    4, 4, 4, 4,  3, 4, 4, 4, -1, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


// Check that "I16" value exists in INFO field
static int read_I16(bcf1_t *rec, const bcf_hdr_t *hdr, std::array<int, 16> &anno) {
    std::vector<int> i16vals;
    //vcf.GetInfoValues("I16", i16vals);
    int *dst = NULL;
    int n_dst = 0;
    int array_size = bcf_get_info_int32(hdr, rec, "I16", &dst, &n_dst);

    if(array_size == 0) {
        return -1;
    }
    if(array_size < 16) {
        return -2;
    }

    //TODO: Should reset array?
    for(int a = 0; a < array_size; a++) {
        anno[a] = dst[a];
    }
    return 0;

}



/**
 * Used by the [pair/snp/indel]Like.cc functions to make alt allele string printable.
 * Removes the N/X allele (X is the pre-BCFv2.2 version of IUPAC N).
 */
inline void formatAltAlleles(std::string &alt) {
    size_t start = alt.find(",X");
    if(start == std::string::npos) {
        start = alt.find(",N");
    }

    if(start != std::string::npos) {
        alt.replace(start, 2, "");
    }
}


void writeToSNPObject(snp_object_t *mom_snp, const bcf_hdr_t *hdr, bcf1_t *rec,
                      int *g, int d,
                      int mq, int &flag, int i, int i0);

void writeToIndelObject(indel_t *mom_indel, const bcf_hdr_t *hdr, bcf1_t *rec,
                        int *g,
                        int d, int mq, int &flag, int i, int i0);




int bcf_2qcall(const bcf_hdr_t *h, bcf1_t *rec, Trio t, qcall_t *mom_snp,
               qcall_t *dad, qcall_t *child, indel_t *mom_indel,
               indel_t *dad_indel, indel_t *child_indel, int &flag); // BCF to QCALL

int bcf2Paired(const bcf_hdr_t *h, bcf1_t *b, Pair p, pair_t *tumor,
               pair_t *normal,
               int &flag);



#endif
