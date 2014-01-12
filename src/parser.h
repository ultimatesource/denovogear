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

#ifndef PARSER_H_
#define PARSER_H_

#include "bcf.h"

#ifndef PED_PARSER_H
#define PED_PARSER_H
#include "pedParser.h"
#endif

#define MAX_QCALL_LINE 2048

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
} qcall_t;

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

// Store Paired Samples
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
} pair_t;

int bcf_2qcall(bcf_hdr_t *h, bcf1_t *b, Trio t, qcall_t* mom,
               qcall_t* dad, qcall_t* child, indel_t* mom_indel,
               indel_t* dad_indel, indel_t* child_indel, int& flag); // BCF to QCALL

int bcf2Paired (bcf_hdr_t *h, bcf1_t *b, Pair p, pair_t* tumor, pair_t* normal, int& flag);


#endif
