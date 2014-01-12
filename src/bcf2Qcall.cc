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
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>

#include "parser.h"
#include "prob1.h"
#include "kstring.h"
#include "time.h"
#include "kseq.h"


#define MIN_READ_DEPTH 10
#define MIN_READ_DEPTH_INDEL 10
#define MIN_MAPQ 40


void writeToSNPObject(qcall_t* mom_snp, bcf1_t *b, bcf_hdr_t* h, int* g, int d, int mq, int& flag, int i, int i0)
{
    strcpy( mom_snp->chr, h->ns[b->tid] );
    mom_snp->pos = b->pos+1;
    mom_snp->ref_base = *b->ref;
    strcpy(mom_snp->alt, b->alt);
    mom_snp->depth = d;
    mom_snp->rms_mapQ = mq;
    strcpy( mom_snp->id, h->sns[i] );
    for (int l1 = 0; l1 < b->n_gi; ++l1) { //CHECK IF PER SAMPLE DEPTH AVAILABLE
      if (b->gi[l1].fmt == bcf_str2int("DP", 2)) {
	mom_snp->depth = ((uint16_t*)b->gi[l1].data)[i];
	//printf("\ndepth1 %d depth2: %d depth3: %d length: %d\n",((uint16_t*)b->gi[l1].data)[0], ((uint16_t*)b->gi[l1].data)[1], ((uint16_t*)b->gi[l1].data)[2], b->gi[l1].len);
      }
    }
    for (int j = 0; j < 10; ++j)
      mom_snp->lk[j] = g[j];
    if (mom_snp->rms_mapQ < MIN_MAPQ || mom_snp->depth < MIN_READ_DEPTH)
      flag =1;
    //printf("\nSNP mom position %d depth %d", b->pos+1, mom_snp->depth);
}


void writeToIndelObject(indel_t* mom_indel, bcf1_t *b, bcf_hdr_t* h, int* g, int d, int mq, int& flag, int i, int i0)
{
  uint8_t *likl = static_cast<uint8_t *>(static_cast<uint8_t *>(b->gi[i0].data) + i * b->gi[i0].len);
  strcpy( mom_indel->chr, h->ns[b->tid] );
  mom_indel->pos = b->pos+1;
  strcpy(mom_indel->ref_base, b->ref);
  strcpy(mom_indel->alt, b->alt);
  mom_indel->depth = d;
  mom_indel->rms_mapQ = mq;
  strcpy( mom_indel->id, h->sns[i] );
  //printf("\nsample: %s, pos %d ",h->sns[i], b->pos+1 );
  for (int l1 = 0; l1 < b->n_gi; ++l1) { //CHECK IF PER SAMPLE DEPTH AVAILABLE
    if (b->gi[l1].fmt == bcf_str2int("DP", 2)) {
      mom_indel->depth = ((uint16_t*)b->gi[l1].data)[i];
      //printf("\ndepth1 %d depth2: %d depth3: %d length: %d\n",((uint16_t*)b->gi[l1].data)[0], ((uint16_t*)b->gi[l1].data)[1], ((uint16_t*)b->gi[l1].data)[2], b->gi[l1].len);
    }
  }
  for (int j = 0; j < 3; ++j) // R/R, R/A and A/A
    mom_indel->lk[j] = likl[j];
  if (mom_indel->rms_mapQ < MIN_MAPQ || mom_indel->depth < MIN_READ_DEPTH)
    flag =1;
  //printf("\nINDEL mom position: %d id: %s depth: %d lik: %d, %d, %d", b->pos+1,  mom_indel->id, d,  likl[0], likl[1], likl[2]);
}



static int8_t nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4, -1, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
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

static int read_I16(bcf1_t *b, int anno[16])
{
	char *p;
	int i;
	if ((p = strstr(b->info, "I16=")) == 0) return -1;
	p += 4;
	for (i = 0; i < 16; ++i) {
		anno[i] = strtol(p, &p, 10);
		if (anno[i] == 0 && (errno == EINVAL || errno == ERANGE)) return -2;
		++p;
	}
	return 0;
}

// Convert BCF to Qcall format - for each line iterate through samples and look for particular trio
// Parsing code adapted from BCFtools code
int bcf_2qcall(bcf_hdr_t *h, bcf1_t *b, Trio t, qcall_t* mom_snp, qcall_t* dad_snp, qcall_t* child_snp, indel_t* mom_indel, indel_t* dad_indel, indel_t* child_indel, int& flag)
{
	int a[4], k, g[10], l, map[4], k1, l1, j, i, i0, anno[16], dp, mq, d_rest, indel = 0, found_trio = 3;
	char *s;
	if (bcf_is_indel(b)) {
		//printf("\nINDEL position %d", b->pos+1);
		indel = 1;
	}
	if (b->ref[1] != 0 || b->n_alleles > 4) { // ref is not a single base ***
		//return 10;
	}
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == bcf_str2int("PL", 2)) break;
	if (i == b->n_gi) return -6; // no PL ***
	if (read_I16(b, anno) != 0) {
	  //return -7; // no I16; FIXME: can be improved ***
	  d_rest = 0;
	}
	else {
	  d_rest = dp = anno[0] + anno[1] + anno[2] + anno[3];
	}
	mq = (int)(sqrt((double)(anno[9] + anno[11]) / dp) + .499);
	i0 = i;
	a[0] = nt4_table[(int)b->ref[0]];
	if (a[0] > 3) { // ref is not A/C/G/T ***
		return 10;
	}
	a[1] = a[2] = a[3] = -2; // -1 has a special meaning
	if (b->alt[0] == 0) return -11; // no alternate allele ***
	map[0] = map[1] = map[2] = map[3] = -2;
	map[a[0]] = 0;
	for (k = 0, s = b->alt, k1 = -1; k < 3 && *s; ++k, s += 2) {
		if (s[1] != ',' && s[1] != 0)  { // ALT is not single base ***
			//return 10;
		}
		a[k+1] = nt4_table[(int)*s];
		if (a[k+1] >= 0) map[a[k+1]] = k+1;
		else k1 = k+1;
		if (s[1] == 0) break;
	}
	for (k = 0; k < 4; ++k)
		if (map[k] < 0) map[k] = k1;


	for (i = 0; i < h->n_smpl; ++i) {	//Iterate through all samples
		int d;
		uint8_t *p = static_cast<uint8_t *>(static_cast<uint8_t *>(b->gi[i0].data) + i * b->gi[i0].len);
		for (j = 0; j < b->gi[i0].len; ++j)
			if (p[j]) break;
		d = (int)((double)d_rest / (h->n_smpl - i) + .499);
		if (d == 0) d = 1;
		if (j == b->gi[i0].len) d = 0;
		d_rest -= d;
		for (k = j = 0; k < 4; ++k) {
			for (l = k; l < 4; ++l) {
				int t, x = map[k], y = map[l];
				if (x > y) t = x, x = y, y = t; // swap
				g[j++] = p[y * (y+1) / 2 + x];
			}
		}

		//found Mom
		if(strcmp(t.mID, h->sns[i]) == 0) {
			found_trio--;
			if(indel == 0) {
			  writeToSNPObject(mom_snp, b, h, g, d, mq, flag, i, i0);
			}
			else {
			  writeToIndelObject(mom_indel, b, h, g, d, mq, flag, i, i0);
			}
		}

		//found Dad
		if(strcmp(t.dID, h->sns[i]) == 0) {
			found_trio--;
			if(indel == 0) {
			  writeToSNPObject(dad_snp, b, h, g, d, mq, flag, i, i0);
			}
			else {
			  writeToIndelObject(dad_indel, b, h, g, d, mq, flag, i, i0);
			}
		}

		//found Child
		if( strcmp(t.cID, h->sns[i]) == 0) {
			found_trio--;
			if( indel == 0 ) {
			  writeToSNPObject(child_snp, b, h, g, d, mq, flag, i, i0);
			}
			else {
			  writeToIndelObject(child_indel, b, h, g, d, mq, flag, i, i0);
			}
		}
	}

	//found entire trio, return
	if ( found_trio == 0 ) {
		return indel;
	} else {
		printf("\n\nUnable to find trio. Code %d:%d ", found_trio, i);
		return -3; //missing member
	}

}
