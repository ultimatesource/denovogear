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
int bcf_2qcall(bcf_hdr_t *h, bcf1_t *b, Trio t, qcall_t* mom_snp, qcall_t* dad_snp, qcall_t* child_snp, indel_t* mom_indel, indel_t* dad_indel, indel_t* child_indel, int& flag)
{
	int a[4], k, g[10], l, map[4], k1, l1, j, i, i0, anno[16], dp, mq, d_rest, indel = 0, found_trio = 3;
	char *s;
	if (bcf_is_indel(b)) { 
		//printf("\nINDEL position %d", b->pos+1);
		indel = 1;
	}
	if (b->ref[1] != 0 || b->n_alleles > 4) { // ref is not a single base ***
		//printf("\nposition %d ref not a single base", b->pos+1);
		//indel = 1;
		//return 10;
	} 
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == bcf_str2int("PL", 2)) break;
	if (i == b->n_gi) return -6; // no PL ***
	if (read_I16(b, anno) != 0) return -7; // no I16; FIXME: can be improved ***
	d_rest = dp = anno[0] + anno[1] + anno[2] + anno[3];
	//if (dp == 0) return -8; // depth is zero ***
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
		if (s[1] != ',' && s[1] != 0)  { // ALT is not single base *** // LINES WITH BOTH SNP, INDEL ??
			//printf("\nposition %d alt not a single base", b->pos+1); 
			//indel =1 ;
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
		if( strcmp( t.mID, h->sns[i] ) == 0 ) {
			found_trio--;
			
			if( indel == 0 ) { // Write to Moms SNP object
				strcpy( mom_snp->chr, h->ns[b->tid] );
				mom_snp->pos = b->pos+1;
				mom_snp->ref_base = *b->ref;
				strcpy(mom_snp->alt, b->alt);
				mom_snp->depth = d;
				mom_snp->rms_mapQ = mq;
				strcpy( mom_snp->id, h->sns[i] );
				for (l1 = 0; l1 < b->n_gi; ++l1) { //CHECK IF PER SAMPLE DEPTH AVAILABLE
					if (b->gi[l1].fmt == bcf_str2int("DP", 2)) {
						mom_snp->depth = ((uint16_t*)b->gi[l1].data)[i];
						//printf("\ndepth1 %d depth2: %d depth3: %d length: %d\n",((uint16_t*)b->gi[l1].data)[0], ((uint16_t*)b->gi[l1].data)[1], ((uint16_t*)b->gi[l1].data)[2], b->gi[l1].len);
					}
				}
				for (j = 0; j < 10; ++j)
					mom_snp->lk[j] = g[j];
				if (mom_snp->rms_mapQ < MIN_MAPQ || mom_snp->depth < MIN_READ_DEPTH) 
					flag =1;
				//printf("\nSNP mom position %d depth %d", b->pos+1, mom_snp->depth); 
			}
			else {
				uint8_t *likl = static_cast<uint8_t *>(static_cast<uint8_t *>(b->gi[i0].data) + i * b->gi[i0].len);
				strcpy( mom_indel->chr, h->ns[b->tid] );
				mom_indel->pos = b->pos+1;
				strcpy(mom_indel->ref_base, b->ref);
				strcpy(mom_indel->alt, b->alt);
				mom_indel->depth = d;
				mom_indel->rms_mapQ = mq;
				strcpy( mom_indel->id, h->sns[i] );
				//printf("\nsample: %s, pos %d ",h->sns[i], b->pos+1 );
				for (l1 = 0; l1 < b->n_gi; ++l1) { //CHECK IF PER SAMPLE DEPTH AVAILABLE
					if (b->gi[l1].fmt == bcf_str2int("DP", 2)) {
						mom_indel->depth = ((uint16_t*)b->gi[l1].data)[i];
						//printf("\ndepth1 %d depth2: %d depth3: %d length: %d\n",((uint16_t*)b->gi[l1].data)[0], ((uint16_t*)b->gi[l1].data)[1], ((uint16_t*)b->gi[l1].data)[2], b->gi[l1].len);
					}
				}
				for (j = 0; j < 3; ++j) // R/R, R/A and A/A
					mom_indel->lk[j] = likl[j];
				if (mom_indel->rms_mapQ < MIN_MAPQ || mom_indel->depth < MIN_READ_DEPTH) 
					flag =1;
				//printf("\nINDEL mom position: %d id: %s depth: %d lik: %d, %d, %d", b->pos+1,  mom_indel->id, d,  likl[0], likl[1], likl[2]); 
			}
		} 
		
		//found Dad
		if( strcmp( t.dID, h->sns[i] ) == 0 ) {
			found_trio--;
			
			if( indel == 0 ) { // Write to Dads SNP object
				strcpy( dad_snp->chr, h->ns[b->tid] );
				dad_snp->pos = b->pos+1;
				dad_snp->ref_base = *b->ref;
				strcpy(dad_snp->alt, b->alt);
				dad_snp->depth = d;
				dad_snp->rms_mapQ = mq;
				strcpy( dad_snp->id, h->sns[i] );
				//printf("\nsample: %s, pos %d ",h->sns[i], b->pos+1 );
				for (l1 = 0; l1 < b->n_gi; ++l1) { //CHECK IF PER SAMPLE DEPTH AVAILABLE
					if (b->gi[l1].fmt == bcf_str2int("DP", 2)) {
						dad_snp->depth = ((uint16_t*)b->gi[l1].data)[i];
						//printf("\ndepth1 %d depth2: %d depth3: %d length: %d\n",((uint16_t*)b->gi[l1].data)[0], ((uint16_t*)b->gi[l1].data)[1], ((uint16_t*)b->gi[l1].data)[2], b->gi[l1].len);
					}
				}
				for (j = 0; j < 10; ++j)
					dad_snp->lk[j] = g[j];
				if (dad_snp->rms_mapQ < MIN_MAPQ || dad_snp->depth < MIN_READ_DEPTH) 
					flag =1;
				//printf("\nSNP dad position %d depth %d", b->pos+1, dad_snp->depth); 
			}
			else {
				uint8_t *likl = static_cast<uint8_t *>(static_cast<uint8_t *>(b->gi[i0].data) + i * b->gi[i0].len);
				strcpy( dad_indel->chr, h->ns[b->tid] );
				dad_indel->pos = b->pos+1;
				strcpy(dad_indel->ref_base, b->ref);
				strcpy(dad_indel->alt, b->alt);
				dad_indel->depth = d;
				dad_indel->rms_mapQ = mq;
				strcpy( dad_indel->id, h->sns[i] );
				//printf("\nsample: %s, pos %d ",h->sns[i], b->pos+1 );
				for (l1 = 0; l1 < b->n_gi; ++l1) { //CHECK IF PER SAMPLE DEPTH AVAILABLE
					if (b->gi[l1].fmt == bcf_str2int("DP", 2)) {
						dad_indel->depth = ((uint16_t*)b->gi[l1].data)[i];
						//printf("\ndepth1 %d depth2: %d depth3: %d length: %d\n",((uint16_t*)b->gi[l1].data)[0], ((uint16_t*)b->gi[l1].data)[1], ((uint16_t*)b->gi[l1].data)[2], b->gi[l1].len);
					}
				}
				for (j = 0; j < 3; ++j) // R/R, R/A and A/A
					dad_indel->lk[j] = likl[j];
				if (dad_indel->rms_mapQ < MIN_MAPQ || dad_indel->depth < MIN_READ_DEPTH) 
					flag =1;
				//printf("\nINDEL dad position: %d id: %s depth: %d lik: %d, %d, %d", b->pos+1,  dad_indel->id, d,  likl[0], likl[1], likl[2]); 
			}
			
		}
		
		//found Child
		if (strcmp( t.cID, h->sns[i]) == 0 ) {
			found_trio--;
			
			if( indel == 0 ) { // Write to Childs SNP object
				strcpy( child_snp->chr, h->ns[b->tid] );
				child_snp->pos = b->pos+1;
				child_snp->ref_base = *b->ref;
				strcpy(child_snp->alt, b->alt);
				child_snp->depth = d;
				child_snp->rms_mapQ = mq;
				strcpy( child_snp->id, h->sns[i] );
				//printf("\nsample: %s, pos %d ",h->sns[i], b->pos+1 );
				for (l1 = 0; l1 < b->n_gi; ++l1) { //CHECK IF PER SAMPLE DEPTH AVAILABLE
					if (b->gi[l1].fmt == bcf_str2int("DP", 2)) {
						child_snp->depth = ((uint16_t*)b->gi[l1].data)[i];
						//printf("\ndepth1 %d depth2: %d depth3: %d length: %d\n",((uint16_t*)b->gi[l1].data)[0], ((uint16_t*)b->gi[l1].data)[1], ((uint16_t*)b->gi[l1].data)[2], b->gi[l1].len);
					}
				}
				for (j = 0; j < 10; ++j) 
					child_snp->lk[j] = g[j];
				if (child_snp->rms_mapQ < MIN_MAPQ || child_snp->depth < MIN_READ_DEPTH) 
					flag =1;
				//printf("\nSNP child position %d depth %d", b->pos+1, child_snp->depth); 
			}
			else {
				uint8_t *likl = static_cast<uint8_t *>(static_cast<uint8_t *>(b->gi[i0].data) + i * b->gi[i0].len);
				strcpy( child_indel->chr, h->ns[b->tid] );
				child_indel->pos = b->pos+1;
				strcpy(child_indel->ref_base, b->ref);
				strcpy(child_indel->alt, b->alt);
				child_indel->depth = d;
				child_indel->rms_mapQ = mq;
				strcpy( child_indel->id, h->sns[i] );
				//printf("\nsample: %s, pos %d ",h->sns[i], b->pos+1 );
				for (l1 = 0; l1 < b->n_gi; ++l1) { //CHECK IF PER SAMPLE DEPTH AVAILABLE
					if (b->gi[l1].fmt == bcf_str2int("DP", 2)) {
						child_indel->depth = ((uint16_t*)b->gi[l1].data)[i];
						//printf("\ndepth1 %d depth2: %d depth3: %d length: %d\n",((uint16_t*)b->gi[l1].data)[0], ((uint16_t*)b->gi[l1].data)[1], ((uint16_t*)b->gi[l1].data)[2], b->gi[l1].len);
					}
				}
				for (j = 0; j < 3; ++j) // R/R, R/A and A/A
					child_indel->lk[j] = likl[j];
				//printf("\n\nINDEL child position: %d id: %s depth: %d lik: %d, %d, %d", b->pos+1, child_indel->id, d,  likl[0], likl[1], likl[2]); 
				//printf("\nref base:%s alt:%s quality:%d", child_indel->ref_base, child_indel->alt, child_indel->rms_mapQ);
			}
		}
		
		
	}
	
	//found entire trio, return
	if ( found_trio == 0 ) { 
		return indel;
	} else {
		printf("\n\nUnable to find trio, exiting Denovogear! ( %d, %d) ", found_trio, i);
		return -3; //missing member	
	}
	
}
