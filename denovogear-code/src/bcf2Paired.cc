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


void writeToSNPObject(pair_t* tumor, bcf1_t *b, bcf_hdr_t* h, int* g, int d, int mq, int& flag, int i, int i0)
{
  strcpy( tumor->chr, h->ns[b->tid] );
  tumor->pos = b->pos+1;
  tumor->ref_base = *b->ref;
  strcpy(tumor->alt, b->alt);
  tumor->depth = d;
  tumor->rms_mapQ = mq;
  strcpy( tumor->id, h->sns[i] );
  for (int l1 = 0; l1 < b->n_gi; ++l1) { //CHECK IF PER SAMPLE DEPTH AVAILABLE
    if (b->gi[l1].fmt == bcf_str2int("DP", 2)) {
      tumor->depth = ((uint16_t*)b->gi[l1].data)[i];
      //printf("\ndepth1 %d depth2: %d depth3: %d length: %d\n",((uint16_t*)b->gi[l1].data)[0], ((uint16_t*)b->gi[l1].data)[1], ((uint16_t*)b->gi[l1].data)[2], b->gi[l1].len);
    }
  }
  for (int j = 0; j < 10; ++j)
    tumor->lk[j] = g[j];
  if (tumor->rms_mapQ < MIN_MAPQ || tumor->depth < MIN_READ_DEPTH) 
    flag =1;
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

// Convert BCF to PairedSample - for each line iterate through samples and look for particular pair
int bcf2Paired (bcf_hdr_t *h, bcf1_t *b, Pair pair1, pair_t* tumor, pair_t* normal, int& flag)
{
	int a[4], k, g[10], l, map[4], k1, l1, j, i, i0, anno[16], dp, mq, d_rest, is_indel = 0;
	int found_pair = 2;// found_pair becomes zero when both samples are found.
	char *s;
	if (bcf_is_indel(b)) { 
		is_indel = 1;
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
				
		//found Tumor
		if( strcmp( pair1.tumorID, h->sns[i] ) == 0 ) {
			found_pair--;
			if( is_indel == 0 ) { // Write to Moms SNP object
			  writeToSNPObject(tumor, b, h, g, d, mq, flag, i, i0);
			}
		} 
		
		//found Normal
		if( strcmp( pair1.normalID, h->sns[i] ) == 0 ) {
			found_pair--;
			if( is_indel == 0 ) { // Write to Moms SNP object
			  writeToSNPObject(normal, b, h, g, d, mq, flag, i, i0);
			}
		}
				
	}
	
	//found entire pair, return
	if ( found_pair == 0 ) { 
	  return is_indel;
	} else {
	  printf("\n\nUnable to find pair, exiting Denovogear! ( %d, %d) ", found_pair, i);
	  return -3; // missing member	
	}
}
