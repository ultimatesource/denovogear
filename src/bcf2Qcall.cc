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
/*
#include "prob1.h"
#include "kstring.h"
#include "time.h"
#include "kseq.h"
*/

#include "htslib/vcf.h"


#define MIN_READ_DEPTH 10
#define MIN_READ_DEPTH_INDEL 10
#define MIN_MAPQ 40

typedef std::vector<std::vector<int>> sample_vals_int;

void writeToSNPObject(qcall_t *mom_snp, const bcf_hdr_t *hdr, bcf1_t *rec, int *g, int d,
                      int mq, int &flag, int i, int i0) {

  strcpy(mom_snp->chr, bcf_hdr_id2name(hdr, rec->rid)); // copy chrom
  mom_snp->pos = rec->pos + 1; // vcf posistion is stored in 0 based

  // Get the ref + alt alleles
  // TODO: just pass in alleles from bcf_2qcall
  char **alleles = rec->d.allele;
  uint32_t n_alleles = rec->n_allele;

  mom_snp->ref_base = alleles[0][0]; // REF
  std::string alt_str;
  for(int a = 1; a < n_alleles; a++) {
    // convert alt alleles into comma seperated string
    // TODO: htslib may already have a function to do this, investigate
    if(a != 1) 
      alt_str += ",";
    alt_str += std::string(alleles[a]);
  }
  strcpy(mom_snp->alt, alt_str.c_str()); // ALT

  mom_snp->rms_mapQ = mq;

  // Get "DP" field for each sample if it exits. Otherwise use depth estimated
  // from I16 fields 9 and 11
  int *res_array = NULL;
  int n_res_array = 0;
  int n_res = bcf_get_format_int32(hdr, rec, "DP", &res_array, &n_res_array);
  if(n_res == 3) {
    mom_snp->depth = res_array[i];
  }
  else {
    mom_snp->depth = d;
  }

  // Get PL liklihoods
  for(int j = 0; j < 10; j++)
    mom_snp->lk[j] = g[j];

  flag = mom_snp->rms_mapQ < MIN_MAPQ || mom_snp->depth < MIN_READ_DEPTH;
    
}


// TODO: Is there any difference between this and the function above? Just merge
void writeToIndelObject(indel_t *mom_indel, const bcf_hdr_t *hdr, bcf1_t *rec, int *g,
                        int d, int mq, int &flag, int i, int i0) {

  /*
    uint8_t *likl = static_cast<uint8_t *>(static_cast<uint8_t *>
                                           (b->gi[i0].data) + i * b->gi[i0].len);
    strcpy(mom_indel->chr, h->ns[b->tid]);
    mom_indel->pos = b->pos + 1;
    strcpy(mom_indel->ref_base, b->ref);
    strcpy(mom_indel->alt, b->alt);
    mom_indel->depth = d;
    mom_indel->rms_mapQ = mq;
    strcpy(mom_indel->id, h->sns[i]);
    //printf("\nsample: %s, pos %d ",h->sns[i], b->pos+1 );
    for(int l1 = 0; l1 < b->n_gi; ++l1) {  //CHECK IF PER SAMPLE DEPTH AVAILABLE
        if(b->gi[l1].fmt == bcf_str2int("DP", 2)) {
            mom_indel->depth = ((uint16_t *)b->gi[l1].data)[i];
            //printf("\ndepth1 %d depth2: %d depth3: %d length: %d\n",((uint16_t*)b->gi[l1].data)[0], ((uint16_t*)b->gi[l1].data)[1], ((uint16_t*)b->gi[l1].data)[2], b->gi[l1].len);
        }
    }
    for(int j = 0; j < 3; ++j) { // R/R, R/A and A/A
        mom_indel->lk[j] = likl[j];
    }
    if(mom_indel->rms_mapQ < MIN_MAPQ || mom_indel->depth < MIN_READ_DEPTH) {
        flag = 1;
    }
    //printf("\nINDEL mom position: %d id: %s depth: %d lik: %d, %d, %d", b->pos+1,  mom_indel->id, d,  likl[0], likl[1], likl[2]);
    */
}


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
static int read_I16(bcf1_t *rec, const bcf_hdr_t *hdr, std::vector<int> &anno)
{
  std::vector<int> i16vals;
  //vcf.GetInfoValues("I16", i16vals);
  int *dst = NULL;
  int n_dst = 0;
  int array_size = bcf_get_info_int32(hdr, rec, "I16", &dst, &n_dst);

  if(array_size == 0)
    return -1;
  if(array_size < 16)
    return -2;

  //TODO: Should reset array?
  for(int a = 0; a < array_size; a++)
    anno.push_back(dst[a]);      
  return 0;

}


// Convert BCF to Qcall format - for each line iterate through samples and look for particular trio
// Parsing code adapted from BCFtools code
int bcf_2qcall(const bcf_hdr_t *hdr, bcf1_t *rec, Trio t, qcall_t *mom_snp,
               qcall_t *dad_snp, qcall_t *child_snp, indel_t *mom_indel, indel_t *dad_indel,
               indel_t *child_indel, int &flag) {
  
  int a[4], k, g[10], l, map[4], k1, l1, j, i, i0, /*anno[16],*/ dp, mq, d_rest,
    /*indel = 0,*/ found_trio = 3;

  // Check if INDEL or SNP
  int indel = (bcf_get_variant_types(rec) == hts::bcf::File::INDEL);
  
  char **alleles = rec->d.allele;
  uint32_t n_alleles = rec->n_allele;
  uint32_t n_samples = bcf_hdr_nsamples(hdr);

  // Make sure reference and alt alleles are only single-bases
  if(strlen(alleles[0]) > 1 || n_alleles > 4) {
    return 10;
  }
  
  // Make sure the PL fields (phred-scaled genotype likihoods) exists, then
  // store in pl_fields
  sample_vals_int pl_fields;
  int *pl_array = NULL;
  int n_pl_array = 0;
  int n_pl = bcf_get_format_int32(hdr, rec, "PL", &pl_array, &n_pl_array);
  if(n_pl == 0) {
    return -6;
  }
  else {
    pl_fields.resize(n_samples);
    int sample_len = n_pl/n_samples;
    for(int a = 0; a < n_pl; a++) {
      // htslib will return array of size (Num Samples)x(PL size)
      int sample_index = a/sample_len;
      pl_fields[sample_index].push_back(pl_array[a]);
    }
  }

  // get I16 values from INFO field
  std::vector<int> anno;
  if(read_I16(rec, hdr, anno) != 0) {
    d_rest = 0;
  }
  else {
    d_rest = dp = anno[0] + anno[1] + anno[2] + anno[3];
  }

  // Calculate map quality from I16 fields 9 and 11
  mq = (int)(sqrt((double)(anno[9] + anno[11]) / dp) + .499);

  // Check that REF is a valid base
  a[0] = nt4_table[(int)alleles[0][0]];
  if(a[0] > 3)
    return 10;

  // Check that ALT alleles exists
  if(rec->n_allele < 2)
    return -11;

  // Map the alternative alleles
  int s;
  a[1] = a[2] = a[3] = -2; // -1 has a special meaning
  map[0] = map[1] = map[2] = map[3] = -2;
  map[a[0]] = 0;
  for(k = 0, s = 1, k1 = -1; k < 3 && s < rec->n_allele; ++k, s++) {
    // skip sites with multiple indel allele
    if(strlen(alleles[s]) > 1)
      return 10;
    
    a[k + 1] = nt4_table[(int)alleles[s][0]];
    if(a[k + 1] >= 0) { 
      map[a[k + 1]] = k + 1; 
    }
    else { 
      k1 = k + 1; 
    }
  }

  // TODO: This doesn't handle the case when there is no 'N' alternative allele
  for(k = 0; k < 4; ++k)
    if(map[k] < 0) { map[k] = k1; }
  
  std::vector<std::string> sample_ids;
  for(int a = 0; a < bcf_hdr_nsamples(hdr); a++)
    sample_ids.push_back(std::string(hdr->samples[a]));
  
  // Iterate through each sample, and for each sample map the values
  // in the PL field into g[], and estimate the depth    
  for(i = 0; i < n_samples; i++) {
    // Go to the first non-zero value in the PL field
    for(j = 0; j < pl_fields[i].size() && pl_fields[i][j]; j++);
    
    // Estimate the depth using I16 fields 1 to 4, divided by num of samples
    int d = (int)((double)d_rest / (n_samples - i) + .4999);
    if(d == 0) {
      d = 1;
    }
    if(j == pl_fields[i].size()) {
      d = 0;
    }
    d_rest -= d;
    
    for(k = j = 0; k < 4; k++) {
      for(l = k; l < 4; l++) { //AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
	int t, x = map[k], y = map[l];
	if(x > y) {
	  t = x;
	  x = y;
	  y = t;
	}
	// see VCF specifications, 'GL' format section
	g[j++] = pl_fields[i][y*(y+1)/2 + x];
      }
    }

    //found Mom
    if(strcmp(t.mID, sample_ids[i].c_str()) == 0) {
      found_trio--;
      if(indel == 0) {
	writeToSNPObject(mom_snp, hdr, rec, g, d, mq, flag, i, i0);
      } else {
	writeToIndelObject(mom_indel, hdr, rec, g, d, mq, flag, i, i0);
      }
    }
  
    //found Dad
    if(strcmp(t.dID, sample_ids[i].c_str()) == 0) {
      found_trio--;
      if(indel == 0) {
	writeToSNPObject(dad_snp, hdr, rec, g, d, mq, flag, i, i0);
      } else {
	writeToIndelObject(dad_indel, hdr, rec, g, d, mq, flag, i, i0);
      }
    }
      
    //found Child
    if(strcmp(t.cID, sample_ids[i].c_str()) == 0) {
      found_trio--;
      if(indel == 0) {
	writeToSNPObject(child_snp, hdr, rec, g, d, mq, flag, i, i0);
      } else {
	writeToIndelObject(child_indel, hdr, rec,  g, d, mq, flag, i, i0);
      }
    }
  }
  

  //found entire trio, return
  if(found_trio == 0) {
    return indel;
  } else {
    printf("\n\nUnable to find trio. Code %d:%d ", found_trio, i);
    return -3; //missing member
  }
  

}
