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


// NOTE: Both bcf2Paired.cc and bcf2Qcall.cc were using a lot of redundant code. To make code
//       managable I move much of into parser.h. 
// TODO: Run more tests with PAIR SNPs and Indels. If working then delete the commented sectionb below
// TODO: Most of the functionality in bcf2paired() and bcf_2qcall() can be merged together 

/**************************************************************************************
#include "prob1.h"
#include "kstring.h"
#include "time.h"
#include "kseq.h"


#define MIN_READ_DEPTH 10
#define MIN_READ_DEPTH_INDEL 10
#define MIN_MAPQ 40

void writeToSNPObject(pair_t *tumor, bcf1_t *rec, bcf_hdr_t *hdr, int *g, int d,
                      int mq, int &flag, int i, int i0) {
  
  strcpy(tumor->chr, bcf_hdr_id2name(hdr, rec->rid));
  tumor->pos = rec->pos + 1;
 
  char **alleles = rec->d.allele;
  uint32_t n_alleles = bcf_hdr_nsamples(hdr);
  tumor->ref_base = alleles[0][0];

  std::string alt_str;
  for(int a = 1; a < n_alleles; a++) {
    if(a != 1)
      alt_str += ",";
    alt_str += std::string(alleles[a]);
  }
  strcpy(tumor->alt, alt_str.c_str());

  tumor->rms_mapQ = mq;

  int *res_array = NULL;
  int n_res_array = 0;
  int n_res = bcf_get_format_int32(hdr, rec, "DP", &res_array, &n_res_array);
  if(n_res == 3) {
    tumor->depth = res_array[i];
  }
  else {
    tumor->depth = d;
  }

  for(int j = 0; j < 10; j++)
    tumor->lk[j] = g[j];

  flag = tumor->rms_mapQ < MIN_MAPQ || tumor->depth < MIN_READ_DEPTH;

}


static int8_t nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
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

static int read_I16(bcf1_t *rec, bcf_hdr_t *hdr, std::vector<int> &anno) {

  // Check that "I16" value exists in INFO field
  std::vector<int> i16vals;
  int *dst = NULL;
  int n_dst = 0;
  int array_size = bcf_get_info_int32(hdr, rec, "I16", &dst, &n_dst);

  if(array_size == 0)
    return -1;
  if(array_size < 16)
    return -2;

  //TODO: Should reset array
  for(int a = 0; a < array_size; a++)
    anno.push_back(dst[a]);      
  return 0;

}
**************************************************************************************/

// TODO: Merge first part of code with bcf2QCall.bcf_2qcall()
// Convert BCF to PairedSample - for each line iterate through samples and look for particular pair
int bcf2Paired(const bcf_hdr_t *hdr, bcf1_t *rec, Pair pair1, pair_t *tumor,
               pair_t *normal, int &flag) {
  int a[4], k, g[10], l, map[4], k1, l1, j, i, i0, /*anno[16],*/ dp, mq, d_rest,
        is_indel = 0;
    int found_pair = 2;// found_pair becomes zero when both samples are found.
    //char *s;


    // Check if record is INDEL
    is_indel = bcf_get_variant_types(rec) == hts::bcf::File::INDEL;
    

    char **alleles = rec->d.allele;
    uint32_t n_alleles = rec->n_allele;
    //TODO: conditional was empty, should we remove or implement?
    /*
    if(b->ref[1] != 0 || b->n_alleles > 4) {  // ref is not a single base ***
        //printf("\nposition %d ref not a single base", b->pos+1);
        //indel = 1;
        //return 10;
    }
    */

    // Make sure the PL fields (phred-scaled genotype likihoods) exists
    sample_vals_int pl_fields;
    int *res_array = NULL;
    int n_res_array = 0;
    int n_res = bcf_get_format_int32(hdr, rec, "PL", &res_array, &n_res_array);
    if(n_res == 0)
      return -6;
    else {
      pl_fields.resize(n_alleles);
      int sample_len = n_res/bcf_hdr_nsamples(hdr);
      for(int a = 0; a < n_res; a++) {
	int sample_index = a/sample_len;
	pl_fields[sample_index].push_back(res_array[a]);
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

    // Calculate map quality
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

      //TODO: conditional block mpty in original code, implement or remove
      /*
      // skip sites with multiple indel allele
      if(strlen(alleles[s]) > 1) {
	indel = 1;
	return 10;
      }
      */
      
      // BUG: Implement conditional above or allow sites with multiple indel alleles
      a[k + 1] = nt4_table[(int)alleles[s][0]];
      if(a[k + 1] >= 0) { 
	map[a[k + 1]] = k + 1; 
      }
      else { 
	k1 = k + 1; 
      }
    }

    for(k = 0; k < 4; ++k)
      if(map[k] < 0) { map[k] = k1; }

    std::vector<std::string> sample_ids;
    for(int a = 0; a < bcf_hdr_nsamples(hdr); a++)
      sample_ids.push_back(std::string(hdr->samples[a]));

    // Iterate through each sample, and for each sample map the values
    // in the PL field into g[], and estimate the depth    
    for(i = 0; i < bcf_hdr_nsamples(hdr); i++) {
      // Go to the first non-zero value in the PL field
      for(j = 0; j < pl_fields[i].size() && pl_fields[i][j]; j++);
      
      // Estimate the depth using I16 fields 1 to 4, divided by num of samples
      int d = (int)((double)d_rest / (bcf_hdr_nsamples(hdr) - i) + .4999);
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
	  g[j++] = pl_fields[i][y*(y+1)/2 + x];
	}
      }
    
      //found Tumor
      if(strcmp(pair1.tumorID, sample_ids[i].c_str()) == 0) {
	found_pair--;
	if(is_indel == 0) {   // Write to Moms SNP object
	  writeToSNPObject(tumor, hdr, rec, g, d, mq, flag, i, i0);
	}
      }

        //found Normal
      if(strcmp(pair1.normalID, sample_ids[i].c_str()) == 0) {
	found_pair--;
	if(is_indel == 0) {   // Write to Moms SNP object
	  writeToSNPObject(normal, hdr, rec, g, d, mq, flag, i, i0);
	}
      }

    }
    

    //found entire pair, return
    if(found_pair == 0) {
        return is_indel;
    } else {
        printf("\n\nUnable to find pair, exiting Denovogear! ( %d, %d) ", found_pair,
               i);
        return -3; // missing member
    }
}
