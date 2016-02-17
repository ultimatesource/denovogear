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

#include <string>
#include <vector>
#include <array>
#include <math.h>

#include "parser.h"
#include "htslib/vcf.h"



void writeToSNPObject(snp_object_t *mom_snp, const bcf_hdr_t *hdr, bcf1_t *rec,
                      int *g, int d,
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
        if(a != 1) {
            alt_str += ",";
        }
        alt_str += std::string(alleles[a]);
    }
    strcpy(mom_snp->alt, alt_str.c_str()); // ALT

    mom_snp->rms_mapQ = mq;
    strcpy(mom_snp->id, hdr->samples[i]); // TODO: Just pass in from calling fnc

    // Get "DP" field for each sample if it exits. Otherwise use depth estimated
    // from I16 fields 9 and 11
    int *res_array = NULL;
    int n_res_array = 0;
    int n_res = bcf_get_format_int32(hdr, rec, "DP", &res_array, &n_res_array);
    if(n_res == 3) {
        mom_snp->depth = res_array[i];
    } else {
        mom_snp->depth = d;
    }


    // Get PL liklihoods
    for(int j = 0; j < 10; j++) {
        mom_snp->lk[j] = g[j];
    }

    flag = mom_snp->rms_mapQ < MIN_MAPQ || mom_snp->depth < MIN_READ_DEPTH;

}

// TODO: Looking at the code the only diff I can find between writeToIndelObject and
// writeToSNPObject is the size of lk[] and ref_base (char vs char*). Should consider
// merging the two functions.
void writeToIndelObject(indel_t *mom_indel, const bcf_hdr_t *hdr, bcf1_t *rec,
                        int *g,
                        int d, int mq, int &flag, int i, int i0, std::vector<uint32_t> &pl_fields) {

    strcpy(mom_indel->chr, bcf_hdr_id2name(hdr, rec->rid)); // copy chrom
    mom_indel->pos = rec->pos + 1; // vcf posistion is stored in 0 based
    char **alleles = rec->d.allele;
    uint32_t n_alleles = rec->n_allele;
    strcpy(mom_indel->ref_base, alleles[0]); // REF
    std::string alt_str;
    for(int a = 1; a < n_alleles; a++) {
        if(a != 1) {
            alt_str += ",";
        }
        alt_str += std::string(alleles[a]);
    }
    strcpy(mom_indel->alt, alt_str.c_str()); // ALT
    mom_indel->rms_mapQ = mq;
    strcpy(mom_indel->id, hdr->samples[i]);

    int *res_array = NULL;
    int n_res_array = 0;
    int n_res = bcf_get_format_int32(hdr, rec, "DP", &res_array, &n_res_array);
    if(n_res == 3) {
        mom_indel->depth = res_array[i];
    } else {
        mom_indel->depth = d;
    }

    // Get PL liklihoods
    for(int j = 0; j < 3; j++) {
        mom_indel->lk[j] = pl_fields[j];
    }

    flag = mom_indel->rms_mapQ < MIN_MAPQ || mom_indel->depth < MIN_READ_DEPTH;

}


// Convert BCF to Qcall format - for each line iterate through samples and look for particular trio
// Parsing code adapted from BCFtools code
// TODO: Merge function with bcf2Paired.bcf2Paired()
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
        //return 10;
    }

    // Make sure the PL fields (phred-scaled genotype likihoods) exists, then
    // store in pl_fields
    sample_vals_int pl_fields;
    int *pl_array = NULL;
    int n_pl_array = 0;
    int n_pl = bcf_get_format_int32(hdr, rec, "PL", &pl_array, &n_pl_array);
    if(n_pl == 0) {
        return -6;
    } else {
        pl_fields.resize(n_samples);
        int sample_len = n_pl / n_samples;
        for(int a = 0; a < n_pl; a++) {
            // htslib will return array of size (Num Samples)x(PL size)
            int sample_index = a / sample_len;
            pl_fields[sample_index].push_back(pl_array[a]);
        }
    }

    // get I16 values from INFO field
    std::array<int, 16> anno;
    if(read_I16(rec, hdr, anno) != 0) {
        d_rest = 0;
    } else {
        d_rest = dp = anno[0] + anno[1] + anno[2] + anno[3];
    }


    // Calculate map quality from I16 fields 9 and 11
    mq = (int)(sqrt((double)(anno[9] + anno[11]) / dp) + .499);

    // Check that REF is a valid base
    a[0] = nt4_table[(int)alleles[0][0]];
    if(a[0] > 3) {
        return 10;
    }

    // Check that ALT alleles exists
    if(rec->n_allele < 2) {
        return -11;
    }


    // Map the alternative alleles
    int s;
    a[1] = a[2] = a[3] = -2; // -1 has a special meaning
    map[0] = map[1] = map[2] = map[3] = -2;
    map[a[0]] = 0;
    for(k = 0, s = 1, k1 = -1; k < 3 && s < rec->n_allele; ++k, s++) {
        // skip sites with multiple indel allele
        if(strlen(alleles[s]) > 1) {
            //return 10;
        }

        a[k + 1] = nt4_table[(int)alleles[s][0]];
        if(a[k + 1] >= 0) {
            map[a[k + 1]] = k + 1;
        } else {
            k1 = k + 1;
        }
    }


    for(k = 0; k < 4; ++k)
        if(map[k] < 0) { map[k] = k1; }

    std::vector<std::string> sample_ids;
    for(int a = 0; a < bcf_hdr_nsamples(hdr); a++) {
        sample_ids.push_back(std::string(hdr->samples[a]));
    }

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
                if(x < 0 || y < 0) {
                	// If PL field is not specified for a given genotype, just assume its likelihood is a close to 0 as possible.
                	g[j++] = MAX_PL;
                }
                else {
                	if(x > y) {
                		t = x;
                		x = y;
                		y = t;
                	}
                	// see VCF specifications, 'GL' format section
                	g[j++] = pl_fields[i][y * (y + 1) / 2 + x];
                }
            }
        }

        //found Mom
        if(strcmp(t.mID, sample_ids[i].c_str()) == 0) {
            found_trio--;
            if(indel == 0) {
                writeToSNPObject(mom_snp, hdr, rec, g, d, mq, flag, i, i0);
            } else {
                writeToIndelObject(mom_indel, hdr, rec, g, d, mq, flag, i, i0, pl_fields[i]);
            }
        }

        //found Dad
        if(strcmp(t.dID, sample_ids[i].c_str()) == 0) {
            found_trio--;
            if(indel == 0) {
                writeToSNPObject(dad_snp, hdr, rec, g, d, mq, flag, i, i0);
            } else {
                writeToIndelObject(dad_indel, hdr, rec, g, d, mq, flag, i, i0, pl_fields[i]);
            }
        }

        //found Child
        if(strcmp(t.cID, sample_ids[i].c_str()) == 0) {
            found_trio--;
            if(indel == 0) {
                writeToSNPObject(child_snp, hdr, rec, g, d, mq, flag, i, i0);
            } else {
                writeToIndelObject(child_indel, hdr, rec,  g, d, mq, flag, i, i0, pl_fields[i]);
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
