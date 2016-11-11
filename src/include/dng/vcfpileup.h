/*
 * Copyright (c) 2014 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
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

#pragma once
#ifndef DNG_VCFPILEUP_H
#define DNG_VCFPILEUP_H

#include <vector>
#include <limits>
#include <unordered_map>

#include <dng/hts/bam.h>
#include <htslib/faidx.h>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/metafunctions.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/algorithm/sort.hpp>

#include <dng/pool.h>
#include <dng/fileio.h>
#include <dng/cigar.h>
#include <dng/read_group.h>

#include "htslib/synced_bcf_reader.h"

namespace dng {
namespace pileup {
namespace vcf {


class VCFPileup {
public:

	typedef std::vector<int32_t> depth_list;

    typedef std::vector<depth_list> data_type;

    typedef std::vector<char> allele_list;

	//typedef void (callback_type)(bcf_hdr_t *, bcf1_t *);
	typedef void (callback_type)(const data_type &data, const allele_list &alleles, const char *chrom, size_t loc);


    template<typename LB>
    VCFPileup(const LB &lb) : libraries_{lb} {

    }

    template<typename Func>
    void operator()(const char *fname, Func func);

    template<typename Func>
    void operator()(const char *fname, bcf_srs_t *rec_reader, Func func);

private:
    ReadGroups::StrSet libraries_;
};

template<typename Func>
void VCFPileup::operator()(const char *fname, Func func) {


}


template<typename Func>
void VCFPileup::operator()(const char *fname, bcf_srs_t *rec_reader, Func func) {
    // TODO? If using multiple input vcf files then we may require scanners to search each file for the same position

    // type erase callback function
    std::function<callback_type> call_back(func);

    // Open the VCF/BCF file
    htsFile *fp = hts_open(fname, "r");
    // Get the header (should be only one input file)
    bcf_hdr_t *hdr = bcf_sr_get_header(rec_reader, 0);


    std::string samples;
    for(auto && str : libraries_) {
        samples += str + ',';
    }
    samples.pop_back();

    bcf_hdr_set_samples(hdr, samples.c_str(), 0);
    int variant_types;

    while(bcf_sr_next_line(rec_reader)) {
    	bcf1_t *rec = bcf_sr_get_line(rec_reader, 0);
    	variant_types = bcf_get_variant_types(rec);
    	// check that the current record is for an SNP or a REF and not an Indel, MNP, or something else
    	if((variant_types != VCF_SNP) && (variant_types != VCF_REF)){
    		continue;
    	}

    	// Won't be able to access ref->d unless we unpack the record first
    	bcf_unpack(rec, BCF_UN_STR);



        // get chrom, position, ref from their fields
        //const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
        //int32_t position = rec->pos;
        uint32_t n_alleles = rec->n_allele;
        uint32_t n_samples = bcf_hdr_nsamples(hdr);
        const char ref_base = *(rec->d.allele[0]);

        const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
        int pos = rec->pos;

        // REF+ALT in the order they appear in the record
        allele_list alleles(n_alleles);
        for(int a = 0; a < n_alleles; ++a) {
        	alleles[a] = *(rec->d.allele[a]);
        }


        // Read all the Allele Depths for every sample into ad array
        int *ad = NULL;
        int n_ad = 0;
        int n_ad_array = 0;
        n_ad = bcf_get_format_int32(hdr, rec, "AD", &ad, &n_ad_array);
        data_type read_depths(n_samples, depth_list(n_alleles));
        for(size_t sample_ndx = 0; sample_ndx < n_samples; ++sample_ndx) {
        	for(size_t allele_ndx = 0; allele_ndx < n_alleles; ++allele_ndx) {
        		read_depths[sample_ndx][allele_ndx] = ad[n_alleles * sample_ndx + allele_ndx];
        	}
        }



        // TODO? Check the QUAL field or PL,PP genotype fields
        // execute func
        call_back(read_depths, alleles, chrom, pos);
    }
}

}
}
} //namespace dng::pileup::vcf

#endif //DNG_VCFPILEUP_H

