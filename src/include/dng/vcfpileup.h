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

namespace dng { namespace pileup { namespace vcf {


class VCFPileup {
public:
        typedef NodeList list_type;
	typedef Node node_type;
	typedef detail::NodePool pool_type;

	typedef std::vector<list_type> data_type;
	typedef void (callback_type)(const data_type&, uint64_t);

	template<typename InFiles, typename Func>
	void operator()(InFiles &range, Func func);

	template<typename RG>
	MPileup(const RG& rg) : pool_(4048), read_groups_(rg) {
	        boost::sort(read_groups_);
	}

	

private:
	char *fname = 0;
	bool iscompressed = false;
};

template<typename InFile, typename Func>
void VCFPileup::operator()(InFile fname, Func func) {
        // TODO: If using multiple input vcf files then we may require scanners to search each file for the same position

        // type erase callback function
        function<callback_type> call_back(func);
   
        // Open the VCF/BCF file
        htsFile *fp = hts_open(fname, "r");
	bcf_hdr_t *hdr = bcf_hdr_read(fp);
	bcf1_t *rec = bcf_init1();

	while(bcf_read1(fp, hdr, rec) >= 0) {
	  // check that the current record is for an SNP and not an Indel, MNP, or something else
	  if(bcf_get_variant_types(rec) != VCF_SNP)
	    continue;

	  // TODO: Check the QUAL field or PL,PP genotype fields

	  // execute func
	  call_back(rec, hdr);
	}
}
	  
}}} //namespace dng::pileup::vcf

#endif //DNG_VCFPILEUP_H

