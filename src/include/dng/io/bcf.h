/*
 * Copyright (c) 2014-2017 Reed A. Cartwright
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
#ifndef DNG_IO_BCF_H
#define DNG_IO_BCF_H

#include <dng/hts/bcf.h>

#include <dng/utility.h>
#include <dng/depths.h>

namespace dng {
namespace io {

class BcfPileup {
public:
    using data_type = dng::pileup::AlleleDepths;

    typedef void (callback_type)(const data_type &, utility::location_t);

    int AddFile(const char* filename);

    // template<typename R>
    // void SelectLibraries(R &range);

    // void ResetLibraries();

    // const libraries_t& libraries() const {
    //     return output_libraries_;
    // }
    // size_t num_libraries() const {
    //     return output_libraries_.names.size();
    // }

    template<typename CallBack>
    void operator()(CallBack func);

    const hts::bcf::SyncedReader& reader() const {
        return reader_;
    }

private:
    hts::bcf::SyncedReader reader_;
};

template<typename CallBack>
void BcfPileup::operator()(CallBack call_back) {
    assert(reader_.num_readers() == 1); // only support one reader

    using utility::location_t;
    using pileup::AlleleDepths;

    int num_samples = bcf_hdr_nsamples(reader_.header(0));

    data_type line(0,0,num_samples);
    std::string alleles;
    alleles.reserve(5);

    while(reader_.NextLine()) {
        bcf1_t *rec = reader_.GetLine(0);
        if(rec == nullptr) {
            continue;
        }
        int variant_types = bcf_get_variant_types(rec);

        // check that the current record is for an SNP or a REF and not an Indel, MNP, or something else
        if((variant_types != VCF_SNP) && (variant_types != VCF_REF)){
            continue;
        }

        // Construct location
        line.location(utility::make_location(rec->rid, rec->pos));

        // REF+ALT in the order they appear in the record
        bcf_unpack(rec, BCF_UN_STR);
        alleles.clear();
        int num_alleles = rec->n_allele;
        for(int a = 0; a < num_alleles; ++a) {
            alleles += rec->d.allele[a];
        }
        line.resize(AlleleDepths::MatchLabel(alleles));

        // Read all the Allele Depths for every sample into ad array
        int *ad = nullptr;
        int n_ad = 0;
        int n_ad_array = 0;
        n_ad = bcf_get_format_int32(reader_.header(0), rec, "AD", &ad, &n_ad_array);
        if(n_ad < 0) {
            // AD tag is missing, so we do nothing at this time
            // TODO: support using calculated genotype likelihoods
            continue;
        }
        if(line.num_nucleotides() == num_alleles) {
            assert(n_ad == line.data_size());
            for(int i=0,k=0;i<num_samples;++i) {
                for(int a=0;a<num_alleles;++a) {
                    line(i,a) = ad[k++];
                }
            }           
        } else {
            assert(n_ad-num_samples == line.data_size());
            for(int i=0,k=0;i<num_samples;++i) {
                k++; // skip reference AD since it is N
                for(int a=1;a<num_alleles;++a) {
                    line(i,a-1) = ad[k++];
                }
            }
        }
        // execute func
        call_back(line);
    }
}

int BcfPileup::AddFile(const char* filename) {
    assert(filename != nullptr);
    assert(reader_.num_readers() == 0); // only support one input file

    return reader_.AddReader(filename);
}

} //namespace io
} //namespace dng

#endif //DNG_IO_BCF_H
