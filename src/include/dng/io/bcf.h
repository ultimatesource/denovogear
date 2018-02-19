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

#include <cstring>

#include <dng/hts/bcf.h>

#include <dng/utility.h>
#include <dng/depths.h>
#include <dng/library.h>
#include <dng/regions.h>

#include <boost/spirit/include/karma_generate.hpp>


namespace dng {
namespace io {

class BcfPileup {
public:
    using data_type = bcf1_t *;

    typedef void (callback_type)(const data_type &, utility::location_t);

    int AddFile(const char* filename);

    template<typename R>
    void SelectLibraries(R &range);

    void ResetLibraries();

    int SetRegions(const regions::contig_fragments_t &frags, bool one_indexed=true);

    const std::vector<regions::contig_t>& contigs() const {
        return contigs_;
    }

    const libraries_t& libraries() const {
         return output_libraries_;
    }
    size_t num_libraries() const {
         return output_libraries_.names.size();
    }

    template<typename CallBack>
    void operator()(CallBack func);

    const hts::bcf::SyncedReader& reader() const {
        return reader_;
    }

    template<typename A>
    static BcfPileup open_and_setup(const A& arg);

private:
    hts::bcf::SyncedReader reader_;

    void ParseSampleLabels(int index);
    void ParseContigs(int index);

    dng::libraries_t input_libraries_;
    dng::libraries_t output_libraries_;
    std::vector<std::string> bcf_samples_;

    std::vector<regions::contig_t> contigs_;
};

template<typename A>
inline
BcfPileup BcfPileup::open_and_setup(const A& arg) {
    if(arg.input.size() > 1) {
        throw std::runtime_error("processing more than one variant file at a time is not supported.");
    }

    BcfPileup mpileup;

    regions::set_regions(arg.region, &mpileup);

    if(mpileup.AddFile(arg.input[0].c_str()) == 0) {
        int errnum = mpileup.reader().handle()->errnum;
        throw std::runtime_error(bcf_sr_strerror(errnum));
    }

    return mpileup;
}

template<typename CallBack>
void BcfPileup::operator()(CallBack call_back) {
    assert(reader_.num_readers() == 1); // only support one reader

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
        // Unpack the record before calling call_back
        bcf_unpack(rec, BCF_UN_STR);

        // execute func
        call_back(rec);
    }
}

inline
int BcfPileup::AddFile(const char* filename) {
    assert(filename != nullptr);
    
    int index = reader_.num_readers();
    assert(index == 0); // only support one input file at this time

    if(reader_.AddReader(filename) == 0) {
        return 0;
    }
    ParseSampleLabels(index);
    ParseContigs(index);

    return 1;
}

// Parse PREFIX/SAMPLE/LIBRARY formatted labels
inline
void BcfPileup::ParseSampleLabels(int index) {

    auto reader = reader_.reader(index);

    bool needs_updating = false;

    int num_samples = bcf_hdr_nsamples(reader->header);
    auto bcf_samples = reader->header->samples;
    for(int i=0;i<num_samples;++i) {
        std::string bcf_sample = bcf_samples[i];
        std::string sample = trim_label_prefix_only_libraries(bcf_sample);
        if(sample.empty()) {
            continue;
        }
        std::string name;

        size_t pos = sample.find(DNG_LABEL_SEPARATOR_CHAR);        
        if(pos != std::string::npos) {
            name = sample.substr(pos+1);
            sample.erase(pos);
        }
        if(name.empty()) {
            name = sample;
        }
        pos = utility::find_position(input_libraries_.names, name);
        if(pos == input_libraries_.names.size()) {
            input_libraries_.names.push_back(std::move(name));
            input_libraries_.samples.push_back(std::move(sample));
            bcf_samples_.push_back(std::move(bcf_sample));
            needs_updating = true;
        } else {
            if(bcf_samples_[pos] != bcf_sample) {
                throw std::runtime_error("Multiple VCF/BCF column names for library '" + name + "': '" +
                    bcf_samples_[pos] + "' and '" + bcf_sample + "'.");
            }
            if(input_libraries_.samples[pos] != sample) {
                throw ("Multiple sample names for library '" + name + "': '" +
                    input_libraries_.samples[pos] + "' and '" + sample + "'.");
            }
        }
    }
    if(needs_updating) {
        ResetLibraries();
    }
}

template<typename R>
void BcfPileup::SelectLibraries(R &range) {
    assert(input_libraries_.names.size() == input_libraries_.samples.size());
    assert(input_libraries_.names.size() == bcf_samples_.size());

    output_libraries_ = {};

    // For every library in range, try to find it in input_libraries_
    std::string selector;
    for(auto it = boost::begin(range); it != boost::end(range); ++it) {
        auto pos = utility::find_position(input_libraries_.names, *it);
        if(pos == input_libraries_.names.size()) {
            // Do nothing if library was not found.
            continue;
        }
        output_libraries_.names.push_back(input_libraries_.names[pos]);
        output_libraries_.samples.push_back(input_libraries_.samples[pos]);

        if(!selector.empty()) {
            selector += ',';
        }
        selector += bcf_samples_[pos];
    }
    if(selector.empty()) {
        if(bcf_hdr_set_samples(reader_.reader(0)->header,nullptr,0) != 0) {
            throw std::runtime_error("Unable to exclude all VCF/BCF columns." );
        }
    } else {
        if(bcf_hdr_set_samples(reader_.reader(0)->header, selector.c_str(),0) != 0) {
            throw std::runtime_error("Unable to select VCF/BCF columns '" + selector + "'." );
        }
    }
}

void BcfPileup::ResetLibraries() {
    output_libraries_ = input_libraries_;
    if(bcf_hdr_set_samples(reader_.reader(0)->header, "-",0) != 0) {
        throw std::runtime_error("Unable to reset VCF/BCF columns." );
    }
}

inline
void BcfPileup::ParseContigs(int index) {
    auto reader = reader_.reader(index);

    const int num_contigs = reader->header->n[BCF_DT_CTG];

    for(int i=0;i<num_contigs;++i) {
        const char *name = reader->header->id[BCF_DT_CTG][i].key;
        auto pos = utility::find_position_if(contigs_,
            [name](const regions::contig_t& contig) { return name == contig.name; });
        int length = reader->header->id[BCF_DT_CTG][i].val->info[0];
        if(pos < contigs_.size()) {
            assert(contigs_[pos].length == length);
            continue;
        }
        contigs_.push_back({name, length});
    }
}

inline
int BcfPileup::SetRegions(const regions::contig_fragments_t &frags, bool one_indexed) {
    // convert ranges to the 1-based synced_bcf_reader format
    namespace karma = boost::spirit::karma;
    std::string str;
    str.reserve(128);
    for(size_t i = 0; i < frags.size(); ++i) {
        if(i != 0) {
            str += ',';
        }
        str += frags[i].contig_name;
        str += ':';
        karma::generate(std::back_inserter(str), frags[i].beg+(!one_indexed));
        str += '-';
        karma::generate(std::back_inserter(str), frags[i].end);
    }
    return reader_.SetRegions(str.c_str());
}

} //namespace io
} //namespace dng

#endif //DNG_IO_BCF_H
