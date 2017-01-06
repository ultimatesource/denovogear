/*
 * Copyright (c) 2015-2017 Reed A. Cartwright
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
#ifndef DNG_IO_AD_H
#define DNG_IO_AD_H

#include <string>
#include <utility>
#include <unordered_map>
#include <numeric>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/functions.hpp>
#include <boost/range/algorithm/find.hpp>

#include <dng/io/file.h>
#include <dng/io/utility.h>
#include <dng/utility.h>
#include <dng/depths.h>
#include <dng/library.h>

#include <dng/detail/unit_test.h>


namespace dng {
namespace io {

/*
@ID FF:TAD VN:0.1 SO:coordinate
@SQ SN:scaffold_1 LN:100
@RG ID:A AD:A
@AD ID:A
@CO This is a comment
@CO POS TYPE A B C D
*/

class Ad : public BinaryFile {
public:
    typedef dng::pileup::AlleleDepths AlleleDepths;

    enum class Format {
        AD = 0, TAD = 1
    };

    struct contig_t {
        std::string name;
        int length;
    };

    explicit Ad(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        Open(filename,mode);
    }
    explicit Ad(const char *filename, std::ios_base::openmode mode = std::ios_base::in) : Ad(std::string{filename},mode) { }

    void Open(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in);

    int ReadHeader();
    int WriteHeader();
    void CopyHeader(const Ad& other);
    void AddHeaderLines(const std::string& lines);

    int Read(AlleleDepths *pline);
    int Write(const AlleleDepths& line);

    const std::vector<contig_t>& contigs() const {
        return contigs_;
    }

    const contig_t& contig(std::vector<contig_t>::size_type pos) const {
        return contigs_[pos];
    }

    const libraries_t& libraries() const {
        return output_libraries_;
    }
    size_t num_libraries() const {
        return output_libraries_.names.size();
    }

    size_t
    AddHeaderContig(std::string name, int length, std::string attributes = {} );

    size_t
    AddHeaderLibrary(std::string name, std::string sample, std::string attributes = {} );

    Format format() const { return format_; }

    template<typename R>
    void SelectLibraries(R &range);

    void ResetLibraries();
    void ResetContigs();

private:
    int ReadHeaderAd();
    int WriteHeaderAd();
    int ReadAd(AlleleDepths *line);
    int WriteAd(const AlleleDepths& line);

    int ReadHeaderTad();
    int WriteHeaderTad();
    int ReadTad(AlleleDepths *line);
    int WriteTad(const AlleleDepths& line);

    std::string HeaderString() const;
    template<typename It>
    void ParseHeaderTokens(It it, It it_last);

    void Clear();

    void FinishHeader();

    struct format_t {
         std::string name;
         uint16_t version;
         std::vector<std::string> attributes;        
    };

    Format format_{Format::AD};
    location_t last_location_{0};
    AlleleDepths::data_t last_data_;

    // Use rollover to trigger counter 
    uint16_t counter_{0};

    // Header Information
    format_t id_;
    std::vector<contig_t> contigs_;
    std::vector<std::string> contig_attributes_;

    struct ad_libraries_t : dng::libraries_t {
        std::vector<std::string> attributes;
    };

    ad_libraries_t input_libraries_;
    ad_libraries_t output_libraries_;

    std::vector<size_t> indexes_;

    std::vector<std::string> extra_headers_;

    std::unordered_map<std::string, int> contig_map_;

    DNG_UNIT_TEST(unittest_dng_io_ad);
};

inline
void Ad::Open(const std::string &filename, std::ios_base::openmode mode) {
    BinaryFile::Open(filename, mode);
    if(boost::iequals(type_label_, "ad")) {
        format_ = Format::AD;
    } else {
        format_ = Format::TAD;
    }
    // Clear header information
    Clear();
}

inline
int Ad::ReadHeader() {
    // Seek to the beginning of the stream
    stream_.seekg(0);
    if(stream_.bad()) {
        return 0;
    }
    stream_.clear();    
    // Clear header information
    Clear();

    int result = 0;
    if(format_ == Format::AD) {
        return ReadHeaderAd();
    } else {
        return ReadHeaderTad();
    }

    return result;
}

inline
int Ad::WriteHeader() {
    // Seek to the beginning of the stream
    stream_.seekp(0);
    if(stream_.bad()) {
        return 0;
    }
    stream_.clear();
    if(format_ == Format::AD) {
        return WriteHeaderAd();
    } else {
        return WriteHeaderTad();
    }
}

inline
void Ad::CopyHeader(const Ad& ad) {
    Clear();

    contigs_ = ad.contigs_;
    contig_attributes_ = ad.contig_attributes_;

    input_libraries_ = ad.input_libraries_;

    output_libraries_ = ad.output_libraries_;

    contig_map_ = ad.contig_map_;
    extra_headers_ = ad.extra_headers_;
}

inline
int Ad::Read(AlleleDepths *pline) {
    // If stream_ has issues, i.e. eof, return 0.
    if(!stream_) {
        return 0;
    }
    if(format_ == Format::AD) {
        return ReadAd(pline);
    } else {
        return ReadTad(pline);
    }
}

inline
int Ad::Write(const AlleleDepths& line) {
    // If stream_ has issues, return 0.
    if(!stream_) {
        return 0;
    }
    if(format_ == Format::AD) {
        return WriteAd(line);
    } else {
        return WriteTad(line);
    }
}

inline size_t
Ad::AddHeaderContig(std::string name, int length, std::string attributes) {
    auto pos = contigs_.size();
    contigs_.push_back({name, length});
    contig_attributes_.push_back(attributes);
    ResetContigs();
    return pos;
}

inline size_t
Ad::AddHeaderLibrary(std::string name, std::string sample, std::string attributes) {
    auto pos = input_libraries_.names.size();
    input_libraries_.names.push_back(std::move(name));
    input_libraries_.samples.push_back(std::move(sample));
    input_libraries_.attributes.push_back(std::move(attributes));
    ResetLibraries();
    return pos;
}

inline
void Ad::ResetLibraries() {
    output_libraries_ = input_libraries_;

    indexes_.resize(input_libraries_.names.size());
    std::iota(indexes_.begin(),indexes_.end(),0);

    last_data_.assign(output_libraries_.names.size(), 0);
    last_location_ = 0;
}

template<typename R>
void Ad::SelectLibraries(R &range) {
    // Clear all output libraries
    indexes_.assign(input_libraries_.names.size(),-1);
    output_libraries_ = {};

    // For every library in range, try to find it in input_libraries_
    size_t k=0;
    for(auto it = boost::begin(range); it != boost::end(range); ++it) {
        auto pos = boost::distance(boost::find<boost::return_begin_found>(input_libraries_.names, *it));
        if(pos == input_libraries_.names.size()) {
            continue;
        }
        output_libraries_.names.push_back(input_libraries_.names[pos]);
        output_libraries_.samples.push_back(input_libraries_.samples[pos]);
        output_libraries_.attributes.push_back(input_libraries_.attributes[pos]);

        indexes_[pos] = k++;
    }

    last_data_.assign(output_libraries_.names.size(), 0);
    last_location_ = 0;
}


inline
void Ad::ResetContigs() {
    contig_map_.clear();

    // Insert contigs into the search tree
    for(int i=0; i < contigs_.size(); ++i) {
        contig_map_.emplace(contigs_[i].name,i);
    }
}

inline
void Ad::FinishHeader() {
    ResetLibraries();
    ResetContigs();
}

}}

#endif
