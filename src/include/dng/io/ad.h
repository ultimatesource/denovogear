/*
 * Copyright (c) 2015-2016 Reed A. Cartwright
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

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <dng/io/file.h>
#include <dng/io/utility.h>
#include <dng/utility.h>
#include <dng/depths.h>

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
        contig_t() {}
        contig_t(std::string n, int l, std::string a = {}) :
            name{std::move(n)}, length{l}, attributes{std::move(a)} {}

        std::string name;
        int length;
        std::vector<std::string> attributes;
    };

    struct library_t {
        library_t() {}
        library_t(std::string n, std::string s, std::string a = {}) :
            name{std::move(n)}, sample{std::move(s)}, attributes{std::move(a)} {}

         std::string name;
         std::string sample;
         std::vector<std::string> attributes;
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

    const std::vector<library_t>& libraries() const {
        return libraries_;
    }

    const library_t& library(std::vector<library_t>::size_type pos) const {
        return libraries_[pos];
    }

    std::vector<contig_t>::size_type
    AddContig(std::string name, int length, std::string attributes = {} ) {
        std::vector<contig_t>::size_type pos = contigs_.size();
        contig_map_.emplace(name, pos);
        contigs_.emplace_back(name, length, attributes);
        return pos;
    }

    std::vector<contig_t>::size_type
    AddLibrary(std::string name, std::string sample, std::string attributes = {} ) {
        std::vector<contig_t>::size_type pos = contigs_.size();
        libraries_.emplace_back(name, sample, attributes);
        num_libraries_ = libraries_.size();
        last_data_.assign(num_libraries_, 0);
        return pos;
    }

    Format format() const { return format_; }

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
    std::vector<library_t> libraries_;
    std::size_t num_libraries_{0};

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
        result = ReadHeaderAd();
    } else {
        result = ReadHeaderTad();
    }
    if(result == 0) {
        return 0;
    }
    // Insert contigs into the search tree
    for(int i=0; i < contigs_.size(); ++i) {
        contig_map_.emplace(contigs_[i].name,i);
    }
    // Save the number of libraries
    num_libraries_ = libraries_.size();
    last_location_ = 0;
    last_data_.assign(num_libraries_,0);
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
    libraries_ = ad.libraries_;
    num_libraries_ = ad.num_libraries_;
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



}}

#endif
