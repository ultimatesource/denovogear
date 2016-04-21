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

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <dng/io/file.h>
#include <dng/io/utility.h>
#include <dng/utility.h>
#include <dng/depths.h>


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
        std::string attributes;
    };

    explicit Ad(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        Open(filename,mode);
    }
    explicit Ad(const char *filename, std::ios_base::openmode mode = std::ios_base::in) : Ad(std::string{filename},mode) { }

    void Open(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        BinaryFile::Open(filename, mode);
        if(boost::iequals(type_label_, "ad")) {
            format_ = Format::AD;
        } else {
            format_ = Format::TAD;
        }
        counter_ = 0;
    }

    int Write(const AlleleDepths& line);
    int WriteHeader();

    int ReadHeader();

    const std::vector<contig_t>& contigs() const {
        return contigs_;
    }

    const contig_t& contig(std::vector<contig_t>::size_type pos) const {
        return contigs_[pos];
    }

    std::vector<contig_t>::size_type
    AddContig(std::string name, int length, std::string attributes = {} ) {
        contigs_.emplace_back(name, length, attributes);
        return contigs_.size()-1;
    }

private:
    int WriteAd(const AlleleDepths& line);
    int WriteTad(const AlleleDepths& line);

    int ReadHeaderTad();
    int ReadHeaderAd();

    // struct format_t {
    //     std::string name;
    //     uint16_t version;
    //     std::string attributes;        
    // };
    // struct library_t {
    //     std::string id;
    //     std::string sample;
    //     std::string attributes;        
    // };

    Format format_{Format::AD};
    location_t last_location_{0};

    // Use rollover to trigger counter 
    uint16_t counter_{0};

    std::vector<contig_t> contigs_;
};

}}

#endif
