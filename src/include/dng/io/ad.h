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

    explicit Ad(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        Open(filename,mode);
    }
    explicit Ad(const char *filename, std::ios_base::openmode mode = std::ios_base::in) : Ad(std::string{filename},mode) { }

    void Open(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        BinaryFile::Open(filename, mode);
        is_binary_ad_ = boost::iequals(type_, "ad");
        counter_ = 0;
    }

    int Write(const AlleleDepths& line);
    int WriteHeader();

    void contigs(std::vector<std::string> names) {
        contig_names_ = std::move(names);
    }
    const std::vector<std::string>& contigs() const {
        return contig_names_;
    }

private:
    int WriteAd(const AlleleDepths& line);
    int WriteTad(const AlleleDepths& line);

    struct contig_t {
        std::string name;
        int length;
        std::string attributes;
    };
    struct library_t {
        std::string name;
        std::string attributes;        
    };
    struct format_t {
        std::string name;
        uint16_t version;
        std::string attributes;        
    };

    bool is_binary_ad_{false};
    location_t last_location_{0};

    // Use rollover to trigger counter 
    uint16_t counter_{0};

    std::vector<std::string> contig_names_;
};

}}

#endif
