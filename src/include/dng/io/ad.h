/*
 * Copyright (c) 2015 Reed A. Cartwright
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


#include <dng/io/utility.h>
#include <dng/utility.h>

namespace dng {
namespace io {

class Ad {
public:
    explicit Ad(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        Open(filename,mode);
    }
    explicit Ad(const char *filename, std::ios_base::openmode mode = std::ios_base::in) : Ad(std::string{filename},mode) { }

    void Open(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        auto split = extract_file_type(filename);
        // assume anything not marked as ad is in text/tad format
        is_binary_ad_ = boost::iequals(split.first, "ad");
        if(is_binary_ad_) {
            mode |= std::ios::bin;
        }
        path_ = split.second;
        file_.open(path_, mode);
    }

    bool is_open() const { return file_.is_open() };
    operator bool() const { return is_open(); }



protected:
    boost::filesystem::path path_;
    boost::filesystem::fstream file_;
    bool is_binary_ad_{false};
};

}}

#endif
