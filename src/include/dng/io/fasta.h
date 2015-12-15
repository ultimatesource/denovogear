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
#ifndef DNG_IO_FASTA_H
#define DNG_IO_FASTA_H

#include <cstdlib>
#include <memory>
#include <utility>
#include <string>

#include <htslib/faidx.h>

#include <boost/filesystem.hpp>

#include <dng/io/utility.h>

namespace dng {
namespace io {

class Fasta {
public:
    explicit Fasta(const char* path);

    const std::string& name() const { return filename_; }
    bool is_open() const { return fai_; }

    operator bool() const { return is_open(); }

    std::pair<const char*,int> FetchSequence(const char* target_name, int first_pos=0, int last_pos=0x7FFFFFFF);
    std::string FetchSequenceAsString(const char* target_name, int first_pos=0, int last_pos=0x7FFFFFFF);

    std::string FetchDictionary() const;

protected:
    std::unique_ptr<faidx_t, void(*)(faidx_t *)> fai_{nullptr, fai_destroy};
    std::unique_ptr<char[], void(*)(void *)> ref_{nullptr, free};
    int ref_length_{0};

    boost::filesystem::path path_;
};

inline Fasta::Fasta(const char* path) : path_{path} {
    if(path_.empty())
        return;
    fai_.reset(fai_load(path_.c_str()));
    if(!fai) {
        throw std::runtime_error("unable to open faidx-indexed reference file '"
                                     + filename_ + "'.");
    }
}

inline 
std::pair<const char*,int> Fasta::FetchSeq(const char* target_name, int first, int last) {
    assert(fai_);
    ref_.reset(faidx_fetch_seq(fai.get(), target_name, first, last, &ref_length_));
    return {ref_.get(), ref_length_};
}

inline
std::string Fasta::FetchSeqAsString(const char* target_name, int first_pos, int last_pos) {
    auto result = FetchSeq(target_name, first_pos, last_pos);
    return {result.first, result.second};
}

inline
std::string Fasta::FetchDictionary() const {
    boost::filesystem::path dict = path_;
    if(dict.extension() == ".gz") {
        dict.replace_extension();
    }
    dict.replace_extension(".dict");
    return slurp_file(boost::filesystem::ifstream{dict});
}

}
} // namespace dng::io

#endif // DNG_IO_FASTA_H
