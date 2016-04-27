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
#include <dng/utility.h>

namespace dng {
namespace io {

class Fasta {
public:
    typedef dng::utility::location_t location_t;

    explicit Fasta(const char* path);

    const std::string& path() const { return path_.native(); }
    bool is_open() const { return (bool)fai_; }

    explicit operator bool() const { return is_open(); }

    const char* FetchSequence(int contig, int first_pos=0, int last_pos=0x80000000);
    std::pair<int, int> FetchRange() const;

    char FetchBase(location_t loc);

    const char* name(int i) const;
    int length(int i) const;
    int count() const;
    bool contains(const char *seq) const;

    std::string FetchDictionary() const;

protected:
    bool UpdateBuffer(int contig, int first, int last);

    std::unique_ptr<faidx_t, void(*)(faidx_t *)> fai_{nullptr, fai_destroy};

    // Buffer to hold fetched sequences so that they are freed
    std::unique_ptr<char[], void(*)(void *)> buffer_{nullptr, free};
    int buffer_length_{0};
    int buffer_contig_{-1};
    int buffer_first_{-1};
    int buffer_last_{-1};

    boost::filesystem::path path_;
};

inline Fasta::Fasta(const char* path) : path_{path} {
    if(path_.empty())
        return;
    fai_.reset(fai_load(path_.c_str()));
    if(!fai_) {
        throw std::runtime_error("unable to open faidx-indexed reference file '"
                                     + path_.native() + "'.");
    }
}

inline
const char* Fasta::name(int i) const {
    assert(fai_);
    return faidx_iseq(fai_.get(), i);
}

inline
int Fasta::length(int i) const {
    assert(fai_);
    return faidx_seq_len(fai_.get(), name(i));
}

inline
int Fasta::count() const {
    assert(fai_);
    return faidx_nseq(fai_.get());
}

inline
bool Fasta::contains(const char *seq) const {
    assert(fai_);
    assert(seq != nullptr);
    return faidx_has_seq(fai_.get(),seq);
}

inline
bool Fasta::UpdateBuffer(int contig, int first, int last) {
    assert(fai_);
    assert(contig >= 0);
    if(buffer_contig_ == contig && buffer_first_ <= first && last <= buffer_last_) {
        return true;
    }
    // faidx_fetch_seq assumes inclusive
    last = last-1;
    if(last < first) {
        first = last;
    }
    // update first and last to reflect htslib rules
    int len = length(contig);
    if(first < 0) {
        first = 0;
    } else if(first >= len) {
        first = len-1;
    }
    if(last < 0) {
        last = 0;
    } else if(last >= len) {
        last = len-1;
    }

    std::cerr << contig << " " << name(contig) << std::endl;
    buffer_.reset(faidx_fetch_seq(fai_.get(), name(contig), first, last, &buffer_length_));
    if(buffer_length_ < 0) {
        buffer_contig_ = -1;
        return false;
    }
    buffer_contig_ = contig;
    buffer_first_ = first;
    buffer_last_ = last+1;
    return true;
}

inline
char Fasta::FetchBase(location_t loc) {
    if(!UpdateBuffer(utility::location_to_contig(loc), 0, 0x80000000)) {
        return 'N';
    }
    int pos = utility::location_to_position(loc);
    if(pos >= buffer_last_ || pos < buffer_first_) {
        return 'N';
    }
    return *(buffer_.get()+pos);
}

inline 
const char* Fasta::FetchSequence(int contig, int first, int last) {
    assert(fai_);
    assert(first <= last);
    if(!UpdateBuffer(contig,first,last) || first < buffer_first_ || (first-buffer_first_) > buffer_length_) {
        return nullptr;
    }
    return buffer_.get()+(first-buffer_first_);
}

inline
std::pair<int, int> Fasta::FetchRange() const {
    return {buffer_first_, buffer_last_};
}

inline
std::string Fasta::FetchDictionary() const {
    boost::filesystem::path dict = path_;
    if(dict.extension() == ".gz") {
        dict.replace_extension();
    }
    dict.replace_extension(".dict");
    boost::filesystem::ifstream file{dict};
    return slurp(file);
}

}} // namespace dng::io

#endif // DNG_IO_FASTA_H
