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
#ifndef CXX_HTS_HTS_H
#define CXX_HTS_HTS_H

#include <htslib/hts.h>
#include <cstdlib>
#include <memory>
#include <cassert>
#include <cstring>
#include <string>
#include <sstream>

namespace hts {

class File {
public:
    File(const char *file, const char *mode) :
        fp_{hts_open(file, mode), hts_close} { }

    int SetFaiFileName(const char *fn) {
        return hts_set_fai_filename(handle(), fn);
    }
    int SetThreads(int n) {
        return hts_set_threads(handle(), n);
    }

    bool is_open() const { return (bool)fp_; }
    bool is_bin() const { return handle()->is_bin; }
    bool is_write() const { return handle()->is_write; }
    bool is_cram() const { return handle()->is_cram; }
    bool is_compressed() const {
        return handle()->format.compression;
    }

    const char *name() const {
        return handle()->fn;
    }
    htsFormat format() const {
        return handle()->format;
    }
    std::string format_description() const {
        std::unique_ptr<char[], void(*)(void *)> s{hts_format_description(&handle()->format), std::free};
        return {s.get()};
    }

// Changed to public so that dng-dnm and dng-call can iterate through each record
//protected:
    htsFile *handle() {
        assert(fp_);
        return fp_.get();
    }
    const htsFile *handle() const {
        assert(fp_);
        return fp_.get();
    }

private:
    std::unique_ptr<htsFile, int(*)(htsFile *)> fp_; // the file handle
};

// convert hts version string to a numeric value
inline
unsigned long version_parse(const char *str) {
    unsigned long v[3] = {0,0,0};
    char *end;
    for(int i=0;i<3;++i) {
        v[i] = std::strtoul(str, &end, 10);
        if(!*end) {
            break;
        }
        if(*end != '.' && *end != '-') {
            return 0; // unexpected format
        }
        str = end+1;
    }
    return v[0]*100*100+v[1]*100+v[2];
}

inline unsigned long version() {
    return version_parse(hts_version());
}

namespace detail {
template<typename S>
std::string make_data_url(const S& str) {
    // data: url format changed in htslib 1.4
    std::string bamstr = "data:";
    if(hts::version() >= 10400 ) {
        bamstr += ",";
    }
    return bamstr+str;
}

} // namespace detail

};

#endif
