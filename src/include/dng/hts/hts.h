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

};

#endif
