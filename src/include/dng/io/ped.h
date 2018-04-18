/*
 * Copyright (c) 2014-2015 Reed A. Cartwright
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
#ifndef DNG_IO_PED_H
#define DNG_IO_PED_H

#include <string>

#include <dng/utility.h>
#include <dng/pedigree.h>
#include <dng/io/utility.h>
#include <dng/io/file.h>

namespace dng {
namespace io {

class Ped : public File {
public:
    using File::File;

    Pedigree Parse();
};

inline
Pedigree Ped::Parse() {
    if(!is_open()) {
        return {};
    }
    return Pedigree::parse_text(io::istreambuf_range(stream_.rdbuf()));
}

inline
Pedigree parse_ped(const std::string &path) {
    if(path.empty()) {
        throw std::invalid_argument("Path to ped file is empty.");
    }
    Ped file(path);
    if(!file.is_open()) {
        throw std::runtime_error("Unable to open ped file '" + path + "'.");
    }
    return file.Parse();
}

}
} // namespace dng::io

#endif // DNG_PEDIGREE_H
