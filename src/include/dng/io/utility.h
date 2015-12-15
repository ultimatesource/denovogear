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
#ifndef DNG_IO_UTILITY_H
#define DNG_IO_UTILITY_H

#include <string>
#include <istream>

namespace dng {
namespace io {

template<typename X, typename T = std::char_traits<X>, typename A = std::allocator<X> >
inline
std::basic_string<X,T,A> slurp_file(std::istream<X,T>& input) {
    std::basic_string<X,T,A> ret;
    if(!input) {
        return {};
    }

    input.seekg(0, std::ios::end);
    ret.resize(input.tellg());
    input.seekg(0, std::ios::beg);
    input.read(&ret[0], ret.size());
    input.close();
    return ret;
}

inline std::string slurp_file(const char *filename, std::ios_base::openmode mode = std::ios_base::in) {
    return slurp_file(std::ifstream{filename,mode});
}

inline std::string slurp_file(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in) {
    return slurp_file(std::ifstream{filename,mode});
}

inline std::wstring wslurp_file(const char *filename, std::ios_base::openmode mode = std::ios_base::in) {
    return slurp_file(std::wifstream{filename,mode});
}

inline std::wstring wslurp_file(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in) {
    return slurp_file(std::wifstream{filename,mode});
}


}
} //namespace dng::io

#endif
