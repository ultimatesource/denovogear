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
#ifndef DNG_IO_UTILITY_H
#define DNG_IO_UTILITY_H

#include <string>
#include <iterator>
#include <iosfwd>
#include <istream>
#include <fstream>

#include <boost/range/iterator_range.hpp>

namespace dng {
namespace io {

template<typename X, typename T = std::char_traits<X>, typename A = std::allocator<X>>
inline
std::basic_string<X,T,A> slurp(std::basic_ifstream<X,T>& input) {
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

inline std::string slurp(const char *filename, std::ios_base::openmode mode = std::ios_base::in) {
    std::ifstream in{filename,mode};
    return slurp(in);
}

inline std::string slurp(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in) {
    std::ifstream in{filename,mode};
    return slurp(in);
}

inline std::wstring wslurp(const char *filename, std::ios_base::openmode mode = std::ios_base::in) {
    std::wifstream in{filename,mode};
    return slurp(in);
}

inline std::wstring wslurp(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in) {
    std::wifstream in{filename,mode};
    return slurp(in);
}

// Helper function that mimics boost::istream_range
template<class Elem, class Traits> inline
boost::iterator_range<std::istreambuf_iterator<Elem, Traits> >
istreambuf_range(std::basic_istream<Elem, Traits> &in) {
    return boost::iterator_range<std::istreambuf_iterator<Elem, Traits>>(
               std::istreambuf_iterator<Elem, Traits>(in),
               std::istreambuf_iterator<Elem, Traits>());
}

inline bool at_slurp(std::string &ss, std::ios_base::openmode mode = std::ios_base::in) {
    if(ss.empty() || ss[0] != '@')
        return false;
    std::ifstream in{ss.c_str()+1, mode};
    ss = slurp(in);
    return true;
}

}
} //namespace dng::io

#endif
