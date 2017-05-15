/*
 * Copyright (c) 2009,2017 Reed A. Cartwright
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
#pragma once
#ifndef DNG_DETAIL_RANGEIO_H
#define DNG_DETAIL_RANGEIO_H

#include <cassert>

#include <vector>
#include <iostream>

namespace dng {

namespace rangeio {
// rangeio adapted from boost/tuple/tuple_io.hpp
// Copyright (C) 2001 Jaakko Jarvi (jaakko.jarvi@cs.utu.fi)
//               2001 Gary Powell (gary.powell@sierra.com)

enum struct manipulator_type {Open, Close, Delimiter};

namespace detail {

inline
int get_stream_index (manipulator_type m) {
    static const int stream_index[]
        = { std::ios::xalloc(), std::ios::xalloc(), std::ios::xalloc() };
    
    assert((int)m < sizeof(stream_index)/sizeof(int));

    return stream_index[(int)m];
}


template<typename CharType, typename CharTrait>
CharType get_manipulator(std::basic_ios<CharType, CharTrait>& i, 
                       manipulator_type m) {
    CharType c = static_cast<CharType>( i.iword(get_stream_index(m)) ); 
    // curly brackets and comma are the default manipulators
    if (!c) {
        switch(m) {
        case manipulator_type::Open:
            c = i.widen('{');
            break;
        case manipulator_type::Close:
            c = i.widen('}');
            break;
        case manipulator_type::Delimiter:
            c = i.widen(',');
            break;
        default:
            assert(false); // should not get here
            break;
        }
    }
    return c;
}

template<typename CharType, typename CharTrait>
void set_manipulator(std::basic_ios<CharType, CharTrait>& i, 
                   manipulator_type m, CharType c) {
    i.iword(get_stream_index(m)) = static_cast<long>(c);
}

template<typename CharType>
struct manipulator_t {
    const manipulator_type type;
    CharType character;

    explicit manipulator_t(manipulator_type m, const char c = 0)
        : type{m}, character{c} {
    }
  
    template<typename CharTrait>
    void set(std::basic_ios<CharType, CharTrait> &io) const {
        set_manipulator(io, type, character);
    }
};

template<typename CharType, typename CharTrait>
inline
std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const manipulator_t<CharType>& m) {
    m.set(o);
    return o;
}

template<typename R>
struct wrapped_t {
    const R& range;
};

template<typename CharType, typename CharTrait, typename R>
inline std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const wrapped_t<R> &v) {
    if(!o.good()) {
        return o;
    }

    auto it = std::begin(v.range);
    auto last = std::end(v.range); 
    
    const CharType l = detail::get_manipulator(o, manipulator_type::Open);
    if(l != 127) { 
        o << l;
    }

    if(it != last) {
        o << *(it++);
    }

    const CharType d = detail::get_manipulator(o, manipulator_type::Delimiter);
    if(d != 127) {
        for(;it != last;++it) {
            o << d << *it;
        }
    } else {
        for(;it != last;++it) {
            o << *it;
        }
    }

    const CharType r = detail::get_manipulator(o, manipulator_type::Close);
    if(r != 127) {
        o << r;
    }
    
    return o;
}


} //namespace detail

template<typename CharType>
inline detail::manipulator_t<CharType> set_open(const CharType c) {
   return {manipulator_type::Open, c};
}

template<typename CharType>
inline detail::manipulator_t<CharType> set_close(const CharType c) {
   return {manipulator_type::Close, c};
}

template<typename CharType>
inline detail::manipulator_t<CharType> set_delimiter(const CharType c) {
   return {manipulator_type::Delimiter, c};
}

template<typename R>
inline detail::wrapped_t<R> wrap(const R& r) {
    return {r};
}

} //namespace dng::rangeio
} //namespace dng
#endif

