/*
 * Copyright (c) 2016 Reed A. Cartwright
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

#include <streambuf>

#include "portable_endian.h"

#pragma once
#ifndef DNG_DETAIL_NTF8_H
#define DNG_DETAIL_NTF8_H

namespace dng { namespace detail { namespace ntf8 {

template<typename CharT, typename Traits>
int get32(std::basic_streambuf<CharT,Traits> *in, int32_t *r) {
    typedef typename std::basic_streambuf<CharT,Traits> buf_type;
    assert(r != nullptr);
    assert(in != nullptr);
    typename buf_type::int_type n = in->sgetc();
    if(n == buf_type::traits_type::eof()) {
        return 0;
    }
    uint8_t x = n;
    if(x < 0x80) {
        // 0bbb bbbb
        in->sbumpc();
        *r = x;
        return 1;
    } else if(x < 0xC0) {
        // 10bb bbbb bbbb bbbb
        uint16_t u = 0;
        if(in->sgetn((char*)&u,2) != 2) {
            return 0;
        }
        *r = be16toh(u) & 0x3FFF;
        return 2;
    } else if(x < 0xE0) {
        // 110b bbbb bbbb bbbb bbbb bbbb
        uint32_t u = 0;
        if(in->sgetn((char*)&u,3) != 3) {
            return 0;
        }
        *r = (be32toh(u) >> 8) & 0x1FFFFF;
        return 3;
    } else if(x < 0xF0) {
        // 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint32_t u = 0;
        if(in->sgetn((char*)&u,4) != 4) {
            return 0;
        }
        *r = be32toh(u) & 0x0FFFFFFF;
        return 4;
    } else {
        assert(x == 0xF0); // If this fails then we may be reading a 64-bit number
        // 1111 0000 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        in->sbumpc();
        uint32_t u = 0;        
        if(in->sgetn((char*)&u,4) != 4) {
            return 0;
        }
        *r = be32toh(u);
        return 5;
    }
}

template<typename CharT, typename Traits>
int get64(std::basic_streambuf<CharT,Traits> *in, int64_t *r) {
    typedef typename std::basic_streambuf<CharT,Traits> buf_type;    
    assert(r != nullptr);
    assert(in != nullptr);
    typename buf_type::int_type n = in->sgetc();
    if(n == buf_type::traits_type::eof()) {
        return 0;
    }
    uint8_t x = n;
    if(x < 0x80) {
        // 0bbb bbbb
        in->sbumpc();
        *r = x;
        return 1;
    } else if(x < 0xC0) {
        // 10bb bbbb bbbb bbbb
        uint16_t u = 0;
        if(in->sgetn((char*)&u,2) != 2) {
            return 0;
        }
        *r = be16toh(u) & 0x3FFF;
        return 2;
    } else if(x < 0xE0) {
        // 110b bbbb bbbb bbbb bbbb bbbb
        uint32_t u = 0;
        if(in->sgetn((char*)&u,3) != 3) {
            return 0;
        }
        *r = (be32toh(u) >> 8) & 0x1FFFFF;
        return 3;
    } else if(x < 0xF0) {
        // 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint32_t u = 0;
        if(in->sgetn((char*)&u,4) != 4) {
            return 0;
        }
        *r = be32toh(u) & 0x0FFFFFFF;
        return 4;
    } else if(x < 0xF8) {
        // 1111 0bbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t u = 0;
        if(in->sgetn((char*)&u,5) != 5) {
            return 0;
        }
        *r = (be64toh(u) >> 24) & 0x07FFFFFFFF;
        return 5;
    } else if(x < 0xFC) {
        // 1111 10bb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t u = 0;
        if(in->sgetn((char*)&u,6) != 6) {
            return 0;
        }
        *r = (be64toh(u) >> 16) & 0x03FFFFFFFFFF;
        return 6;
    } else if(x < 0xFE) {
        // 1111 110b bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t u = 0;
        if(in->sgetn((char*)&u,7) != 7) {
            return 0;
        }
        *r = (be64toh(u) >> 8) & 0x01FFFFFFFFFFFF;
        return 7;
    } else if(x < 0xFF) {
        // 1111 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t u = 0;
        if(in->sgetn((char*)&u,8) != 8) {
            return 0;
        }
        *r = be64toh(u) & 0x00FFFFFFFFFFFFFF;
        return 8;
    } else {
        // 1111 1111 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        in->sbumpc();
        uint64_t u = 0;        
        if(in->sgetn((char*)&u,8) != 8) {
            return 0;
        }
        *r = be64toh(u);
        return 9;
    }
}

}}} // dng::detail::ntf8
#endif
