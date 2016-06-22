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

// Google Protocol buffers provides an implementation of variable-length int.
// We will reimplement it here.
//
// The encoding operates on unsigned integers of up to 64 bits in length.
// Each byte of the encoded value has the format:
// * bits 0-6: Seven bits of the number being encoded.
// * bit 7: Zero if this is the last byte in the encoding (in which
//   case all remaining bits of the number are zero) or 1 if
//   more bytes follow.
// The first byte contains the least-significant 7 bits of the number, the
// second byte (if present) contains the next-least-significant 7 bits,
// and so on.  So, the binary number 1011000101011 would be encoded in two
// bytes as "10101011 00101100".
//
// In theory, varint could be used to encode integers of any length.
// However, for practicality we set a limit at 64 bits.  The maximum encoded
// length of a number is thus 10 bytes.

#include <streambuf>
#include <utility>
#include <cstdint>
#include <limits.h>

#pragma once
#ifndef DNG_DETAIL_VARINT_H
#define DNG_DETAIL_VARINT_H

namespace dng { namespace detail { namespace varint {

// Check our assumptions
static_assert(CHAR_BIT == 8, "Sizeof char/byte is not 8 bits. Your machine is not supported." );

typedef std::basic_streambuf<char> bytebuf_t;

std::pair<uint64_t,bool> get_fallback(bytebuf *in, uint64_t first_byte);

inline
std::pair<uint64_t,bool> get(bytebuf *in) {
    assert(in != nullptr);
    // grab the first character
    bytebuf_t::int_type n = in->sbumpc();
    // if you have reached the end of stream, return error
    if(buf_type::traits_type::eq_int_type(n, buf_type::traits_type::eof())) {
        return {0,false};
    }
    // Convert back to a char and save in a 64-bit num.
    uint64_t u = buf_type::traits_type::to_char_type(n);
    // if MSB is not set, return the result
    if(!(u & 0x80)) {
        return {u,true};
    }
    // continue processing
    return varint::get_fallback(u);
}

}}} // dng::detail::varint
#endif


