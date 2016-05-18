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
#ifndef DNG_SEQ_H
#define DNG_SEQ_H

#include <htslib/hts.h>

#include <dng/detail/unit_test.h>

namespace dng {
namespace seq {

// convert a nucleotide character into a 4-bit representation
inline uint8_t encode_base(unsigned char x) {
	assert(0 <= seq_nt16_table[x] && seq_nt16_table[x] < 16);
    return seq_nt16_table[x]; // from htslib
}

// convert a 4-bit nucleotide into a character
inline char decode_base(uint8_t x) {
    assert(0 <= x && x < 16);	
    return seq_nt16_str[x]; // from htslib
}

// convert a 4-bit nucleotide into a index in [0-4]
inline int base_index(uint8_t x) {
    static constexpr int table[] = {4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4};
    assert(0 <= x && x < 16);
    std::size_t y = table[x];
    assert(0 <= y && y < 5);
    return y;
}

// convert an index into a 4-bit nucleotide
inline uint8_t indexed_base(std::size_t x) {
    static constexpr uint8_t table[] = {1, 2, 4, 8, 15};
    assert(0 <= x && x < 5);
    uint8_t y = table[x];
    assert(0 <= y && y < 16);
    return y;
}

// convert an index into an ASCII character
inline char indexed_char(std::size_t x) {
    static constexpr char table[] = "ACGTNNNN";
    assert(0 <= x && x < 5);
    return table[x];
}

// convert a 4-bit nucleotide into a index in [0-4]
inline int char_index(uint8_t x) {
    return base_index(encode_base(x));
}

}
} // namespace dng::seq

#endif

