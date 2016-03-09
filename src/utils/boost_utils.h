/*
 * Copyright (c) 2016 Steven H. Wu
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
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
#ifndef DENOVOGEAR_BOOST_UTILS_H
#define DENOVOGEAR_BOOST_UTILS_H

#include <boost/range/iterator_range_core.hpp>
namespace utils {


// Helper function that mimics boost::istream_range
template<class Elem, class Traits> inline
boost::iterator_range <std::istreambuf_iterator<Elem, Traits>> istreambuf_range(
        std::basic_istream <Elem, Traits> &in) {
    return boost::iterator_range < std::istreambuf_iterator < Elem, Traits >> (
            std::istreambuf_iterator<Elem, Traits>(in),
                    std::istreambuf_iterator<Elem, Traits>());
}


} // namespace utils

#endif //DENOVOGEAR_BOOST_UTILS_H
