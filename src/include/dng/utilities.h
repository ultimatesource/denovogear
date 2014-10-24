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
#ifndef DNG_UTILITIES_H
#define DNG_UTILITIES_H

#include <boost/algorithm/string/predicate.hpp>
#include <tuple>

namespace dng { namespace util {

template<class A, class B, std::size_t N>
std::size_t key_switch(A &ss, const B (&key)[N]) {
	using boost::algorithm::istarts_with;
	for(std::size_t i=0;i<N;++i) {
		if(istarts_with(key[i], ss))
			return i;
	}
	return (std::size_t)-1;
}

template<class A, class B, std::size_t N>
const B& key_switch_tuple(A &ss, const B (&key)[N], const B& default_value) {
	using boost::algorithm::istarts_with;
	for(std::size_t i=0;i<N;++i) {
		if(istarts_with(std::get<0>(key[i]), ss))
			return key[i];
	}
	return default_value;
}

}};

#endif
