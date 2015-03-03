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
#ifndef CXX_HTS_EXTRA_H
#define CXX_HTS_EXTRA_H

#include <string>
#include <algorithm>

namespace hts { namespace extra {

// TODO: Should this stay here?  Or move out of the hts wrapper?
// TODO: Undestand .hidden files ?
std::pair<std::string,std::string> extract_file_type(const std::string& path) {
	if(path.empty())
		return {};

	std::string::size_type x=path.length()-1;

	for(std::string::size_type u=0; u < path.length(); ++u) {
		if(path[u] == ':' && u > 1)
			return {path.substr(0,u), path.substr(u+1)};
		if(path[u] == '.')
			x = u;
	}
	return {path.substr(x+1), path};
}

}}

#endif