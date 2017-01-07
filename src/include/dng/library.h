/*
 * Copyright (c) 2017 Reed A. Cartwright
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
#ifndef DNG_LIBRARY_H
#define DNG_LIBRARY_H

#include <vector>
#include <string>

#include <dng/detail/graph.h>

namespace dng {

struct libraries_t {
    std::vector<std::string> names;
    std::vector<std::string> samples;
};

#define DNG_LABEL_PREFIX_GERMLINE "GL"
#define DNG_LABEL_PREFIX_SOMATIC  "SM"
#define DNG_LABEL_PREFIX_LIBRARY  "LB"
#define DNG_LABEL_SEPARATOR "/"
#define DNG_LABEL_SEPARATOR_CHAR '/'

}

#endif // DNG_LIBRARY_H
