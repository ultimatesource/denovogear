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
#ifndef DNG_TASK_H
#define DNG_TASK_H

#include <vector>

namespace dng {

namespace task {
struct arg_t {
    std::vector< std::string > input;
};
}

template<typename A>
class Task {
public:
    typedef A   argument_type;
    typedef int result_type;

    int operator()(argument_type &) {
        return EXIT_SUCCESS;
    }
};

} // namespace dng

#endif // DNG_TASK_H

