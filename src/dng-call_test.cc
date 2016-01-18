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

#include <iostream>

#include <dng/app.h>
#include <dng/task/call.h>

#ifdef DNG_DEVEL
#   include <boost/timer/timer.hpp>
#endif

// http://www.boost.org/development/requirements.html
// http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml

typedef dng::CommandLineApp<dng::task::Call> CallApp;

int main(int argc, char *argv[]) {
#ifdef DNG_DEVEL
    boost::timer::auto_cpu_timer measure_speed(std::cerr);
#endif
    try {
        return CallApp(argc, argv)();
    } catch(std::exception &e) {
        std::cerr << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}
