/*
 * Copyright (c) 2016-2017 Reed A. Cartwright
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

#define BOOST_TEST_MODULE dng::multithread

#include <dng/multithread.h>

#include "../testing.h"

#include <atomic>
#include <iostream>

BOOST_AUTO_TEST_CASE(test_basicpool) {
	using namespace std;
	using namespace dng::multithread;

	atomic<int> result{0};

	auto f = [&result](int x) {
		result += x;
	};

	{
		BasicPool<int> pool(f,4);
		for(int i=0;i<=100;++i) {
			pool.Enqueue(i);
		}
	}
	BOOST_CHECK(result == (100*101)/2);
}