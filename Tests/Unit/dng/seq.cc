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

#define BOOST_TEST_MODULE dng::seq

#include <dng/seq.h>

#include "../testing.h"

#include <vector>
#include <string>

using namespace std;

BOOST_AUTO_TEST_CASE(test_char_index)
{
  BOOST_CHECK_EQUAL(dng::seq::char_index('N'), 4);
  BOOST_CHECK_EQUAL(dng::seq::char_index('A'), 0);
  BOOST_CHECK_EQUAL(dng::seq::char_index('C'), 1);
  BOOST_CHECK_EQUAL(dng::seq::char_index('G'), 2);
  BOOST_CHECK_EQUAL(dng::seq::char_index('T'), 3);
}

BOOST_AUTO_TEST_CASE(test_indexed_char)
{
  BOOST_CHECK_EQUAL(dng::seq::indexed_char(4), 'N');
  BOOST_CHECK_EQUAL(dng::seq::indexed_char(0), 'A');
  BOOST_CHECK_EQUAL(dng::seq::indexed_char(1), 'C');
  BOOST_CHECK_EQUAL(dng::seq::indexed_char(2), 'G');
  BOOST_CHECK_EQUAL(dng::seq::indexed_char(3), 'T');
}
