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

#include <cstdlib>
#include <fstream>
#include <iterator>
#include <iosfwd>

#include <iostream>

#include <boost/range/iterator_range.hpp>

#include <dng/task/call.h>
#include <dng/pedigree.h>
#include <dng/pedigree_peeler.h>

using namespace dng::task;

// Helper function that mimics boost::istream_range
template<class Elem, class Traits> inline
    boost::iterator_range<std::istreambuf_iterator<Elem, Traits> >
istreambuf_range(std::basic_istream<Elem, Traits>& in)
{
    return boost::iterator_range<std::istreambuf_iterator<Elem, Traits>>(
        std::istreambuf_iterator<Elem, Traits>(in),
        std::istreambuf_iterator<Elem, Traits>());
}

int Call::operator()(const Call::argument_type &arg) {
	using namespace std;
	
	dng::Pedigree pedigree;

	// Parse pedigree from file	
	if(!arg.ped.empty()) {
		ifstream ped_file(arg.ped);
		if(!ped_file.is_open()) {
			throw std::runtime_error(
				"unable to open pedigree file '" + arg.ped + "'."
			);
		}
		pedigree.Parse(istreambuf_range(ped_file));
	} else {
		throw std::runtime_error("pedigree file was not specified.");
	}
	
	dng::PedigreePeeler peeler;
	peeler.Initialize(arg.theta, arg.mu);
	peeler.Construct(pedigree);
	
	dng::PedigreePeeler::IndividualBuffer buf;
	buf.resize(7, dng::PedigreePeeler::Vector10d::Zero());
	buf[1][1] = 1.0;
	buf[2][1] = 1.0;
	buf[3][1] = 1.0;
	buf[4][0] = 1.0;
	buf[5][0] = 1.0;
	buf[6][0] = 1.0;
	
	double d = peeler.CalculateLogLikelihood(buf);
	cout << exp(d) << endl;
	
	return EXIT_SUCCESS;
}
