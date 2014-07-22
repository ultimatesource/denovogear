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
#ifndef DNG_APP_CALL_H
#define DNG_APP_CALL_H

#include <dng/app.h>

namespace dng { namespace app {

namespace call {

struct arg_t : public base::arg_t {
	// use X-Macros to specify argument variables
#	define XM(lname, sname, desc, type, def) type XV(lname) ;
#	include "call.xmh"
#	undef XM

	// This function connects variables in this structure to command line
	// arguments.
	void add_to(po::options_description &desc) {
		desc.add_options()
#		define XM(lname, sname, desc, type, def) ( \
			XS(lname) IFD(sname, "," BOOST_PP_STRINGIZE sname), \
			po::value< type >(&arg.XV(lname))->default_value(def), \
			desc )	
#		include "call.xmh"
#		undef XM
		;
		base::arg_t::add_to(desc);		
	}
};

};

class Call : public Base<call::arg_t> {
public:
	typedef Base<call::arg_t> base_t;
	Call(int argc, char* argv[]) : base_t(argc,argv) {
	}

	virtual int Run() {
		return EXIT_SUCCESS;
	}
};

}}

#endif
