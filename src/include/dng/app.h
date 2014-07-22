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
#ifndef DNG_APP_H
#define DNG_APP_H

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/logic/tribool.hpp>
#include <boost/logic/tribool_io.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// Support the use of tribools in boost::program_options
namespace boost {
void validate(boost::any& v, const std::vector<std::string>& xs, boost::tribool*, int) {
    using namespace boost::program_options;
    validators::check_first_occurrence(v);
	std::string s(validators::get_single_string(xs, true));

    for (size_t i = 0; i < s.size(); ++i)
        s[i] = char(tolower(s[i]));

    if (s.empty() || s == "on" || s == "yes" || s == "1" || s == "true")
		v = boost::any(boost::tribool(true));
    else if (s == "off" || s == "no" || s == "0" || s == "false")
		v = boost::any(boost::tribool(false));
    else if (s == "null" || s == "maybe" || s == "2" || s == "indeterminate")
		v = boost::any(boost::tribool(boost::indeterminate));
    else
        boost::throw_exception(validation_error(validation_error::invalid_option_value, s));
}
#if !defined(BOOST_NO_STD_WSTRING)
void validate(boost::any& v, const std::vector<std::wstring>& xs, boost::tribool*, int) {
    using namespace boost::program_options;
    validators::check_first_occurrence(v);
	std::wstring s(validators::get_single_string(xs, true));

    for (size_t i = 0; i < s.size(); ++i)
        s[i] = char(tolower(s[i]));

    if (s.empty() || s == L"on" || s == L"yes" || s == L"1" || s == L"true")
		v = boost::any(boost::tribool(true));
    else if (s == L"off" || s == L"no" || s == L"0" || s == L"false")
		v = boost::any(boost::tribool(false));
    else if (s == L"null" || s == L"maybe" || s == L"2" || s == L"indeterminate")
		v = boost::any(boost::tribool(boost::indeterminate));
    else
        boost::throw_exception(validation_error(validation_error::invalid_option_value));
}
#endif
}

namespace boost { namespace program_options {
template<>
typed_value<bool>* value(bool* v) {
	return bool_switch(v);
}
template<>
typed_value<boost::tribool>* value(boost::tribool* v) {
	typed_value<boost::tribool>* r = new typed_value<boost::tribool>(v);
    r->implicit_value(true, "on");
	return r;
}
}}

namespace dng { namespace app {

namespace base { struct arg_t {
	std::string arg_file;
	bool version;
	bool help;
	std::vector< std::string > input;
		
	void add_to(po::options_description &desc) {
		desc.add_options()
			("help", po::value<bool>(&help)->default_value(false, "off"),
				"display usage information")
			("version", po::value<bool>(&version)->default_value(false, "off"),
				"display version information")
			("arg-file", po::value<std::string>(&arg_file)->default_value(""),
				"read command-line arguments from a file")
		;		
	}
};}

/******************************************************************************
 * class dng::app::Base<Arg>                                                  *
 ******************************************************************************/
template<typename Arg>
class Base {
public:
	typedef Arg arg_t;
	
	arg_t arg;
	
	Base(int argc, char* argv[]) : desc("Allowed Options")  {
		using namespace std;
		runname = argv[0];
		
		arg.add_to(desc);
		indesc.add_options()
			("input", po::value< vector<string> >(&arg.input), "input files")
		;
		indesc.add(desc);
		pdesc.add("input", -1);
		po::store(po::command_line_parser(argc, argv).options(indesc).positional(pdesc).run(), vm);
		po::notify(vm);
		if(!arg.arg_file.empty()) {
			if(arg.arg_file == "-") {
				po::store(po::parse_config_file(cin, desc), vm);	
			} else {
				std::ifstream ifs(arg.arg_file.c_str());
				if(!ifs.is_open()) {
					throw std::runtime_error(
						"unable to open argument file '" + arg.arg_file + "'."
					);
				}
				po::store(po::parse_config_file(ifs, desc), vm);
			}
			po::notify(vm);
		}
	}
	
	int operator()() {
		using namespace std;
		// TODO: Split this up and allow customization
		if(arg.help || arg.input.empty()) {
			//cerr << endl << VERSION_MSG << endl << endl;
			string usage_name(runname);
			if(runname.substr(0,4) == "dng-")
				usage_name[3] = ' ';
			cerr << "Usage:\n  "
				 << usage_name << " [options] input1 input2 input3 ..."
				 << endl << endl;
			cerr << desc << endl;
			return EXIT_SUCCESS;
		}
		if(arg.version) {
			// TODO this
		}
		

		return this->Run();
	}
	
	virtual int Run() = 0;

protected:	
	// TODO: change naming
	std::string runname;
	po::options_description desc, indesc;
	po::positional_options_description pdesc;
	po::variables_map vm;
};

}} // namespace dng::app


#endif // DNG_APP_H

