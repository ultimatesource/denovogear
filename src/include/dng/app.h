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

#pragma once
#ifndef DNG_APP_H
#define DNG_APP_H

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/filesystem.hpp>

#include <boost/logic/tribool.hpp>
#include <boost/logic/tribool_io.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "version.h"

// Support the use of tribools in boost::program_options
namespace boost {
void validate(boost::any &v, const std::vector<std::string> &xs,
              boost::tribool *, int) {
    using namespace boost::program_options;
    validators::check_first_occurrence(v);
    std::string s(validators::get_single_string(xs, true));

    for(size_t i = 0; i < s.size(); ++i) {
        s[i] = char(tolower(s[i]));
    }

    if(s.empty() || s == "on" || s == "yes" || s == "1" || s == "true") {
        v = boost::any(boost::tribool(true));
    } else if(s == "off" || s == "no" || s == "0" || s == "false") {
        v = boost::any(boost::tribool(false));
    } else if(s == "null" || s == "maybe" || s == "2" || s == "indeterminate") {
        v = boost::any(boost::tribool(boost::indeterminate));
    } else {
        boost::throw_exception(validation_error(validation_error::invalid_option_value,
                                                s));
    }
}
#if !defined(BOOST_NO_STD_WSTRING)
void validate(boost::any &v, const std::vector<std::wstring> &xs,
              boost::tribool *, int) {
    using namespace boost::program_options;
    validators::check_first_occurrence(v);
    std::wstring s(validators::get_single_string(xs, true));

    for(size_t i = 0; i < s.size(); ++i) {
        s[i] = char(tolower(s[i]));
    }

    if(s.empty() || s == L"on" || s == L"yes" || s == L"1" || s == L"true") {
        v = boost::any(boost::tribool(true));
    } else if(s == L"off" || s == L"no" || s == L"0" || s == L"false") {
        v = boost::any(boost::tribool(false));
    } else if(s == L"null" || s == L"maybe" || s == L"2" || s == L"indeterminate") {
        v = boost::any(boost::tribool(boost::indeterminate));
    } else {
        boost::throw_exception(validation_error(
                                   validation_error::invalid_option_value));
    }
}
#endif
}

namespace boost {
namespace program_options {
template<>
typed_value<bool> *value(bool *v) {
    return bool_switch(v);
}
template<>
typed_value<boost::tribool> *value(boost::tribool *v) {
    typed_value<boost::tribool> *r = new typed_value<boost::tribool>(v);
    r->implicit_value(true, "on");
    return r;
}
}
}

namespace dng {

namespace detail {
template<class T>
void add_app_args(po::options_description &desc, T &arg) {
    return;
}

}

template<typename T>
class CommandLineApp {
public:
    typedef T task_type;

    struct arg_t : public task_type::argument_type {
        bool help;
        bool version;
        std::string arg_file;

        std::string run_name;
        std::string run_path;
    } arg;

    CommandLineApp(int argc, char *argv[]) : ext_desc_("Allowed Options")  {
        using namespace std;
        using detail::add_app_args;

        boost::filesystem::path bin_path(argv[0]);
        arg.run_name = bin_path.filename().generic_string();
        arg.run_path = bin_path.parent_path().generic_string();

        add_app_args(ext_desc_, static_cast<typename task_type::argument_type &>(arg));

        ext_desc_.add_options()
        ("version", po::value<bool>(&arg.version)->default_value(false, "off"),
         "display version information")
        ("help", po::value<bool>(&arg.help)->default_value(false, "off"),
         "display usage informaiton")
        ("arg-file", po::value<std::string>(&arg.arg_file)->default_value(""),
         "read command-line arguments from a file")
        ;

        int_desc_.add_options()
        ("input", po::value< vector<string> >(&arg.input), "input files")
        ;
        int_desc_.add(ext_desc_);
        pos_desc_.add("input", -1);

        po::store(po::command_line_parser(argc, argv)
                  .options(int_desc_).positional(pos_desc_).run(), vm_);
        po::notify(vm_);

        if(!arg.arg_file.empty()) {
            if(arg.arg_file == "-") {
                po::store(po::parse_config_file(cin, ext_desc_), vm_);
            } else {
                std::ifstream ifs(arg.arg_file.c_str());
                if(!ifs.is_open()) {
                    throw std::runtime_error(
                        "unable to open argument file '" + arg.arg_file + "'."
                    );
                }
                po::store(po::parse_config_file(ifs, ext_desc_), vm_);
            }
            po::notify(vm_);
        }
    }

    int operator()() {
        using namespace std;
        if(arg.version) {
            return CmdVersion();
        }
        if(arg.help || arg.input.empty()) {
            return CmdHelp();
        }
        return task_(arg);
    }

protected:
    virtual int CmdHelp() const {
        using namespace std;
        string usage_name(arg.run_name);
        if(usage_name.substr(0, 4) == "dng-") {
            usage_name[3] = ' ';
        }
        cerr << "Usage:\n  "
             << usage_name << " [options] input1 input2 input3 ..."
             << endl << endl;
        cerr << ext_desc_ << endl;
        return EXIT_SUCCESS;
    }
    virtual int CmdVersion() const {
        using namespace std;
        string usage_name(arg.run_name);
        if(usage_name.substr(0, 4) == "dng-") {
            usage_name[3] = ' ';
        }

        cerr << usage_name << " v" PACKAGE_VERSION << "\n";
        cerr << "Copyright (c) 2014-2015 Reed A. Cartwright, Kael Dai, et al.\n";
        return EXIT_SUCCESS;
    }

    po::options_description ext_desc_, int_desc_;
    po::positional_options_description pos_desc_;
    po::variables_map vm_;

    task_type task_;
};

} // namespace dng


#endif // DNG_APP_H

