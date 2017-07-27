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
#include <string>

#include <boost/logic/tribool.hpp>
#include <boost/logic/tribool_io.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace dng {

namespace task {

struct bool_switch_t {
    bool value;

    explicit bool_switch_t(bool b) : value(b) { }
    bool_switch_t() { }

    operator bool() {
        return value;
    }
    bool_switch_t& operator =(bool b) {
        value = b;
        return *this;
    }
};

inline
std::ostream& operator<<(std::ostream& os, bool_switch_t b) {
    os << b.value;
    return os;
}

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

namespace dng {
namespace task {
// Support the use of bool_switch_t and tribool in boost::program_options
inline
void validate(boost::any &v, const std::vector<std::string> &xs,
              bool_switch_t *, int) {
    boost::program_options::validate(v, xs, (bool*)0, 0);
    bool b = boost::any_cast<bool>(v);
    v = boost::any(bool_switch_t{b});
}

#if !defined(BOOST_NO_STD_WSTRING)
inline
void validate(boost::any &v, const std::vector<std::wstring> &xs,
              bool_switch_t *, int) {
    boost::program_options::validate(v, xs, (bool*)0, 0);
    bool b = boost::any_cast<bool>(v);
    v = boost::any(bool_switch_t{b});
}
#endif
}
}

namespace boost {
inline
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
inline
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
inline
typed_value<::dng::task::bool_switch_t>*
value(::dng::task::bool_switch_t *v) {
    using bool_switch_t = ::dng::task::bool_switch_t;  
    typed_value<bool_switch_t>* r = new typed_value<bool_switch_t>(v);
    r->default_value(bool_switch_t{false});
    r->zero_tokens();
    return r;
}

template<>
inline
typed_value<bool>* value(bool *v) {
    typed_value<bool>* r = new typed_value<bool>(v);
    r->implicit_value(true, "on");
    return r;
}

template<>
inline
typed_value<boost::tribool>* value(boost::tribool *v) {
    typed_value<boost::tribool>* r = new typed_value<boost::tribool>(v);
    r->implicit_value(true, "on");
    return r;
}
}
}

#endif // DNG_TASK_H

