/*
 * Copyright (c) 2014-2016 Reed A. Cartwright
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

#include <dng/utility.h>

#include <boost/filesystem/convenience.hpp>

namespace dng { namespace utility {

// extracts extension and filename from both file.ext and ext:file.foo
// returns {ext, filename.ext}
// trims whitespace as well
std::pair<std::string, std::string> extract_file_type(const std::string &path) {
    if(path.empty())
        return {};
    std::locale loc;

    auto last = path.length();
    decltype(last) first = 0;
    while(first < path.length() && std::isspace(path[first],loc)) {
        ++first;
    }
    if(first == last) {
        return {};
    }
    while(std::isspace(path[last-1],loc)) {
        --last;
    }

    auto x = last-1;
    for(auto u = first; u < last; ++u) {
        if(path[u] == ':' && u > first+1) { // u > 1 skips windows drive letters
            return {path.substr(first, u-first), path.substr(u+1,last-(u+1))};
        }
        if(path[u] == '.' && u > first) { // u > 0 skips unix .hidden files
            x = u;
        }
    }
    return {path.substr(x + 1, last - (x + 1)), path.substr(first,last-first)};
}

static const std::string file_category_keys[] = {
    "bam","sam","cram",
    "bcf","vcf",
    "ad","tad"
};

FileCat file_category(const std::string &ext) {
    switch(key_switch_iequals(ext, file_category_keys)) {
    case 0:
    case 1:
    case 2:
        return FileCat::Sequence;
    case 3:
    case 4:
        return FileCat::Variant;
    case 5:
    case 6:
    	return FileCat::Pileup;
    default:
    	break;
    };
    return FileCat::Unknown;
}


FileCat input_category(const std::string &in, FileCatSet mask, FileCat def) {

    std::string ext = extract_file_type(in).first;
    if(ext == "gz" || ext == "gzip" || ext == "bgz") {
    	std::string extr = boost::filesystem::change_extension(in, "").string();
    	ext = extract_file_type(extr).first;
    }

    FileCat cat = file_category(ext);
    if(cat == FileCat::Unknown)
        cat = def;
    if(mask & cat) {
        return cat;
    } else {
        throw std::runtime_error("Argument error: file type '" + ext + "' not supported. Input file was '" + in + "'.");
    }
    return FileCat::Unknown;
}

std::pair<std::string, std::string> timestamp() {
    using namespace std;
    using namespace std::chrono;
    std::string buffer(127, '\0');
    auto now = system_clock::now();
    auto now_t = system_clock::to_time_t(now);
    size_t sz = strftime(&buffer[0], 127, "%FT%T%z",
                         localtime(&now_t));
    buffer.resize(sz);
    auto epoch = std::chrono::duration_cast<std::chrono::milliseconds>(
                     now.time_since_epoch());
    return {buffer, to_string(epoch.count())};
}

std::string vcf_timestamp() {
    auto stamp = timestamp();
    return "Date=" + stamp.first + ",Epoch=" + stamp.second;
}

} }
