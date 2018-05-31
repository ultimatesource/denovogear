/*
 * Copyright (c) 2014-2018 Reed A. Cartwright
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
#include <boost/algorithm/string/trim.hpp>

namespace dng { namespace utility {

// extracts extension and filename from both file.ext and ext:file.foo
// trims whitespace as well
file_type_t extract_file_type(std::string path) {
    using boost::algorithm::trim;
    using boost::algorithm::ends_with;
    const auto &npos = std::string::npos;

    auto is_compressed = [](const std::string &path) -> std::string {
        for(auto && s : {".gz", ".gzip", ".bgz"}) {
            if(ends_with(path, s)) {
                return {s+1};
            }
        }
        return {};
    };

    // Remove whitespace
    trim(path);

    // Format ext:path ???
    auto colon = path.find_first_of(':');
    if(colon != npos) {
        auto ext = path.substr(0, colon);
        path.erase(0,colon+1);
        // Format ext.gz:path ???
        auto gz = is_compressed('.'+ext);
        if(!gz.empty()) {
            // Format gz:path
            if(gz.size() == ext.size()) {
                ext.erase();
            } else {
                auto gz_pos = ext.size() - (gz.size()+1);
                ext.erase(gz_pos);
            }
        }
        return {path, ext, gz};
    }
    // Format path.gz ???
    auto gz = is_compressed(path);
    auto gz_pos = path.size() - (gz.empty() ? 0 : gz.size()+1);
    auto ext_pos = path.find_last_of('.', gz_pos-1);
    if(ext_pos == npos || ext_pos == 0) {
        return {path, {}, gz};
    }
    auto ext = path.substr(ext_pos+1, gz_pos-(ext_pos+1));
    return {path, ext, gz};
}

static const std::string file_category_keys[] = {
    "bam","sam","cram",
    "bcf","vcf"
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
    default:
    	break;
    };
    return FileCat::Unknown;
}

FileCat input_category(const std::string &in, FileCatSet mask, FileCat def) {
    std::string ext = extract_file_type(in).type_ext;

    FileCat cat = file_category(ext);
    if(cat == FileCat::Unknown)
        cat = def;
    if(mask & cat) {
        return cat;
    } else {
        throw std::invalid_argument("file type '" + ext + "' not supported. Input file was '" + in + "'.");
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
