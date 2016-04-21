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

#pragma once
#ifndef DNG_IO_FILE_H
#define DNG_IO_FILE_H

#include <boost/filesystem.hpp>
#include <tuple>
#include <iostream>

#include <dng/utility.h>

namespace dng {
namespace io {

class BinaryFile {
public:
    BinaryFile() { }

    void Open(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        mode |= std::ios::binary;
        std::tie(type_label_, path_) = utility::extract_file_type(filename);
        buffer_ = nullptr;
        if(path_.empty()) {
            // don't do anything, just setup type
        } else if(path_ != "-") {
            // if path is not "-" open a file
            file_.open(path_, mode);
            // if file is open, associate it with the stream
            if(file_.is_open()) {
                buffer_ = file_.rdbuf();
            }
        } else if((mode & std::ios_base::in) == (mode & std::ios_base::out)) {
            // can't do anything if both are set or none are set
        } else if(mode & std::ios_base::in) {
            buffer_ = std::cin.rdbuf();
        } else {
            buffer_ = std::cout.rdbuf();
        }
        Attach();
    }
    void Attach(std::streambuf *buffer) {
        is_open_ = (buffer != nullptr);
        stream_.rdbuf(buffer);
    }
    void Attach() {
        Attach(buffer_);
    }

    bool is_open() const { return is_open_; };
    operator bool() const { return is_open(); }

protected:
    std::iostream stream_{nullptr};
    std::streambuf *buffer_{nullptr};
    boost::filesystem::path path_;
    std::string type_label_;

private:
    bool is_open_{false};
    boost::filesystem::fstream file_;
};

} // namespace io
} // namespace dng

#endif // DNG_IO_FILE_H
