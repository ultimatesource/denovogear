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
#include <fstream>

#include <dng/utility.h>

namespace dng {
namespace io {

class File {
public:
    File() = default;
    
    File(File&& other);

    File& operator=(File&& other);

    explicit File(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        Open(filename, mode);
    }

    bool Open(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        std::tie(type_label_, path_) = utility::extract_file_type(filename);
        buffer_ = nullptr;
        file_.reset();
        if(path_.empty()) {
            // don't do anything, just setup type
        } else if(path_ != "-") {
            // if path is not "-" open a file
            file_.reset(new std::fstream);
            file_.get()->open(path_.c_str(), mode);
            // if file is open, associate it with the stream
            if(file_.get()->is_open()) {
                buffer_ = file_.get()->rdbuf();
            }
        } else if((mode & std::ios_base::in) == (mode & std::ios_base::out)) {
            // can't do anything if both are set or none are set
        } else if(mode & std::ios_base::in) {
            buffer_ = std::cin.rdbuf();
        } else {
            buffer_ = std::cout.rdbuf();
        }
        return Attach();
    }
    bool Attach(std::streambuf *buffer) {
        is_open_ = (buffer != nullptr);
        stream_.rdbuf(buffer);
        stream_.unsetf(std::ios_base::skipws);
        return is_open_;
    }
    bool Attach() {
        return Attach(buffer_);
    }

    bool is_open() const { return is_open_; }
    operator bool() const { return is_open(); }

    std::string path() const { return path_.native(); }

protected:
    std::iostream stream_{nullptr};
    std::streambuf *buffer_{nullptr};
    boost::filesystem::path path_;
    std::string type_label_;

private:
    bool is_open_{false};
    
    std::unique_ptr<std::fstream> file_;
};

inline
File::File(File&& other) : 
    buffer_{std::move(other.buffer_)},
    path_{std::move(other.path_)},
    type_label_{std::move(other.type_label_)},
    file_{std::move(other.file_)}
{
    std::streambuf *buffer = other.stream_.rdbuf();
    Attach(buffer);
    other.Attach(nullptr);
} 

inline
File& File::operator=(File&& other) {
    if(this == &other) {
        return *this;
    }

    buffer_ = std::move(other.buffer_);
    path_ = std::move(other.path_);
    type_label_ = std::move(other.type_label_);

    file_ = std::move(other.file_);

    std::streambuf *buffer = other.stream_.rdbuf();
    Attach(buffer);
    other.Attach(nullptr);

    return *this;
}

class BinaryFile : public File {
public:
    BinaryFile() = default;

    BinaryFile(BinaryFile&&) = default;
    BinaryFile& operator=(BinaryFile&&) = default;
    
    explicit BinaryFile(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        Open(filename, mode);
    }

    void Open(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        File::Open(filename, mode | std::ios::binary);
    }
};

} // namespace io
} // namespace dng

#endif // DNG_IO_FILE_H
