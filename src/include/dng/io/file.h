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
        auto file_type = utility::extract_file_type(filename);
        path_ = file_type.path;
        type_label_ = file_type.type_ext;
        compression_ = file_type.compress_ext;

        if(path_.empty()) {
            // don't do anything, just setup type
        } else if(path_ != "-") {
            // if path is not "-" open a file
            buffer_.reset(new std::filebuf);
            std::filebuf *p = static_cast<std::filebuf*>(buffer_.get());
            p->open(path_.c_str(), mode);
            // if file is open, associate it with the stream
            if(p->is_open()) {
                return Attach(buffer_.get());
            }
        } else if((mode & std::ios_base::in) == (mode & std::ios_base::out)) {
            // can't do anything if both are set or none are set
        } else if(mode & std::ios_base::in) {
            return Attach(std::cin.rdbuf());
        } else {
            return Attach(std::cout.rdbuf());
        }
        return Attach(nullptr);
    }
    bool Attach(std::streambuf *buffer) {
        stream_.rdbuf(buffer);
        stream_.unsetf(std::ios_base::skipws);
        return is_open();
    }

    bool is_open() const { return stream_.rdbuf() != nullptr; }
    operator bool() const { return is_open(); }

    std::string path() const { return path_.native(); }

protected:
    std::iostream stream_{nullptr};
    boost::filesystem::path path_;
    std::string type_label_;
    std::string compression_;

private:
    std::unique_ptr<std::streambuf> buffer_;    
};

inline
File::File(File&& other) : 
    buffer_{std::move(other.buffer_)},
    path_{std::move(other.path_)},
    type_label_{std::move(other.type_label_)}
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
