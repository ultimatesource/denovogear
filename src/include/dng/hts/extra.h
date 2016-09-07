/*
 * Copyright (c) 2015 Reed A. Cartwright
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
#ifndef CXX_HTS_EXTRA_H
#define CXX_HTS_EXTRA_H

#include <string>
#include <algorithm>

#include "htslib/sam.h"
#include "htslib/vcf.h"

namespace hts {
namespace extra {

// VCF header lacks a function to get sequence lengths
// So we will extract the contig lines from the input header
inline std::vector<std::string> extract_contigs(const bcf_hdr_t *hdr) {
    if(hdr == nullptr)
        return {};
    // Read text of header
    int len;
    std::unique_ptr<char[], void(*)(void *)> str{bcf_hdr_fmt_text(hdr, 0, &len), free};
    if(!str)
        return {};
    std::vector<std::string> contigs;

    // parse ##contig lines
    const char *text = str.get();
    if(strncmp(text, "##contig=", 9) != 0) {
        text = strstr(text, "\n##contig=");
    } else {
        text = text - 1;
    }
    const char *end;
    for(; text != nullptr; text = strstr(end, "\n##contig=")) {
        for(end = text + 10; *end != '\n' && *end != '\0'; ++end)
            /*noop*/;
        if(*end != '\n') {
            return contigs;    // bad header, return what we have.
        }
        contigs.emplace_back(text + 1, end);
    }

    return contigs;
}

// Build a list of all of the possible contigs to add to the vcf header
inline std::vector<std::pair<std::string, uint32_t>> parse_contigs(
const bam_hdr_t *hdr) {
    if(hdr == nullptr)
        return {};
    std::vector<std::pair<std::string, uint32_t>> contigs;
    uint32_t n_targets = hdr->n_targets;
    for(size_t a = 0; a < n_targets; a++) {
        if(hdr->target_name[a] == nullptr) {
            continue;
        }
        contigs.emplace_back(hdr->target_name[a], hdr->target_len[a]);
    }
    return contigs;
}


}
}

#endif
