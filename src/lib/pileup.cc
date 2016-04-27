/*
 * Copyright (c) 2014-2015 Reed A. Cartwright
 * Copyright (c) 2015 Kael Dai
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *           Kael Dai <kdai1@asu.edu>
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

#include <cstdlib>
#include <fstream>
#include <iterator>
#include <iosfwd>

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <sstream>
#include <string>
#include <array>

#include <boost/range/iterator_range.hpp>
#include <boost/range/algorithm/replace.hpp>
#include <boost/range/algorithm/max_element.hpp>

#include <boost/algorithm/string.hpp>

#include <dng/task/pileup.h>
#include <dng/io/fasta.h>
#include <dng/pedigree.h>
#include <dng/fileio.h>
#include <dng/pileup.h>
#include <dng/read_group.h>
#include <dng/seq.h>
#include <dng/utility.h>
#include <dng/hts/bcf.h>
#include <dng/hts/extra.h>
#include <dng/vcfpileup.h>

#include "version.h"

using namespace std;
using namespace dng;
using namespace dng::task;

// The main loop for dng-pileup application
// argument_type arg holds the processed command line arguments
int Pileup::operator()(Pileup::argument_type &arg) {
    using namespace std;

    // Open Reference
    io::Fasta reference{arg.fasta.c_str()};

    // quality thresholds
    int min_qual = arg.min_basequal;

    // Open input files
    vector<hts::bam::File> bamdata;
    for(auto && str : arg.input) {
        bamdata.emplace_back(str.c_str(), "r");
        if(bamdata.back().is_open()) {
            continue;
        }
        throw std::runtime_error("unable to open bam/sam/cram file '" + str + "' for reading.");
    }

    // replace arg.region with the contents of a file if needed
    io::at_slurp(arg.region);
    // parse region
    if(!arg.region.empty()) {
        for(auto && f : bamdata) {
            auto r = regions::bam_parse(arg.region, f);
            f.regions(std::move(r));
        }
    }    

    // Add each genotype/sample column
    dng::ReadGroups rgs;
    rgs.ParseHeaderText(bamdata, arg.rgtag);

    // Pileup data
    std::vector<depth_t> read_depths(rgs.libraries().size());

    const int min_basequal = arg.min_basequal;
    auto filter_read = [min_basequal](
    dng::BamPileup::data_type::value_type::const_reference r) -> bool {
        return (r.is_missing
        || r.qual.first[r.pos] < min_basequal
        || seq::base_index(r.aln.seq_at(r.pos)) >= 4);
    };

    // Read header from first file
    dng::BamPileup mpileup{rgs.groups(), arg.min_qlen};
    mpileup(bamdata, [&](const dng::BamPileup::data_type & data, utility::location_t loc) {
        // Calculate target position and fetch sequence name
        char ref_base = reference.FetchBase(loc);
        size_t ref_index = seq::char_index(ref_base);

        // reset all depth counters
        read_depths.assign(read_depths.size(), {});

        // construct values for depth sorting
        typedef array<int, 3> key_t;
        key_t total_depths[5] = {{0,0,0}, {0,0,1}, {0,0,2}, {0,0,3}, {0,0,4}};
        // make sure the ref base will come first
        total_depths[ref_index][0] = 1;

        // pileup on read counts
        for(std::size_t u = 0; u < data.size(); ++u) {
            for(auto && r : data[u]) {
                if(filter_read(r)) {
                    continue;
                }
                std::size_t base = seq::base_index(r.aln.seq_at(r.pos));
                total_depths[base][1] += 1;
                read_depths[rgs.library_from_id(u)].counts[ base ] += 1;
                // detect overflow of our signed number
                assert(read_depths[rgs.library_from_id(u)].counts[ base ] >= 0);
            }
        }
        sort(&total_depths[0],&total_depths[5],greater<key_t>());


    });

    return EXIT_SUCCESS;
}