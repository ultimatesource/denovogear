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

#include <string>
#include <array>
#include <climits>

#include <boost/range/iterator_range_core.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/for_each.hpp>

#include <dng/task/pileup.h>
#include <dng/io/fasta.h>
#include <dng/io/ad.h>
#include <dng/pileup.h>
#include <dng/read_group.h>
#include <dng/seq.h>
#include <dng/utility.h>
#include <dng/depths.h>
#include <dng/matrix.h>

#include "version.h"

using namespace std;
using namespace dng;
using namespace dng::task;

void tab_append(std::string *out, const std::string &in, char d='\t') {
    assert(out != nullptr);
    if(out->empty()) {
        *out = in;
    } else {
        *out += d;
        *out += in;
    }
}

// The main loop for dng-pileup application
// argument_type arg holds the processed command line arguments
int Pileup::operator()(Pileup::argument_type &arg) {
    using namespace std;
    using dng::pileup::AlleleDepths;

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

    // Outputfile
    AlleleDepths line;
    io::Ad output{arg.output, std::ios_base::out};
    string dict = reference.FetchDictionary();
    output.AddHeaderLines(dict);

    // Construct Library Information
    {
        vector<array<string,3>> libs;
        libs.resize(rgs.libraries().size());
        for(size_t u = 0; u < rgs.data().size(); ++u) {
            auto &lib = libs[rgs.library_from_id(u)];
            auto &rg = rgs.data().get<rg::nx>()[u];
            if(!lib[0].empty()) {
                // we have see this library
                if(lib[1] != rg.sample) {
                    throw(runtime_error("@AD ID:" + lib[0] + " is connected to multiple samples."));
                }
                lib[2] += "\tRG:" + rg.id;
                continue;
            }
            lib[0] = rg.library;
            lib[1] = rg.sample;
            lib[2] = "RG:" + rg.id;
        }
        for(auto && l : libs) {
            output.AddLibrary(std::move(l[0]), std::move(l[1]), std::move(l[2]));
        }
    }
    output.WriteHeader();

    // Do pileup
    dng::BamPileup mpileup{rgs.groups(), arg.min_qlen};
    mpileup(bamdata, [&](const dng::BamPileup::data_type & data, utility::location_t loc) {
        // Calculate target position and fetch sequence name
        char ref_base = reference.FetchBase(loc);
        size_t ref_index = seq::char_index(ref_base);

        // reset all depth counters
        read_depths.assign(read_depths.size(), {});

        // construct values for depth sorting
        std::array<unsigned long long,5> total_depths;
        total_depths.fill(0);
        // pileup on read counts
        for(std::size_t u = 0; u < data.size(); ++u) {
            for(auto && r : data[u]) {
                if(filter_read(r)) {
                    continue;
                }
                std::size_t base = seq::base_index(r.aln.seq_at(r.pos));
                total_depths[base] += 1;
                read_depths[rgs.library_from_id(u)].counts[ base ] += 1;
                // detect overflow of our signed number
                assert(read_depths[rgs.library_from_id(u)].counts[ base ] >= 0);
            }
        }
        // shift and pack the nucleotide number at the lowest bits
        // since depths must fit into a signed int, this is safe, check anyways
        static_assert(sizeof(int)*CHAR_BIT+8 <= sizeof(unsigned long long)*CHAR_BIT,
            "Data model may produce unexpected results.");
        for(int i=0;i<total_depths.size();++i) {
            total_depths[i] = (total_depths[i] << 8) | (i & 0xFF);
        }
        // make sure the ref base will come first
        total_depths[ref_index] |= (1ULL << (sizeof(unsigned long long)*CHAR_BIT-1));
        // sort in decreasing order
        boost::sort(total_depths, greater<uint64_t>());
        // find the color of the site
        auto first = total_depths.begin();
        auto last = boost::find_if(total_depths, [](unsigned long long x) { return ((x >> 8) == 0);});
        int first_is_N = 0;
        if(first == last) {
            return; // no data at this site
        } else if( (*first & 0xFF) >= 4 ) {
            first_is_N = 64;
            if(++first == last) {
                return; // no data at this site
            }
        }
        auto rng = boost::make_iterator_range(first,last);
        boost::for_each(rng, [](unsigned long long &x) { x &= 0xFF; });
        int color = AlleleDepths::MatchIndexes(rng)+first_is_N;
        assert(0 <= color && color < 128);


    });

    return EXIT_SUCCESS;
}