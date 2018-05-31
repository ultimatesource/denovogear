/*
 * Copyright (c) 2016-2018 Reed A. Cartwright
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
// Workaround for using boost::fusion::invoke and c++11 lambdas before boost 1.58
#include <boost/version.hpp>
#if BOOST_VERSION < 105800
#   define BOOST_RESULT_OF_USE_TR1_WITH_DECLTYPE_FALLBACK 1
#endif

#include <cstdlib>
#include <fstream>
#include <iterator>
#include <iosfwd>

#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <string>
#include <vector>
#include <stack>
#include <numeric>

#include <boost/range/iterator_range.hpp>
#include <boost/range/algorithm/replace.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/replace_if.hpp>

#include <boost/algorithm/string.hpp>

#include <dng/task/loglike.h>
#include <dng/probability.h>
#include <dng/fileio.h>
#include <dng/seq.h>
#include <dng/utility.h>
#include <dng/hts/bcf.h>
#include <dng/mutation.h>
#include <dng/stats.h>
#include <dng/io/utility.h>
#include <dng/io/fasta.h>
#include <dng/io/ped.h>
#include <dng/io/bam.h>
#include <dng/io/bcf.h>
#include <dng/depths.h>
#include <dng/multithread.h>

#include <htslib/faidx.h>
#include <htslib/khash.h>

#include "version.h"

using namespace std;
using namespace dng;
using namespace dng::task;

using utility::make_array;

// Sub-tasks
namespace {
int process_bam(LogLike::argument_type &arg);
int process_bcf(LogLike::argument_type &arg);
} // anon namespace

// The main loop for dng-loglike application
// argument_type arg holds the processed command line arguments
int task::LogLike::operator()(task::LogLike::argument_type &arg) {
    using utility::FileCat;
    using utility::FileCatSet;
    // if input is empty default to stdin.
    if(arg.input.empty()) {
        arg.input.emplace_back("-");
    }

    // Check that all input formats are of same category
    auto it = arg.input.begin();
    FileCat mode = utility::input_category(*it, FileCat::Sequence|FileCat::Variant, FileCat::Sequence);
    for(++it; it != arg.input.end(); ++it) {
        if(utility::input_category(*it, FileCat::Sequence|FileCat::Variant, FileCat::Sequence) != mode) {
            throw std::invalid_argument("Mixing sam/bam/cram and vcf/bcf input files is not supported.");
        }
    }
    // Execute sub tasks based on input type
    if(mode == FileCat::Variant) {
        // vcf, bcf
        return process_bcf(arg);
    } else if(mode == FileCat::Sequence) {
        return process_bam(arg);
    } else {
        throw std::invalid_argument("Unknown input data file type.");
    }
    return EXIT_FAILURE;
}

namespace {

void output_loglike_results(std::ostream &o, double total, double observed) {
    // output results
    o << setprecision(std::numeric_limits<double>::max_digits10)
         << "log_likelihood\t" << total << "\n";
    o << setprecision(std::numeric_limits<double>::max_digits10)
         << "log_hidden\t" << (total-observed) << "\n";
    o << setprecision(std::numeric_limits<double>::max_digits10)
         << "log_observed\t" << observed << "\n";
}

int process_bam(LogLike::argument_type &arg) {
    // Open Reference
    if(arg.fasta.empty()){
        throw std::invalid_argument("Path to reference file must be specified with --fasta when processing bam/sam/cram files.");
    }
    io::Fasta reference{arg.fasta.c_str()};

    // Open input files
    auto mpileup = io::BamPileup::open_and_setup(arg);

    auto relationship_graph = create_relationship_graph(arg, &mpileup);

    Probability model{relationship_graph, get_model_parameters(arg)};

    const size_t num_nodes = relationship_graph.num_nodes();
    const size_t library_start = relationship_graph.library_nodes().first;
    const int min_basequal = arg.min_basequal;
    auto filter_read = [min_basequal](
    decltype(mpileup)::data_type::value_type::const_reference r) -> bool {
        return (r.is_missing
        || r.base_qual() < min_basequal
        || seq::base_index(r.base()) >= 4);
    };

    decltype(mpileup)::Alleles count_alleles(mpileup.num_libraries());

    // Treat sequence_data and variant data separately
    dng::stats::ExactSum sum_data;
    dng::stats::ExactSum sum_scale;
    auto h = mpileup.header();
    
    mpileup([&](const decltype(mpileup)::data_type & data, utility::location_t loc) {
        // Calculate target position and fetch sequence name
        int contig = utility::location_to_contig(loc);
        int position = utility::location_to_position(loc);

        // Calculate reference base
        assert(0 <= contig && contig < h->n_targets);
        char ref_base = reference.FetchBase(h->target_name[contig],position);
        size_t ref_index = seq::char_index(ref_base);

        auto read_depths = count_alleles(data, ref_index, filter_read);
        size_t n_sz = read_depths.shape()[1];
        if(n_sz == 0) {
            return;
        }

        double loglike = model.CalculateLLD(read_depths, n_sz);
        sum_data += loglike;
        sum_scale += model.work().ln_scale;
    });

    output_loglike_results(cout, sum_data.result(), sum_scale.result()/M_LN10);

    return EXIT_SUCCESS;
}

int process_bcf(LogLike::argument_type &arg) {
    // Read input data
    auto mpileup = io::BcfPileup::open_and_setup(arg);

    auto relationship_graph = create_relationship_graph(arg, &mpileup);

    // Read header from first file
    const bcf_hdr_t *header = mpileup.reader().header(0); // TODO: fixthis
    const int num_libs = mpileup.num_libraries();

    Probability model{relationship_graph, get_model_parameters(arg)};

    const size_t num_nodes = relationship_graph.num_nodes();
    const size_t library_start = relationship_graph.library_nodes().first;

    // allocate space for ad. bcf_get_format_int32 uses realloc internally
    int n_ad_capacity = num_libs*5;
    auto ad = hts::bcf::make_buffer<int>(n_ad_capacity);

    // Treat sequence_data and variant data separately
    dng::stats::ExactSum sum_data;
    dng::stats::ExactSum sum_scale;
    
    // run calculation based on the depths at each site.
    mpileup([&](const decltype(mpileup)::data_type & rec) {
        // Read all the Allele Depths for every sample into an AD array
        const int n_ad = hts::bcf::get_format_int32(header, rec, "AD", &ad, &n_ad_capacity);
        if(n_ad <= 0) {
            // AD tag is missing, so we do nothing at this time
            // TODO: support using calculated genotype likelihoods
            return;
        }
        const int n_sz = n_ad / num_libs;
        assert(n_ad % num_libs == 0);

        // replace "missing" and "end" values with 0
        boost::replace_if(boost::make_iterator_range(ad.get(), ad.get()+n_ad),
            [](int a) { return a < 0; }, 0 );

        pileup::allele_depths_ref_t read_depths(ad.get(), make_array(num_libs,n_sz));

        double loglike = model.CalculateLLD(read_depths,n_sz);
        sum_data += loglike;
        sum_scale += model.work().ln_scale;
    });

    // output results
    output_loglike_results(cout, sum_data.result(), sum_scale.result()/M_LN10);

    return EXIT_SUCCESS;
}

} // anon namespce
