/*
 * Copyright (c) 2016-2017 Reed A. Cartwright
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
#include <dng/io/ad.h>
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
int process_bam(LogLike::argument_type &arg);
int process_ad(LogLike::argument_type &arg);
int process_bcf(LogLike::argument_type &arg);

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
    FileCat mode = utility::input_category(*it, FileCat::Sequence|FileCat::Pileup|FileCat::Variant, FileCat::Sequence);
    for(++it; it != arg.input.end(); ++it) {
        if(utility::input_category(*it, FileCat::Sequence|FileCat::Pileup|FileCat::Variant, FileCat::Sequence) != mode) {
            throw std::invalid_argument("Mixing sam/bam/cram, vcf/bcf, and tad/ad input files is not supported.");
        }
    }
    // Execute sub tasks based on input type
    if(mode == FileCat::Pileup) {
        return process_ad(arg);
    } else if(mode == FileCat::Variant) {
        // vcf, bcf
        return process_bcf(arg);
    } else if(mode == FileCat::Sequence) {
        return process_bam(arg);
    } else {
        throw std::invalid_argument("Unknown input data file type.");
    }
    return EXIT_FAILURE;
}

void output_loglike_results(std::ostream &o, double hidden, double observed) {
    // output results
    o << setprecision(std::numeric_limits<double>::max_digits10)
         << "log_likelihood\t" << hidden+observed << "\n";
    o << setprecision(std::numeric_limits<double>::max_digits10)
         << "log_hidden\t" << hidden << "\n";
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

    LogProbability do_loglike(relationship_graph, get_model_parameters(arg));

    // Pileup data
    dng::pileup::allele_depths_t read_depths(make_array(mpileup.num_libraries(), 5u));

    const size_t num_nodes = relationship_graph.num_nodes();
    const size_t library_start = relationship_graph.library_nodes().first;
    const int min_basequal = arg.min_basequal;
    auto filter_read = [min_basequal](
    decltype(mpileup)::data_type::value_type::const_reference r) -> bool {
        return (r.is_missing
        || r.base_qual() < min_basequal
        || seq::base_index(r.base()) >= 4);
    };

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

        // reset all depth counters
        std::fill_n(read_depths.data(), read_depths.num_elements(), 0);

        // pileup on read counts
        for(std::size_t u = 0; u < data.size(); ++u) {
            for(auto && r : data[u]) {
                if(filter_read(r)) {
                    continue;
                }
                std::size_t base = seq::base_index(r.base());
                assert(read_depths[u][(ref_index+base) % 5] < 65535);

                read_depths[u][(ref_index+base) % 5] += 1;
            }
        }
        int num_alts = (ref_index < 4) ? 3 : 4;
        auto loglike = do_loglike(read_depths, num_alts, ref_index < 4);
        sum_data += loglike.log_data;
        sum_scale += loglike.log_scale;
    });

    output_loglike_results(cout, sum_data.result(), sum_scale.result());

    return EXIT_SUCCESS;
}

int process_bcf(LogLike::argument_type &arg) {
    // Read input data
    auto mpileup = io::BcfPileup::open_and_setup(arg);

    auto relationship_graph = create_relationship_graph(arg, &mpileup);

    // Read header from first file
    const bcf_hdr_t *header = mpileup.reader().header(0); // TODO: fixthis
    const int num_libs = mpileup.num_libraries();

    LogProbability do_loglike (relationship_graph, get_model_parameters(arg));

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

        auto loglike = do_loglike(read_depths,n_sz,true);
        sum_data += loglike.log_data;
        sum_scale += loglike.log_scale;
    });

    // output results
    output_loglike_results(cout, sum_data.result(), sum_scale.result());

    return EXIT_SUCCESS;
}

int process_ad(LogLike::argument_type &arg) {
    using namespace std;
    // Parse pedigree from file
    Pedigree ped = io::parse_ped(arg.ped);

    // Open Reference
    io::Fasta reference{arg.fasta.c_str()};

    // replace arg.region with the contents of a file if needed
    //auto region_ext = io::at_slurp(arg.region);
    if(!arg.region.empty()) {
        throw std::invalid_argument("--region not supported when processing ad/tad file.");
    }

    // Open input files
    if(arg.input.size() != 1) {
        throw std::runtime_error("can only process one ad/tad file at a time.");
    }
    
    io::Ad input{arg.input[0], std::ios_base::in};
    if(!input) {
        throw std::runtime_error("unable to open input file '" + arg.input[0] + "'.");
    }
    if(input.ReadHeader() == 0) {
        throw std::runtime_error("unable to read header from '" + input.path() + "'.");
    }

    // Construct peeling algorithm from parameters and pedigree information
    InheritanceModel model = inheritance_model(arg.model);

    RelationshipGraph relationship_graph;
    if(!relationship_graph.Construct(ped, input.libraries(), model, arg.mu, arg.mu_somatic, arg.mu_library,
        arg.normalize_somatic_trees)) {
        throw std::invalid_argument("Unable to construct peeler for pedigree; "
                            "possible non-zero-loop pedigree.");
    }
    // Select libraries in the input that are used in the pedigree
    input.SelectLibraries(relationship_graph.library_names());

    for(auto && line : relationship_graph.BCFHeaderLines()) {
        cout << line << "\n";
    }

    // Construct function object to calculate log likelihoods
    LogProbability::params_t params = get_model_parameters(arg);

    LogProbability calculate (relationship_graph, params);

    stats::ExactSum sum_data, sum_scale;

    // using a single thread to read data and process it
    if(arg.threads == 0) {
        pileup::AlleleDepths line;
        line.data().reserve(4*input.num_libraries());
        dng::pileup::allele_depths_t read_depths(make_array(line.num_libraries(),5u));

        // read each line of data into line and process it
        while(input.Read(&line)) {
            int ref_is_N = (line.color() >= 64);
            size_t sz = line.num_nucleotides() + ref_is_N;
            read_depths.resize(make_array(line.num_libraries(),sz));
            for(int i=0;i<line.num_libraries();++i) {
                for(int a=ref_is_N;a<sz;++a) {
                    read_depths[i][a] = line(i,a-ref_is_N);
                }
            }

            auto loglike = calculate(read_depths, sz-1, !ref_is_N);
            sum_data += loglike.log_data;
            sum_scale += loglike.log_scale;
        }
    } else {
        // Use a single thread to read data from the input
        // and multiple threads to process it.
        // Data will be processed in batches and a finite
        // number of batches will be shared between the main
        // and worker threads. std::move will be used to
        // efficiently copy the batches between threads.

        // force batch_size to be positive
        if(arg.batch_size <= 0) {
            arg.batch_size = 10000;
        }
        // construct vectors to hold our batches
        size_t num_batches = arg.threads+3;
        size_t batch_size = arg.batch_size;
        typedef vector<pileup::AlleleDepths> batch_t;
        stack<batch_t> batches;
        for(size_t u=0; u<num_batches; ++u) {
            batches.emplace(batch_size,pileup::AlleleDepths{});
        }

        // synchronization objects for the batch stack
        mutex batch_mutex;
        condition_variable batch_cond;

        // construct a lambda function which will be used for the worker threads
        // each thread needs to have its own, mutable copy of calculate
        // because calculate stores working data in a member object
        auto batch_calculate = [&,calculate](batch_t& reads) mutable {
                stats::ExactSum sum_d, sum_s;
                dng::pileup::allele_depths_t read_depths;

                for(auto && line : reads) {
                    int ref_is_N = (line.color() >= 64);
                    size_t sz = line.num_nucleotides() + ref_is_N;
                    read_depths.resize(make_array(line.num_libraries(),sz));
                    for(int i=0;i<line.num_libraries();++i) {
                        for(int a=ref_is_N;a<sz;++a) {
                            read_depths[i][a] = line(i,a-ref_is_N);
                        }
                    }
                    auto loglike = calculate(read_depths, sz-1, !ref_is_N);
                    sum_d += loglike.log_data;
                    sum_s += loglike.log_scale;
            }
            {
                // use batch_mutex to save results and return our data object
                lock_guard<mutex> lock(batch_mutex);
                sum_data += sum_d;
                sum_scale += sum_s;
                batches.emplace(std::move(reads));
            }
            // notify the reader thread that we have returned our object
            batch_cond.notify_one();
        };

        // construct worker thread pool to run batch_calculate
        multithread::BasicPool<batch_t&> worker_pool(batch_calculate,arg.threads);

        // loop until we hit the end-of-file
        bool eof = false;
        while(!eof) {
            // setup a batch
            batch_t b;
            {
                // try to pull a batch object off of the stack
                std::unique_lock<std::mutex> lock(batch_mutex);
                batch_cond.wait(lock, [&](){return !batches.empty();});
                b = std::move(batches.top());
                batches.pop();
            }
            // resize the batch
            b.resize(batch_size);
            // read up to batch_size number of reads into the batch
            // resize if needed
            for(size_t u=0; u < batch_size; ++u) {
                if(!input.Read(&b[u])) {
                    b.resize(u);
                    eof = true;
                    break;
                }
            }
            // schedule this batch to run
            worker_pool.Enqueue(std::move(b));
        }
    }

    // output results
    output_loglike_results(cout, sum_data.result(), sum_scale.result());

    return EXIT_SUCCESS;
}
