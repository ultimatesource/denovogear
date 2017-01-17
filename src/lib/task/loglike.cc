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

#include <boost/algorithm/string.hpp>

#include <dng/task/loglike.h>
#include <dng/probability.h>
#include <dng/fileio.h>
#include <dng/pileup.h>
#include <dng/read_group.h>
#include <dng/seq.h>
#include <dng/utility.h>
#include <dng/hts/bcf.h>
#include <dng/hts/extra.h>
#include <dng/vcfpileup.h>
#include <dng/mutation.h>
#include <dng/stats.h>
#include <dng/io/utility.h>
#include <dng/io/fasta.h>
#include <dng/depths.h>
#include <dng/multithread.h>

#include <htslib/faidx.h>
#include <htslib/khash.h>

#include "version.h"

using namespace std;
using namespace dng;
using namespace dng::task;

// Helper function for writing the vcf header information
void cout_add_header_text(task::LogLike::argument_type &arg) {
    using namespace std;
    string line{"##DeNovoGearCommandLine=<ID=dng-loglike,Version="
                PACKAGE_VERSION ","};
    line += utility::vcf_timestamp();
    line += ",CommandLineOptions=\"";

#define XM(lname, sname, desc, type, def) \
    line += utility::vcf_command_line_text(XS(lname),arg.XV(lname)) + ' ';
#   include <dng/task/loglike.xmh>
#undef XM
    for(auto && a : arg.input) {
        line += a + ' ';
    }

    line.pop_back();
    line += "\">";

    std::cout << line << "\n";
 }


// Sub-tasks
int process_bam(LogLike::argument_type &arg);
int process_ad(LogLike::argument_type &arg);

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
    FileCat mode = utility::input_category(*it, FileCat::Sequence|FileCat::Pileup, FileCat::Sequence);
    for(++it; it != arg.input.end(); ++it) {
        if(utility::input_category(*it, FileCat::Sequence|FileCat::Pileup, FileCat::Sequence) != mode) {
            throw runtime_error("Argument error: mixing sam/bam/cram and tad/ad input files is not supported.");
        }
    }
    // Execute sub tasks based on input type
    if(mode == FileCat::Pileup) {
        return process_ad(arg);
    }
    return process_bam(arg);
}

int process_bam(LogLike::argument_type &arg) {
    using namespace hts::bcf;

    // Parse pedigree from file
    io::Pedigree ped = io::parse_pedigree(arg.ped);

    // Open Reference
    io::Fasta reference{arg.fasta.c_str()};

    // Parse Nucleotide Frequencies
    std::array<double, 4> freqs = utility::parse_nuc_freqs(arg.nuc_freqs);

    // quality thresholds
    int min_qual = arg.min_basequal;

    // Print header to output
    cout_add_header_text(arg);

    // replace arg.region with the contents of a file if needed
    io::at_slurp(arg.region);

    // Open input files
    dng::ReadGroups rgs;
    vector<hts::bam::File> bamdata;
    for(auto && str : arg.input) {
        bamdata.emplace_back(str.c_str(), "r", arg.fasta.c_str(), arg.min_mapqual, arg.header.c_str());
        if(!bamdata.back().is_open()) {
            throw std::runtime_error("Error: Unable to open input file '" + str + "'.");
        }
        // add regions
        if(!arg.region.empty()) {
            bamdata.back().regions(regions::bam_parse_region(arg.region,bamdata.back()));
        }
        // Add each genotype/sample column
        rgs.ParseHeaderText(bamdata, arg.rgtag);
    }

    // Construct peeling algorithm from parameters and pedigree information
    InheritanceModel model = inheritance_model(arg.model);


    dng::RelationshipGraph graph;
    if(!graph.Construct(ped, rgs.GetLibraries(), model, arg.mu, arg.mu_somatic, arg.mu_library)) {
        throw std::runtime_error("Error: Unable to construct peeler for pedigree; "
                                 "possible non-zero-loop pedigree.");
    }
    // Select libraries in the input that are used in the pedigree
    rgs.SelectLibraries(graph.library_names());

    if(arg.gamma.size() < 2) {
        throw std::runtime_error("Error: Unable to construct genotype-likelihood model; "
                                 "Gamma needs to be specified at least twice to change model from default.");
    }

    for(auto && line : graph.BCFHeaderLines()) {
        std::cout << line << "\n";
    }

    LogProbability calculate (graph,
        { arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1] } );

    // Pileup data
    pileup::RawDepths read_depths(rgs.libraries().size());

    const int min_basequal = arg.min_basequal;
    auto filter_read = [min_basequal](
    dng::BamPileup::data_type::value_type::const_reference r) -> bool {
        return (r.is_missing
        || r.qual.first[r.pos] < min_basequal
        || seq::base_index(r.aln.seq_at(r.pos)) >= 4);
    };

    // Treat sequence_data and variant data separately
    dng::stats::ExactSum sum_data;
    dng::stats::ExactSum sum_scale;
    const bam_hdr_t *h = bamdata[0].header();
    dng::BamPileup mpileup{rgs.groups(), arg.min_qlen};
    mpileup(bamdata, [&](const dng::BamPileup::data_type & data, location_t loc) {
        // Calculate target position and fetch sequence name
        int contig = utility::location_to_contig(loc);
        int position = utility::location_to_position(loc);

        // Calculate reference base
        assert(0 <= contig && contig < h->n_targets);
        char ref_base = reference.FetchBase(h->target_name[contig],position);
        size_t ref_index = seq::char_index(ref_base);

        // reset all depth counters
        read_depths.assign(read_depths.size(), {});

        // pileup on read counts
        for(std::size_t u = 0; u < data.size(); ++u) {
            for(auto && r : data[u]) {
                if(filter_read(r)) {
                    continue;
                }
                std::size_t base = seq::base_index(r.aln.seq_at(r.pos));
                assert(read_depths[rgs.library_from_id(u)].counts[ base ] < 65535);

                read_depths[rgs.library_from_id(u)].counts[ base ] += 1;
            }
        }
        auto loglike = calculate(read_depths, ref_index);
        sum_data += loglike.log_data;
        sum_scale += loglike.log_scale;
    });
    // output results
    cout << "log_likelihood\tlog_hidden\tlog_observed\n";
    cout << setprecision(std::numeric_limits<double>::max_digits10)
        << sum_data.result()+sum_scale.result() << "\t"
        << sum_data.result() << "\t" << sum_scale.result() << "\n";

    return EXIT_SUCCESS;
}

int process_ad(LogLike::argument_type &arg) {
    using namespace std;
    // Parse pedigree from file
    io::Pedigree ped = io::parse_pedigree(arg.ped);

    // Open Reference
    io::Fasta reference{arg.fasta.c_str()};

    // Parse Nucleotide Frequencies
    array<double, 4> freqs = utility::parse_nuc_freqs(arg.nuc_freqs);

    // Print header to output
    cout_add_header_text(arg);

    // replace arg.region with the contents of a file if needed
    io::at_slurp(arg.region);

    // Open input files
    if(arg.input.size() != 1) {
        throw std::runtime_error("Argument Error: can only process one ad/tad file at a time.");
    }
    
    io::Ad input{arg.input[0], std::ios_base::in};
    if(!input) {
        throw std::runtime_error("Argument Error: unable to open input file '" + arg.input[0] + "'.");
    }
    if(input.ReadHeader() == 0) {
        throw std::runtime_error("Argument Error: unable to read header from '" + input.path() + "'.");
    }

    // Construct peeling algorithm from parameters and pedigree information
    InheritanceModel model = inheritance_model(arg.model);

    RelationshipGraph graph;
    if(!graph.Construct(ped, input.libraries(), model, arg.mu, arg.mu_somatic, arg.mu_library)) {
        throw runtime_error("Error: Unable to construct peeler for pedigree; "
                            "possible non-zero-loop pedigree.");
    }
    // Select libraries in the input that are used in the pedigree
    input.SelectLibraries(graph.library_names());

    if(arg.gamma.size() < 2) {
        throw runtime_error("Error: Unable to construct genotype-likelihood model; "
                            "Gamma needs to be specified at least twice to change model from default.");
    }

    for(auto && line : graph.BCFHeaderLines()) {
        cout << line << "\n";
    }

    // Construct function object to calculate log likelihoods
    LogProbability calculate (graph,
        { arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1] } );

    stats::ExactSum sum_data, sum_scale;

    // using a single thread to read data and process it
    if(arg.threads == 0) {
        pileup::AlleleDepths line;
        line.data().reserve(4*input.num_libraries());
        // read each line of data into line and process it
        while(input.Read(&line)) {
            auto loglike = calculate(line);
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
            for(auto && r : reads) {
                auto loglike = calculate(r);
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
    cout << "log_likelihood\tlog_hidden\tlog_observed\n";
    cout << setprecision(std::numeric_limits<double>::max_digits10)
        << sum_data.result()+sum_scale.result() << "\t"
        << sum_data.result() << "\t" << sum_scale.result() << "\n";

    return EXIT_SUCCESS;
}
