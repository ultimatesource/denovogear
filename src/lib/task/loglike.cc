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

bool is_missing_data(int i){
    return (i==hts::bcf::int32_missing);
}

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

int process_bam(LogLike::argument_type &arg) {
    // Parse pedigree from file
    Pedigree ped = io::parse_ped(arg.ped);

    // Open Reference
    if(arg.fasta.empty()){
    	throw std::invalid_argument("Path to reference file must be specified with --fasta when processing bam/sam/cram files.");
    }
    io::Fasta reference{arg.fasta.c_str()};

    // Parse Nucleotide Frequencies
    std::array<double, 4> freqs = utility::parse_nuc_freqs(arg.nuc_freqs);

    // Fetch the selected regions
    auto region_ext = io::at_slurp(arg.region); // replace arg.region with the contents of a file if needed

    // Open input files
    using BamPileup = dng::io::BamPileup;
    BamPileup mpileup{arg.min_qlen, arg.rgtag};
    for(auto && str : arg.input) {
        // TODO: We put all this logic here to simplify the BamPileup construction
        // TODO: However it might make sense to incorporate it into BamPileup::AddFile
        hts::bam::File input{str.c_str(), "r", arg.fasta.c_str(), arg.min_mapqual, arg.header.c_str()};
        if(!input.is_open()) {
            throw std::runtime_error("Unable to open input file '" + str + "'.");
        }
        // add regions
        // TODO: make this work with region_ext
        if(!arg.region.empty()) {
            if(arg.region.find(".bed") != std::string::npos) {
                input.regions(regions::bam_parse_bed(arg.region, input));
            } else {
                input.regions(regions::bam_parse_region(arg.region, input));
            }
        }
        mpileup.AddFile(std::move(input));
    }

    // Construct peeling algorithm from parameters and pedigree information
    RelationshipGraph relationship_graph;
    if (!relationship_graph.Construct(ped, mpileup.libraries(), inheritance_model(arg.model),
                                      arg.mu, arg.mu_somatic, arg.mu_library,
                                      arg.normalize_somatic_trees)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                 "possible non-zero-loop relationship_graph.");
    }

    // Select libraries in the input that are used in the pedigree
    mpileup.SelectLibraries(relationship_graph.library_names());

    // Print header to output
    cout_add_header_text(arg);

    for(auto && line : relationship_graph.BCFHeaderLines()) {
        std::cout << line << "\n";
    }

    LogProbability calculate (relationship_graph,
            { arg.theta, freqs, arg.ref_weight,
            {arg.lib_overdisp, arg.lib_error, arg.lib_bias} } );

    // Pileup data
    dng::pileup::RawDepths read_depths(mpileup.num_libraries());

    const size_t num_nodes = relationship_graph.num_nodes();
    const size_t library_start = relationship_graph.library_nodes().first;
    const int min_basequal = arg.min_basequal;
    auto filter_read = [min_basequal](
    BamPileup::data_type::value_type::const_reference r) -> bool {
        return (r.is_missing
        || r.base_qual() < min_basequal
        || seq::base_index(r.base()) >= 4);
    };

    // Treat sequence_data and variant data separately
    dng::stats::ExactSum sum_data;
    dng::stats::ExactSum sum_scale;
    auto h = mpileup.header();
    
    mpileup([&](const BamPileup::data_type & data, utility::location_t loc) {
        // Calculate target position and fetch sequence name
        int contig = utility::location_to_contig(loc);
        int position = utility::location_to_position(loc);

        // Calculate reference base
        assert(0 <= contig && contig < h->n_targets);
        char ref_base = reference.FetchBase(h->target_name[contig],position);
        size_t ref_index = seq::char_index(ref_base);

        // reset all depth counters
        boost::fill(read_depths, pileup::depth_t{});

        // pileup on read counts
        for(std::size_t u = 0; u < data.size(); ++u) {
            for(auto && r : data[u]) {
                if(filter_read(r)) {
                    continue;
                }
                std::size_t base = seq::base_index(r.base());
                assert(read_depths[u].counts[ base ] < 65535);

                read_depths[u].counts[ base ] += 1;
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
    Pedigree ped = io::parse_ped(arg.ped);

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
    if(!graph.Construct(ped, input.libraries(), model, arg.mu, arg.mu_somatic, arg.mu_library,
        arg.normalize_somatic_trees)) {
        throw runtime_error("Error: Unable to construct peeler for pedigree; "
                            "possible non-zero-loop pedigree.");
    }
    // Select libraries in the input that are used in the pedigree
    input.SelectLibraries(graph.library_names());

    for(auto && line : graph.BCFHeaderLines()) {
        cout << line << "\n";
    }

    // Construct function object to calculate log likelihoods
    LogProbability calculate (graph,
        { arg.theta, freqs, arg.ref_weight,
        {arg.lib_overdisp, arg.lib_error, arg.lib_bias} } );

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

int process_bcf(LogLike::argument_type &arg) {
    // Parse the pedigree file
    Pedigree ped = io::parse_ped(arg.ped);

    // Parse Nucleotide Frequencies
    auto freqs = utility::parse_nuc_freqs(arg.nuc_freqs);

    // Read input data
    if(arg.input.size() > 1) {
        throw std::runtime_error("dng call can only handle one variant file at a time.");
    }

    using dng::io::BcfPileup;
    BcfPileup mpileup;

    // Fetch the selected regions
    auto region_ext = io::at_slurp(arg.region); // replace arg.region with the contents of a file if needed
    if(!arg.region.empty()) {
        auto ranges = regions::parse_ranges(arg.region);
        if(ranges.second == false) {
            throw std::runtime_error("unable to parse the format of '--region' argument.");
        }
        mpileup.SetRegions(ranges.first);
    }

    if(mpileup.AddFile(arg.input[0].c_str()) == 0) {
        int errnum = mpileup.reader().handle()->errnum;
        throw std::runtime_error(bcf_sr_strerror(errnum));
    }

    // Construct peeling algorithm from parameters and pedigree information
    InheritanceModel model = inheritance_model(arg.model);

    dng::RelationshipGraph relationship_graph;
    if (!relationship_graph.Construct(ped, mpileup.libraries(), model,
                                      arg.mu, arg.mu_somatic, arg.mu_library,
                                      arg.normalize_somatic_trees)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                 "possible non-zero-loop relationship_graph.");
    }
    mpileup.SelectLibraries(relationship_graph.library_names());

    // Read header from first file
    const bcf_hdr_t *header = mpileup.reader().header(0); // TODO: fixthis
    const int num_libs = mpileup.num_libraries();

    // Print header to output
    cout_add_header_text(arg);

    for(auto && line : relationship_graph.BCFHeaderLines()) {
        std::cout << line << "\n";
    }

    LogProbability calculate (relationship_graph,
            { arg.theta, freqs, arg.ref_weight,
            {arg.lib_overdisp, arg.lib_error, arg.lib_bias} } );

    const size_t num_nodes = relationship_graph.num_nodes();
    const size_t library_start = relationship_graph.library_nodes().first;

    using pileup::AlleleDepths;

    AlleleDepths data;

    // allocate space for ad. bcf_get_format_int32 uses realloc internally
    int n_ad_capacity = num_libs*5;
    std::unique_ptr<int[],decltype(std::free)*> ad{
        reinterpret_cast<int*>(std::malloc(sizeof(int)*n_ad_capacity)), std::free};
    if(!ad) {
        throw std::bad_alloc{};
    }

    // Treat sequence_data and variant data separately
    dng::stats::ExactSum sum_data;
    dng::stats::ExactSum sum_scale;
    
    // run calculation based on the depths at each site.
    mpileup([&](const BcfPileup::data_type & rec) {
        data.location(rec->rid, rec->pos);

        // Identify site
        std::string allele_str;
        std::vector<char> allele_indexes;
        std::vector<int>  allele_list;
        const int num_alleles = rec->n_allele;
        int first_is_n = 0;
        for(int a = 0; a < num_alleles; ++a) {
            int n = AlleleDepths::MatchAlleles(rec->d.allele[a]);
            if(0 <= n && n <= 3) {
                allele_indexes.push_back(n);
                allele_list.push_back(a);
            } else if(a == 0 && n == 4) {
                first_is_n = 64;
            }
        }
        int color = AlleleDepths::MatchIndexes(allele_indexes);
        assert(color != -1);
        data.resize(color+first_is_n, rec->n_sample);

        // Read all the Allele Depths for every sample into an AD array
        int *pad = ad.get();
        int n_ad = bcf_get_format_int32(header, rec, "AD", &pad, &n_ad_capacity);
        if(n_ad == -4) {
            throw std::bad_alloc{};
        } else if(pad != ad.get()) {
            // update pointer
            ad.release();
            ad.reset(pad);
        }
        if(n_ad <= 0) {
            // AD tag is missing, so we do nothing at this time
            // TODO: support using calculated genotype likelihoods
            return;
        }

        assert(n_ad >= data.data_size());
        for(int i=0;i<num_libs;++i) {
            int offset = i*num_alleles;
            for(int a=0; a<allele_indexes.size();++a) {
                data(i,a) = ad[offset+allele_list[a]];
            }
        }

        // Replace missing data with 0
        boost::range::replace_if(data.data(), [](int n) {return n==hts::bcf::int32_missing;}, 0);

        auto loglike = calculate(data);
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
