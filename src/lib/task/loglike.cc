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
#include <dng/relationship_graph.h>
#include <dng/fileio.h>
#include <dng/pileup.h>
#include <dng/read_group.h>
#include <dng/likelihood.h>
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

class LogProbability {
public:
    struct params_t {
        double theta;
        std::array<double, 4> nuc_freq;
        double ref_weight;

        dng::genotype::DirichletMultinomialMixture::params_t params_a;
        dng::genotype::DirichletMultinomialMixture::params_t params_b;
    };

    struct stats_t {
        double log_data;
        double log_scale;        
    };

    LogProbability(RelationshipGraph pedigree, params_t params);

    stats_t operator()(const std::vector<depth_t> &depths, int ref_index);
    stats_t operator()(const pileup::AlleleDepths &depths, const std::vector<size_t>& indexes);

protected:
    dng::RelationshipGraph pedigree_;
    params_t params_;
    dng::peel::workspace_t work_; // must be declared after pedigree_ (see constructor)

    dng::TransitionVector full_transition_matrices_;
    dng::TransitionVector transition_matrices_[dng::pileup::AlleleDepths::type_info_table_length];
    double prob_monomorphic_[4];

    // Model genotype likelihoods as a mixture of two dirichlet multinomials
    // TODO: control these with parameters
    dng::genotype::DirichletMultinomialMixture genotype_likelihood_;

    dng::GenotypeArray genotype_prior_[5]; // Holds P(G | theta)
};

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
            throw std::runtime_error("unable to open input file '" + str + "'.");
        }
        // add regions
        if(!arg.region.empty()) {
            bamdata.back().regions(regions::bam_parse_region(arg.region,bamdata.back()));
        }
        // Add each genotype/sample column
        rgs.ParseHeaderText(bamdata, arg.rgtag);
    }

    // Construct peeling algorithm from parameters and pedigree information
    dng::RelationshipGraph pedigree;
    if(!pedigree.Construct(ped, rgs, arg.mu, arg.mu_somatic, arg.mu_library)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                 "possible non-zero-loop pedigree.");
    }
    if(arg.gamma.size() < 2) {
        throw std::runtime_error("Unable to construct genotype-likelihood model; "
                                 "Gamma needs to be specified at least twice to change model from default.");
    }

    for(auto && line : pedigree.BCFHeaderLines()) {
        std:cout << line << "\n";
    }

    LogProbability calculate (pedigree,
        { arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1] } );

    // Pileup data
    std::vector<depth_t> read_depths(rgs.libraries().size());

    const size_t num_nodes = pedigree.num_nodes();
    const size_t library_start = pedigree.library_nodes().first;
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

    // quality thresholds
    int min_qual = arg.min_basequal;

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
    // Construct ReadGroups
    ReadGroups rgs;
    rgs.ParseLibraries(input.libraries());

    // Construct peeling algorithm from parameters and pedigree information
    RelationshipGraph pedigree;
    if(!pedigree.Construct(ped, rgs, arg.mu, arg.mu_somatic, arg.mu_library)) {
        throw runtime_error("Unable to construct peeler for pedigree; "
                            "possible non-zero-loop pedigree.");
    }
    // Since pedigree may have removed libraries, map libraries to positions
    vector<size_t> library_to_index;
    library_to_index.resize(rgs.libraries().size());
    for(size_t u=0; u < input.libraries().size(); ++u) {
        size_t pos = rg::index(rgs.libraries(), input.library(u).name);
        if(pos == -1) {
            continue;
        }
        library_to_index[pos] = u;
    }

    if(arg.gamma.size() < 2) {
        throw runtime_error("Unable to construct genotype-likelihood model; "
                            "Gamma needs to be specified at least twice to change model from default.");
    }

    for(auto && line : pedigree.BCFHeaderLines()) {
        cout << line << "\n";
    }

    // Construct function object to calculate log likelihoods
    LogProbability calculate (pedigree,
        { arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1] } );

    stats::ExactSum sum_data, sum_scale;

    // using a single thread to read data and process it
    if(arg.threads == 0) {
        pileup::AlleleDepths line;
        line.data().reserve(4*input.libraries().size());
        // read each line of data into line and process it
        while(input.Read(&line)) {
            auto loglike = calculate(line,library_to_index);
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
                auto loglike = calculate(r,library_to_index);
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

// pedigree_ will be initialized before work_, so we will reference it.
LogProbability::LogProbability(RelationshipGraph pedigree, params_t params) :
    pedigree_{std::move(pedigree)},
    params_(std::move(params)), genotype_likelihood_{params.params_a, params.params_b},
    work_{pedigree_.CreateWorkspace()} {

    using namespace dng;

    // Use a parent-independent mutation model, which produces a
    // beta-binomial
    genotype_prior_[0] = population_prior(params_.theta, params_.nuc_freq, {params_.ref_weight, 0, 0, 0});
    genotype_prior_[1] = population_prior(params_.theta, params_.nuc_freq, {0, params_.ref_weight, 0, 0});
    genotype_prior_[2] = population_prior(params_.theta, params_.nuc_freq, {0, 0, params_.ref_weight, 0});
    genotype_prior_[3] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, params_.ref_weight});
    genotype_prior_[4] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, 0});

    // Calculate mutation matrices
    full_transition_matrices_.assign(work_.num_nodes, {});
 
    for(size_t child = 0; child < work_.num_nodes; ++child) {
        auto trans = pedigree_.transitions()[child];
        if(trans.type == RelationshipGraph::TransitionType::Germline) {
            auto dad = f81::matrix(trans.length1, params_.nuc_freq);
            auto mom = f81::matrix(trans.length2, params_.nuc_freq);

            full_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom);
        } else if(trans.type == RelationshipGraph::TransitionType::Somatic ||
                  trans.type == RelationshipGraph::TransitionType::Library) {
            auto orig = f81::matrix(trans.length1, params_.nuc_freq);

            full_transition_matrices_[child] = mitosis_diploid_matrix(orig);
        } else {
            full_transition_matrices_[child] = {};
        }
    }

    // Extract relevant subsets of matrices
    for(size_t color = 0; color < dng::pileup::AlleleDepths::type_info_table_length; ++color) {
        int width = dng::pileup::AlleleDepths::type_info_gt_table[color].width;
        // Resize our subsets to the right width
        transition_matrices_[color].resize(full_transition_matrices_.size());
        // enumerate over all children
        for(size_t child = 0; child < work_.num_nodes; ++child) {
            auto trans = pedigree_.transitions()[child];
            if(trans.type == RelationshipGraph::TransitionType::Germline) {
                // resize transition matrix to w*w,w
                transition_matrices_[color][child].resize(width*width,width);
                // Assume column major order which is the default
                for(int a = 0; a < width; ++a) {
                    int ga = dng::pileup::AlleleDepths::type_info_gt_table[color].indexes[a];
                    for(int b = 0; b < width; ++b) {
                        int gb = dng::pileup::AlleleDepths::type_info_gt_table[color].indexes[b];
                        // copy correct value from the full matrix to the subset matrix
                        int x = a*width+b;
                        int gx = ga*10+gb;
                        for(int y = 0; y < width; ++y) {
                            int gy = dng::pileup::AlleleDepths::type_info_gt_table[color].indexes[y];
                            // copy correct value from the full matrix to the subset matrix
                            transition_matrices_[color][child](x,y) = full_transition_matrices_[child](gx,gy); 
                        }
                    }
                }
            } else if(trans.type == RelationshipGraph::TransitionType::Somatic ||
                      trans.type == RelationshipGraph::TransitionType::Library) {
                transition_matrices_[color][child].resize(width,width);
                // Assume column major order which is the default
                for(int x = 0; x < width; ++x) {
                    int gx = dng::pileup::AlleleDepths::type_info_gt_table[color].indexes[x];
                    for(int y = 0; y < width; ++y) {
                        int gy = dng::pileup::AlleleDepths::type_info_gt_table[color].indexes[y];
                        // copy correct value from the full matrix to the subset matrix
                        transition_matrices_[color][child](x,y) = full_transition_matrices_[child](gx,gy); 
                    }
                }
            } else {
                transition_matrices_[color][child] = {};
            }
        }
    }

    // Precalculate monomorphic histories (first 4 colors)
    size_t num_libraries = work_.library_nodes.second - work_.library_nodes.first;
    pileup::AlleleDepths depths{0,0,num_libraries,pileup::AlleleDepths::data_t(num_libraries, 0)};
    std::vector<size_t> indexes(num_libraries,0);
    std::iota(indexes.begin(),indexes.end(),0);

    for(int color=0; color<4;++color) {
        // setup monomorphic prior
        depths.color(color);
        GenotypeArray prior(1);
        prior(0) = genotype_prior_[color](pileup::AlleleDepths::type_info_gt_table[color].indexes[0]);
        work_.SetFounders(prior);
        double scale = genotype_likelihood_(depths, indexes, work_.lower.begin()+work_.library_nodes.first);
        double logdata = pedigree_.PeelForwards(work_, transition_matrices_[color]);
        prob_monomorphic_[color] = exp(logdata);
    }
}

// returns 'log10 P(Data ; model)-log10 scale' and log10 scaling.
LogProbability::stats_t LogProbability::operator()(const std::vector<depth_t> &depths,
                               int ref_index) {
    // calculate genotype likelihoods and store in the lower library vector
    double scale = 0.0, stemp;
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::tie(work_.lower[work_.library_nodes.first + u], stemp) =
            genotype_likelihood_(depths[u], ref_index);
        scale += stemp;
    }

    // Set the prior probability of the founders given the reference
    work_.SetFounders(genotype_prior_[ref_index]);

    // Calculate log P(Data ; model)
    double logdata = pedigree_.PeelForwards(work_, full_transition_matrices_);

    return {logdata/M_LN10, scale/M_LN10};
}

// Calculate the probability of a depths object considering only indexes
// TODO: make indexes a property of the pedigree
LogProbability::stats_t LogProbability::operator()(const pileup::AlleleDepths &depths, const std::vector<size_t>& indexes) {
    int ref_index = depths.type_info().reference;
    int color = depths.color();
    size_t gt_width = depths.type_gt_info().width;

    double scale, logdata;
    // For monomorphic sites we have pre-calculated the peeling part
    if(color < 4) {
        assert(gt_width == 1 && ref_index == color);
        // Calculate genotype likelihoods
        scale = genotype_likelihood_(depths, indexes, work_.lower.begin()+work_.library_nodes.first);

        // Multiply our pre-calculated peeling results with the genotype likelihoods
        logdata = prob_monomorphic_[color];
        for(auto it = work_.lower.begin()+work_.library_nodes.first;
            it != work_.lower.begin()+work_.library_nodes.second; ++it) {
            logdata *= (*it)(0);
        }
        // convert to a log-likelihood
        logdata = log(logdata);
    } else {
        // Set the prior probability of the founders given the reference
        GenotypeArray prior(gt_width);
        for(int i=0;i<gt_width;++i) {
            prior(i) = genotype_prior_[ref_index](depths.type_gt_info().indexes[i]);
        }
        work_.SetFounders(prior);
    
        // Calculate genotype likelihoods
        scale = genotype_likelihood_(depths, indexes, work_.lower.begin()+work_.library_nodes.first);

         // Calculate log P(Data ; model)
        logdata = pedigree_.PeelForwards(work_, transition_matrices_[color]);
    }
    return {logdata/M_LN10, scale/M_LN10};
}
