//
// Created by steven on 1/11/16.
//
//#define BOOST_TEST_MODULE "VT_TEST"

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

#include <boost/range/iterator_range.hpp>
#include <boost/range/algorithm/replace.hpp>
#include <boost/range/algorithm/max_element.hpp>

#include <boost/algorithm/string.hpp>

#include <dng/task/call.h>
#include <dng/pedigree.h>
#include <dng/fileio.h>
#include <dng/pileup.h>
#include <dng/read_group.h>
#include <dng/likelihood.h>
#include <dng/seq.h>
#include <dng/utilities.h>
#include <dng/hts/bcf.h>
#include <dng/hts/extra.h>
#include <dng/vcfpileup.h>
#include <dng/mutation.h>
#include <dng/stats.h>

#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <dng/app.h>

#include "version.h"



#include "vcf_helper.h"
#include "find_mutation.h"
#include "find_mutation_getter.h"
#include "../../test/boost_test_helper.h"
//#include <boost/test/unit_test.hpp>

using namespace dng::task;
using namespace dng;




//TODO: Replace all test_** with BOOST_AUTO_TEST_CASE(Test_**)
//TODO: Replace main with some sort of setup() method in boost

void test_basic_parameterts(const FindMutationsGetter &find_mutation) {
    //Default values form call.xmh

//    BOOST_CHECK_EQUAL(find_mutation.getMin_prob_(), 0.01);
//    AssertTrue(find_mutation.getMin_prob_(), 0.01);
    assert(find_mutation.getMin_prob_() == 0.01);

    //getParams_()
    assert(find_mutation.getParams_().theta == 0.001);//test_param_1.theta);
    assert(find_mutation.getParams_().ref_weight == 1);//test_param_1.ref_weight);

    std::array<double, 4> expect_freqs{0.3, 0.2, 0.2, 0.3};
    auto freqs1 = find_mutation.getParams_().nuc_freq;
    for (int j = 0; j < 4; ++j) {
        assert(freqs1[j] == expect_freqs[j]);
    }
    std::vector<std::array<double, 4>> expect_gamma{{0.98, 0.0005, 0.0005, 1.04},
                                                    {0.02, 0.075,  0.005,  1.18}};
    auto gamma_a = find_mutation.getParams_().params_a;
    assert(gamma_a.pi == expect_gamma[0][0]);
    assert(gamma_a.phi == expect_gamma[0][1]);
    assert(gamma_a.epsilon == expect_gamma[0][2]);
    assert(gamma_a.omega == expect_gamma[0][3]);

    auto gamma_b = find_mutation.getParams_().params_b;
    assert(gamma_b.pi == expect_gamma[1][0]);
    assert(gamma_b.phi == expect_gamma[1][1]);
    assert(gamma_b.epsilon == expect_gamma[1][2]);
    assert(gamma_b.omega == expect_gamma[1][3]);

    auto event = find_mutation.getEvent_();
    assert(event.size() == 5);
    for (int k = 0; k < 5; ++k) {
        assert(event[k] == 0);
    }
}

void test_prior(const FindMutationsGetter &find_mutation) {

    std::vector<std::array<double, 10>> expected_prior {
//    calculated with the following parameters
//    double theta = 0.001;
//    std::array<double, 4> expect_freqs{0.3, 0.2, 0.2, 0.3};

/* R Code
prior0<- c(0,0,0,0)
freq<- c(0.3,0.2,0.2,0.3)
theta <- 0.001
tempComb<- combn(4,2)
allComb<-cbind(c(1,1), tempComb[,1:3], c(2,2), tempComb[,4:5], c(3,3), tempComb[,6], c(4,4))
result<- vector(mode="list", length=5)
for(i in 1:5){
    prior<- prior0
    if(i<5){
        prior[i]<- prior[i]+1
    }
    alpha <- freq*theta + prior
    alpha_sum <- sum(alpha)
    scale<- alpha_sum * (1+alpha_sum)
    result[[i]]<- apply(allComb, 2, function(x){
        if(x[1]==x[2]){
            return(alpha[x[1]]*(1+alpha[x[1]])/scale)
        }
        else{
            return(2*alpha[x[1]]*alpha[x[2]]/scale)
        }
    })
}

s<- sapply(result, function(x){
    paste(x, sep="", collapse=", ")
})
cat("{", paste(s, collapse="}, \n{"), "}\n" )
*/
            { 0.998951118846171, 0.000199760259730275, 0.000199760259730275, 0.000299640389595412, 9.98701448476561e-05, 3.99400699250774e-08, 5.99101048876161e-08, 9.98701448476561e-05, 5.99101048876161e-08, 0.000149820194797706},
            {0.000149820194797706, 0.000299610434542968, 5.99101048876161e-08, 8.98651573314242e-08, 0.998801318621409, 0.000199740289695312, 0.000299610434542968, 9.98701448476561e-05, 5.99101048876161e-08, 0.000149820194797706},
            {0.000149820194797706, 5.99101048876161e-08, 0.000299610434542968, 8.98651573314242e-08, 9.98701448476561e-05, 0.000199740289695312, 5.99101048876161e-08, 0.998801318621409, 0.000299610434542968, 0.000149820194797706},
            {0.000149820194797706, 5.99101048876161e-08, 5.99101048876161e-08, 0.000299640389595412, 9.98701448476561e-05, 3.99400699250774e-08, 0.000199760259730275, 9.98701448476561e-05, 0.000199760259730275, 0.998951118846171},
            {0.29979020979021, 0.00011988011988012, 0.00011988011988012, 0.00017982017982018, 0.19984015984016, 7.99200799200799e-05, 0.00011988011988012, 0.19984015984016, 0.00011988011988012, 0.29979020979021 }
    };

    auto pArray = find_mutation.getGenotype_prior_();
    for (int i = 0; i < 5; ++i) {
        auto array = pArray[i];
        for (int k = 0; k < 10; ++k) {
            assert(expected_prior[i][k] - array[k] < ABS_TEST_THRESHOLD);
        }
    }


}

void test_full_transition(const FindMutationsGetter &find_mutation) {

    std::array<double, 4> freqs{0.3, 0.2, 0.2, 0.3};
    double mu = 1e-8;
    auto dad = f81::matrix(3*mu, freqs);
    auto mom = f81::matrix(3*mu, freqs);

    auto exp_germline_full = meiosis_diploid_matrix(dad, mom);
    auto exp_germline_nomut = meiosis_diploid_matrix(dad, mom, 0);
    auto exp_germline_posmut = exp_germline_full - exp_germline_nomut;
    auto exp_germline_onemut = meiosis_diploid_matrix(dad, mom, 1);
    auto exp_germline_mean = meiosis_diploid_mean_matrix(dad, mom);



    auto orig = f81::matrix(2*mu, freqs);

    auto exp_somatic_full = mitosis_diploid_matrix(orig);
    auto exp_somatic_nomut = mitosis_diploid_matrix(orig, 0);
    auto exp_somatic_posmut = exp_somatic_full - exp_somatic_nomut;
    auto exp_somatic_onemut = mitosis_diploid_matrix(orig, 1);
    auto exp_somatic_mean = mitosis_diploid_mean_matrix(orig);

    auto full_matrices = find_mutation.getFull_transition_matrices_();
    auto nomut_matrices = find_mutation.getNomut_transition_matrices_();
    auto posmut_matrices = find_mutation.getPosmut_transition_matrices_();
    auto onemut_matrices = find_mutation.getOnemut_transition_matrices_();
    auto mean_matrices = find_mutation.getMean_matrices_();
    
    std::vector<TransitionMatrix> exp_full {{},{},exp_germline_full, exp_somatic_full, exp_somatic_full};
    std::vector<TransitionMatrix> exp_nomut {{},{},exp_germline_nomut, exp_somatic_nomut, exp_somatic_nomut};
    std::vector<TransitionMatrix> exp_posmut {{},{},exp_germline_posmut, exp_somatic_posmut, exp_somatic_posmut};
    std::vector<TransitionMatrix> exp_onemut {{},{},exp_germline_onemut, exp_somatic_onemut, exp_somatic_onemut};
    std::vector<TransitionMatrix> exp_mean {{},{},exp_germline_mean, exp_somatic_mean, exp_somatic_mean};

    for (int i = 0; i < 5; ++i) {
        AssertEigenMatrixNear(full_matrices[i],  exp_full[i]);
        AssertEigenMatrixNear(nomut_matrices[i], exp_nomut[i]);
        AssertEigenMatrixNear(posmut_matrices[i], exp_posmut[i]);
        AssertEigenMatrixNear(onemut_matrices[i], exp_onemut[i]);
        AssertEigenMatrixNear(mean_matrices[i], exp_mean[i]);
    }

}

void test_operator(const FindMutationsGetter &find_mutation) {
    FindMutations::stats_t stats;
    int ref_index = 0;
    std::vector<depth_t> read_depths(3);
    for (int j = 0; j < read_depths.size(); ++j) {
//        read_depths[j] = {0,0,0,0};
        read_depths[j].counts[0] = 10+j;
    }

//    find_mutation(read_depths, ref_index, &stats);


    std::vector<std::string> expect_gamma{"0.98, 0.0005, 0.0005, 1.04",
                                                    "0.02, 0.075,  0.005,  1.18"};
//    FindMutations::params_t(
//    { arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1] }  );
    dng::genotype::DirichletMultinomialMixture genotype_likelihood_ {
            dng::genotype::DirichletMultinomialMixture::params_t {expect_gamma[0]},
            dng::genotype::DirichletMultinomialMixture::params_t {expect_gamma[1]}
    };




    Pedigree pedigree = find_mutation.getPedigree_();

    auto family = pedigree.family_members_;
    dng::peel::workspace_t workspace;//  = pedigree.CreateWorkspace();
    
    workspace.Resize(5);
//    workspace.founder_nodes = std::make_pair(first_founder_, first_nonfounder_);
//    workspace.germline_nodes = std::make_pair(first_founder_, first_somatic_);
//    workspace.somatic_nodes = std::make_pair(first_somatic_, first_library_);
//    workspace.library_nodes = std::make_pair(first_library_, num_nodes_);
    workspace.founder_nodes = std::make_pair(0, 2);
    workspace.germline_nodes = std::make_pair(0, 2);
    workspace.somatic_nodes = std::make_pair(2, 2);
    workspace.library_nodes = std::make_pair(2, 5);

    double scale = 0.0, stemp;
    for(std::size_t u = 0; u < read_depths.size(); ++u) {
        std::tie(workspace.lower[2  + u], stemp) =
                genotype_likelihood_(read_depths[u], ref_index);
//        scale += stemp;
    }
    workspace.SetFounders(find_mutation.getGenotype_prior_()[ref_index]);//FIXME: lazy hack!

    auto numut_matrices = find_mutation.getNomut_transition_matrices_();
    workspace.CleanupFast();
    dng::peel::up(workspace, family[0], numut_matrices);
    dng::peel::to_father(workspace, family[1], numut_matrices);
    dng::peel::up(workspace, family[2], numut_matrices);

    double ret = log((workspace.lower[0] * workspace.upper[0]).sum());
    std::cout << "Expected: " << ret << std::endl;



    auto full_matrices = find_mutation.getFull_transition_matrices_();
    workspace.CleanupFast();
    dng::peel::up(workspace, family[0], full_matrices );
    dng::peel::to_father(workspace, family[1], full_matrices );
    dng::peel::up(workspace, family[2], full_matrices );

    double ret_full = log((workspace.lower[0] * workspace.upper[0]).sum());
    std::cout << "Expected full: " << ret_full << std::endl;


//    peeling_ops_;
//    // The modified, "faster" operations
//    std::vector<decltype(peel::op::NUM)> peeling_functions_ops_;
    // Peel pedigree one family at a time
//    for(std::size_t i = 0; i < peeling_functions_.size(); ++i) {
//    std::cout << "PeelForward: index: "<<  i << "\tfunction: " ;
//    //            std::cout << "\t" << *(peeling_functions_[i]) ;
//    //            std::cout << "\t" << (peeling_functions_[i]) ;
//    //
//    //            std::cout << "\t" << (peeling_ops_[i]) ;
//    std::cout << std::endl;
//    //            fprintf(std::cout, )
//    (*peeling_functions_[i])(workspace, family_members_[i], mat);
//    }
//    // Sum over roots
//    double ret = 0.0;
//    for(auto r : roots_) {
//    ret += log((workspace.lower[r] * workspace.upper[r]).sum());
//    }
//    bool FindMutations::operator()(const std::vector<depth_t> &depths,
//                                   int ref_index, stats_t *stats) {

}


// The main loop for dng-call application
// argument_type arg holds the processed command line arguments
int dng::task::Call::operator()(dng::task::Call::argument_type &arg) {
    using namespace std;
    using namespace hts::bcf;
    using dng::util::lphred;
    using dng::util::phred;

    std::cout << arg.mu << "\t" << arg.mu_somatic << "\t" << arg.mu_library << std::endl;


    // Parse pedigree from file
    dng::io::Pedigree ped;
    if (!arg.ped.empty()) {
        ifstream ped_file(arg.ped);
        if (!ped_file.is_open()) {
            throw std::runtime_error(
                    "unable to open pedigree file '" + arg.ped + "'.");
        }
        ped.Parse(istreambuf_range(ped_file));
    } else {
        throw std::runtime_error("pedigree file was not specified.");
    }

    // Open Reference
    unique_ptr<char[], void (*)(void *)> ref{nullptr, free};
    int ref_sz = 0, ref_target_id = -1;
    unique_ptr<faidx_t, void (*)(faidx_t *)> fai{nullptr, fai_destroy};
    if (!arg.fasta.empty()) {
        fai.reset(fai_load(arg.fasta.c_str()));
        if (!fai)
            throw std::runtime_error("unable to open faidx-indexed reference file '"
                                     + arg.fasta + "'.");
    }

    // Parse Nucleotide Frequencies
    std::array<double, 4> freqs;
    // TODO: read directly into freqs????  This will need a wrapper that provides an "insert" function.
    // TODO: include the size into the pattern, but this makes it harder to catch the second error.
    // TODO: turn all of this into a template function that returns array<double,4>?
    {
        auto f = util::parse_double_list(arg.nuc_freqs, ',', 4);
        if (!f.second) {
            throw std::runtime_error("Unable to parse nuc-freq option. "
                                             "It must be a comma separated list of floating-point numbers.");
        }
        if (f.first.size() != 4) {
            throw std::runtime_error("Wrong number of values passed to nuc-freq. "
                                             "Expected 4; found " + std::to_string(f.first.size()) + ".");
        }
        std::copy(f.first.begin(), f.first.end(), &freqs[0]);
    }

    // quality thresholds
    int min_qual = arg.min_basequal;
    double min_prob = arg.min_prob;

    // Open input files
    dng::ReadGroups rgs;
    vector<hts::File> indata;
    vector<hts::bam::File> bamdata;
    vector<hts::bcf::File> bcfdata;
    for (auto &&str : arg.input) {
        indata.emplace_back(str.c_str(), "r");
        if (indata.back().is_open()) {
            continue;
        }
        throw std::runtime_error("unable to open input file '" + str + "'.");
    }

    // Check to see if all inputs are of the same type
    const htsFormatCategory cat = indata[0].format().category;
    for (auto &&f : indata) {
        if (f.format().category == cat) {
            continue;
        }
        throw std::runtime_error("mixing sequence data and variant data as input is not supported.");
    }

    // Begin writing VCF header
    auto out_file = vcf_get_output_mode(arg);
    hts::bcf::File vcfout(out_file.first.c_str(), out_file.second.c_str());
    vcf_add_header_text(vcfout, arg);

    if (cat == sequence_data) {
        // Wrap input in hts::bam::File
        for (auto &&f : indata) {
            bamdata.emplace_back(std::move(f), arg.region.c_str(), arg.fasta.c_str(),
                                 arg.min_mapqual, arg.header.c_str());
        }

        // Read header from first file
        const bam_hdr_t *h = bamdata[0].header();

        // Add contigs to header
        for (auto &&contig : parse_contigs(h)) {
            vcfout.AddContig(contig.first.c_str(), contig.second);
        }

        // Add each genotype/sample column
        rgs.ParseHeaderText(bamdata, arg.rgtag);
    } else if (cat == variant_data) {
        cout << "CALL VCF DATA\n" << endl;
        bcfdata.emplace_back(std::move(indata[0]));
        // Read header from first file
        const bcf_hdr_t *h = bcfdata[0].header();



        // Add contigs to header
        for (auto &&contig : extract_contigs(h)) {
            vcfout.AddHeaderMetadata(contig.c_str());
//            cout << contig.c_str() << endl;
        }
        // Add each genotype/sample column
        rgs.ParseSamples(bcfdata[0]);
    } else {
        throw runtime_error("unsupported file category.");
    }

    // Construct peeling algorithm from parameters and pedigree information
    dng::Pedigree pedigree;
    cout << "MU: " << arg.mu << "\t" << arg.mu_somatic << "\t" << arg.mu_library << std::endl;
    if (!pedigree.Construct(ped, rgs, arg.mu, arg.mu_somatic, arg.mu_library)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                         "possible non-zero-loop pedigree.");
    }
    if (arg.gamma.size() < 2) {
        throw std::runtime_error("Unable to construct genotype-likelihood model; "
                                         "Gamma needs to be specified at least twice to change model from default.");
    }

    for (auto &&line : pedigree.BCFHeaderLines()) {
        vcfout.AddHeaderMetadata(line.c_str());
//        cout << line.c_str() << endl;
    }
    pedigree.PrintTable(std::cout);
    pedigree.PrintMachine(std::cout);
    for(auto && line : pedigree.BCFHeaderLines()) {
        cout << line.c_str() << endl;
    }
//    pedigree.PrintStates(std::cout);

    std::cout << "Founder, nonfounder, somatic, library, num_nodes" << std::endl;
    cout << pedigree.first_founder_ << endl;    // start of founder germline
    cout << pedigree.first_nonfounder_ << endl; // start of non-founder germline
    cout << pedigree.first_somatic_ << endl;    // start of somatic nodes
    cout << pedigree.first_library_ << endl;    // start of libraries
    cout << pedigree.num_nodes_ << endl;        // total number of nodes

    cout << "Roots\t"<< pedigree.roots_.size() << endl;

    // Pedigree Structure
    cout << "Labels: " << pedigree.labels_.size() << endl;
    cout << pedigree.transitions_.size() << endl;

    // The original, simplified peeling operations
    cout << pedigree.peeling_ops_.size() << endl;
    for (auto a : pedigree.peeling_ops_) {
        std::cout << "op: "<< a << std::endl;
    }

    // The modified, "faster" operations
    cout << pedigree.peeling_functions_ops_.size() << endl;
    for (auto a : pedigree.peeling_functions_ops_) {
        std::cout << "func op: "<< a << std::endl;
    }

    // Array of functions that will be called to perform the peeling
    cout << pedigree.peeling_functions_.size() << endl;
    cout << pedigree.peeling_reverse_functions_.size() << endl;
    for (auto a : pedigree.peeling_functions_) {
        std::cout << "functions: "<< a << std::endl;
    }



    // The arguments to a peeling operation
    cout << pedigree.family_members_.size() << endl;


//    FindMutations calculate{min_prob, pedigree,
//                            {arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1]}};

    FindMutations::params_t test_param_1{arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1]};
    // Pileup data
    FindMutationsGetter find_mutation{min_prob, pedigree, test_param_1};
    test_basic_parameterts(find_mutation);
    test_prior(find_mutation);

    test_full_transition(find_mutation);

    test_operator(find_mutation);

    std::cout << "ready to exit: " << EXIT_SUCCESS << std::endl;
    return EXIT_SUCCESS;

}





int main(int argc, char *argv[]) {


//BOOST_AUTO_TEST_CASE(Test_FM){


    try {
//        return CallApp(argc, argv)();
//        dng::CommandLineApp<dng::task::Call> a (argc, argv) ;

//        char *argv[] = {"test", "-p", "arg2", NULL};
//        int argc = sizeof(argv) / sizeof(char*) - 1;
        int argc=4;
        char *argv[argc+1];
        argv[0] = (char*) "test";
        argv[1] = (char*) "-p";
        argv[2] = (char*) "testdata/sample_5_3/ceu.ped"; //"pedFile";
        argv[3] = (char*) "testdata/sample_5_3/test1.vcf"; //test1.bam
        dng::CommandLineApp<dng::task::Call> a (argc, argv) ;
        a();

    } catch(std::exception &e) {
        std::cerr << e.what() << std::endl;
    }
    return EXIT_FAILURE;

}

//
// get unknow error message
//unknown location(0): fatal error in "Test_FM": signal: illegal operand; address of failing instruction: 0x0061ef82
//
//*** 1 failure detected in test suite "VT_TEST"

//BOOST_AUTO_TEST_CASE(Test_FM){
//    try {
////        return CallApp(argc, argv)();
////        dng::CommandLineApp<dng::task::Call> a (argc, argv) ;
//
////        char *argv[] = {"test", "-p", "arg2", NULL};
////        int argc = sizeof(argv) / sizeof(char*) - 1;
//        int argc=4;
//        char *argv[argc+1];
//        char const *p = "abc";
//        argv[0] = (char*) "test";
//        argv[1] = (char*) "-p";
//        argv[2] = "testdata/sample_5_3/ceu.ped"; //"pedFile";
//        argv[3] = "testdata/sample_5_3/test1.vcf"; //test1.bam
//        dng::CommandLineApp<dng::task::Call> a (argc, argv) ;
//        std::cout << "DIV" << std::endl;
//        BOOST_CHECK_EQUAL(0.01, 0.011);
//        a();
//    } catch(std::exception &e) {
//        std::cerr << e.what() << std::endl;
//    }
////    return EXIT_FAILURE;
//
//}
//
