//
// Created by steven on 1/15/16.
//

//#define BOOST_TEST_MODULE "dng::util::seq.h"

//#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>

#include "dng/seq.h"

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
//#include <boost/test/unit_test.hpp>
#include "find_mutation_getter.h"
#include "assert_helper.h"



using namespace dng::task;
using namespace dng;

using namespace std;

void test1(const FindMutationsGetter &find_mutation){

    FindMutations::stats_t stats;
    int ref_index = 0;
    std::vector<depth_t> read_depths(3);
    for (int j = 0; j < read_depths.size(); ++j) {
//        read_depths[j] = {0,0,0,0};
        read_depths[j].counts[0] = 10+j;
    }




    std::vector<std::string> expect_gamma{"0.98, 0.0005, 0.0005, 1.04",
                                          "0.02, 0.075,  0.005,  1.18"};
    dng::genotype::DirichletMultinomialMixture genotype_likelihood_ {
            dng::genotype::DirichletMultinomialMixture::params_t {expect_gamma[0]},
            dng::genotype::DirichletMultinomialMixture::params_t {expect_gamma[1]}
    };


    Pedigree pedigree = find_mutation.getPedigree_();

    auto families = pedigree.family_members_;


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
    dng::peel::up(workspace, families[0], numut_matrices);
    dng::peel::to_father(workspace, families[1], numut_matrices);
    dng::peel::up(workspace, families[2], numut_matrices);

    double ret = log((workspace.lower[0] * workspace.upper[0]).sum());
    std::cout << "Expected: " << ret << std::endl;


}



void TestSumOverChild() {



    int child_offset = 2;
    for (int t = 0; t < 100; ++t) {
        peel::family_members_t family1 {0,1}; //0 parent, 1 child
        dng::peel::workspace_t workspace;

        int num_child = std::rand()%10+1;

        int total_family_size = num_child + child_offset;
        vector<TransitionMatrix> m;
        vector<GenotypeArray> g;
        m.resize(total_family_size);
        g.resize(total_family_size);
        for (int k = 2; k < 2+num_child; ++k) {
                                                                                                                                                                                                                                                                                                                    family1.push_back(k);
            m[k] = TransitionMatrix::Random(10, 10);
            g[k] = GenotypeArray::Random();
        }

        GenotypeArray expected = GenotypeArray::Ones();
        for (int k = 2; k < total_family_size; ++k) {
            GenotypeArray temp_array = GenotypeArray::Zero();
            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 10; ++j) {
                    temp_array[i] += m[k](i, j) * g[k][j];
                }
            }
            expected *= temp_array;
        }

        workspace.Resize(total_family_size);
        dng::TransitionVector full_matrix;
        full_matrix.resize(total_family_size);
        full_matrix[0] = {};
        full_matrix[1] = {};
        for (int k = 2; k < total_family_size; ++k) {
            full_matrix[k] = m[k];
            workspace.lower[k] = g[k];
        }


        GenotypeArray result = dng::peel::sum_over_child(workspace, family1, full_matrix);
        for (int i = 0; i < 10; ++i) {
            AssertNear(expected[i], result[i]);
        }
    }



}


void TestUpCore(){

    peel::family_members_t family1 {0,1}; //0 parent, 1 child
    dng::peel::workspace_t workspace;

    for (int t = 0; t < 10; ++t) {
        const TransitionMatrix m = TransitionMatrix::Random(10, 10);
        const GenotypeArray g = GenotypeArray::Random();
        GenotypeArray expected = GenotypeArray::Zero();
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                expected[i] += m(i, j) * g[j];
            }
        }

        workspace.Resize(2);
        dng::TransitionVector full_matrix{2};
        full_matrix[0] = {};
        full_matrix[1] = m;
        workspace.lower[1] = g;

        GenotypeArray result = dng::peel::up_core(workspace, family1, full_matrix);
        for (int i = 0; i < 10; ++i) {
            AssertNear(expected[i], result[i]);
        }
    }
}

void TestUp(){

    peel::family_members_t family1 {0,1}; //0 parent, 1 child
    dng::peel::workspace_t workspace;

    for (int t = 0; t < 10; ++t) {
        const TransitionMatrix m = TransitionMatrix::Random(10, 10);
        const GenotypeArray g0 = GenotypeArray::Random();
        const GenotypeArray g1 = GenotypeArray::Random();
        GenotypeArray expected = GenotypeArray::Zero();
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                expected[i] += m(i, j) * g1[j];
            }
            expected[i] *= g0[i];
        }

        workspace.Resize(2);
        dng::TransitionVector full_matrix{2};
        full_matrix[0] = {};
        full_matrix[1] = m;
        workspace.lower[0] = g0;
        workspace.lower[1] = g1;

        dng::peel::up(workspace, family1, full_matrix);
        GenotypeArray result = workspace.lower[0];
        for (int i = 0; i < 10; ++i) {
            AssertNear(expected[i], result[i]);
        }
    }
}



void TestToFather(){



//// Family Order: Father, Mother, Child1, Child2, ...
//    void dng::peel::to_father(workspace_t &work, const family_members_t &family,
//                              const TransitionVector &mat) {
//        std::cout << " to father: " << family[0] << "\t" << family[1] << "\t" << family[2] << std::endl;
//        assert(family.size() >= 3);
//        auto dad = family[0];
//        auto mom = family[1];
//        // Sum over children
//        work.paired_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
//        for(std::size_t i = 3; i < family.size(); ++i) {
//            work.paired_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
//        }
//        // Include Mom
//        work.paired_buffer.resize(10, 10);
//        work.lower[dad] *= (work.paired_buffer.matrix() * (work.upper[mom] *
//                                                           work.lower[mom]).matrix()).array();
//        work.paired_buffer.resize(100, 1); //Might not need this, from the website: Assignment is the action of copying a matrix into another, using operator=. Eigen resizes the matrix on the left-hand side automatically so that it matches the size of the matrix on the right-hand size. For example:
//    }

    peel::family_members_t family1 {0, 1, 2}; //0 dad, 1 , 2 child
    dng::peel::workspace_t workspace;

    for (int t = 0; t < 10; ++t) {
        const TransitionMatrix m = TransitionMatrix::Random(100, 10);
        const GenotypeArray g0 = GenotypeArray::Random();
        const GenotypeArray g1 = GenotypeArray::Random();
        const GenotypeArray g2 = GenotypeArray::Random();
        const GenotypeArray u0 = GenotypeArray::Random();
        const GenotypeArray u1 = GenotypeArray::Random();

        PairedGenotypeArray all_child = PairedGenotypeArray::Zero(100, 1);
        for (int i = 0; i < 100; ++i) {
            for (int j = 0; j < 10; ++j) {
                all_child(i, 0) += m(i, j) * g2[j];
            }
        }
        all_child.resize(10, 10);

        GenotypeArray ga_mum;
        for (int j = 0; j < 10; ++j) {
            ga_mum[j] = u1[j] * g1[j];
        }

        GenotypeArray expected = GenotypeArray::Zero();
        GenotypeArray expected_fast = GenotypeArray::Zero();
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                expected[i] += all_child(i, j) * ga_mum[j];
            }
            expected_fast[i] = expected[i];
            expected[i] *= g0[i];
        }
//        std::cout << expected.sum() << std::endl;

        workspace.Resize(3);
        dng::TransitionVector full_matrix{3};
        full_matrix[0] = {};
        full_matrix[1] = {};
        full_matrix[2] = m;
        workspace.lower[0] = g0;
        workspace.lower[1] = g1;
        workspace.lower[2] = g2;
        workspace.upper[0] = u0;
        workspace.upper[1] = u1;
        workspace.upper[2] = {};

        dng::peel::to_father(workspace, family1, full_matrix);
        GenotypeArray result = workspace.lower[0];
        for (int i = 0; i < 10; ++i) {
            AssertNear(expected[i], result[i]);
        }

        dng::peel::to_father_fast(workspace, family1, full_matrix);
        GenotypeArray result_fast = workspace.lower[0];
        for (int i = 0; i < 10; ++i) {
            AssertNear(expected_fast[i], result_fast[i]);
        }


    }
}



void TestToMother(){



//// Family Order: Father, Mother, Child1, Child2, ...
//    void dng::peel::to_father(workspace_t &work, const family_members_t &family,
//                              const TransitionVector &mat) {
//        std::cout << " to father: " << family[0] << "\t" << family[1] << "\t" << family[2] << std::endl;
//        assert(family.size() >= 3);
//        auto dad = family[0];
//        auto mom = family[1];
//        // Sum over children
//        work.paired_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
//        for(std::size_t i = 3; i < family.size(); ++i) {
//            work.paired_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
//        }
//        // Include Mom
//        work.paired_buffer.resize(10, 10);
//        work.lower[dad] *= (work.paired_buffer.matrix() * (work.upper[mom] *
//                                                           work.lower[mom]).matrix()).array();
//        work.paired_buffer.resize(100, 1); //Might not need this, from the website: Assignment is the action of copying a matrix into another, using operator=. Eigen resizes the matrix on the left-hand side automatically so that it matches the size of the matrix on the right-hand size. For example:
//    }

    peel::family_members_t family1 {0, 1, 2}; //0 dad, 1 , 2 child
    dng::peel::workspace_t workspace;

    for (int t = 0; t < 10; ++t) {
        const TransitionMatrix m = TransitionMatrix::Random(100, 10);
        const GenotypeArray g0 = GenotypeArray::Random();
        const GenotypeArray g1 = GenotypeArray::Random();
        const GenotypeArray g2 = GenotypeArray::Random();
        const GenotypeArray u0 = GenotypeArray::Random();
        const GenotypeArray u1 = GenotypeArray::Random();

        PairedGenotypeArray all_child = PairedGenotypeArray::Zero(100, 1);
        for (int i = 0; i < 100; ++i) {
            for (int j = 0; j < 10; ++j) {
                all_child(i, 0) += m(i, j) * g2[j];
            }
        }
        all_child.resize(10, 10);

        GenotypeArray ga_dad;
        for (int j = 0; j < 10; ++j) {
            ga_dad[j] = u0[j] * g0[j];
        }

        GenotypeArray expected = GenotypeArray::Zero();
        GenotypeArray expected_fast = GenotypeArray::Zero();
        auto temp_child = all_child.transpose();
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                expected[i] += temp_child(i, j) * ga_dad[j];
            }
            expected_fast[i] = expected[i];
            expected[i] *= g1[i];
        }
//        std::cout << expected.sum() << std::endl;

        workspace.Resize(3);
        dng::TransitionVector full_matrix{3};
        full_matrix[0] = {};
        full_matrix[1] = {};
        full_matrix[2] = m;
        workspace.lower[0] = g0;
        workspace.lower[1] = g1;
        workspace.lower[2] = g2;
        workspace.upper[0] = u0;
        workspace.upper[1] = u1;
        workspace.upper[2] = {};

        dng::peel::to_mother(workspace, family1, full_matrix);
        GenotypeArray result = workspace.lower[1];
        for (int i = 0; i < 10; ++i) {
            AssertNear(expected[i], result[i]);
        }

        dng::peel::to_mother_fast(workspace, family1, full_matrix);
        GenotypeArray result_fast = workspace.lower[1];
        for (int i = 0; i < 10; ++i) {
            AssertNear(expected_fast[i], result_fast[i]);
        }


    }
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


    // The arguments to a peeling operation
    cout << pedigree.family_members_.size() << endl;


//    FindMutations calculate{min_prob, pedigree,
//                            {arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1]}};

    FindMutations::params_t test_param_1{arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1]};
    // Pileup data
    FindMutationsGetter find_mutation{min_prob, pedigree, test_param_1};

    test1(find_mutation);

    std::cout << "ready to exit: " << EXIT_SUCCESS << std::endl;
    return EXIT_SUCCESS;

}





int main(int argc, char *argv[]) {
//BOOST_AUTO_TEST_CASE(Test_FM){
//    TestUpCore();
//    TestUp();
    TestToFather();
    TestToMother();
//    TestSumOverChild();
//    try {
////        return CallApp(argc, argv)();
////        dng::CommandLineApp<dng::task::Call> a (argc, argv) ;
//
////        char *argv[] = {"test", "-p", "arg2", NULL};
////        int argc = sizeof(argv) / sizeof(char*) - 1;
//        int argc=4;
//        char *argv[argc+1];
//        argv[0] = (char*) "test";
//        argv[1] = (char*) "-p";
//        argv[2] = (char*) "testdata/sample_5_3/ceu.ped"; //"pedFile";
//        argv[3] = (char*) "testdata/sample_5_3/test1.vcf"; //test1.bam
////        dng::CommandLineApp<dng::task::Call> a (argc, argv) ;
////        a();
//
//    } catch(std::exception &e) {
//        std::cerr << e.what() << std::endl;
//    }
    std::cout << "Tests done"<< std::endl;
    return EXIT_FAILURE;

}

//BOOST_AUTO_TEST_CASE(Test_char_index)
//{
//    BOOST_CHECK(dng::seq::char_index('N') == 4);
//    BOOST_CHECK(dng::seq::char_index('A') == 0);
//    BOOST_CHECK(dng::seq::char_index('C') == 1);
//    BOOST_CHECK(dng::seq::char_index('G') == 2);
//    BOOST_CHECK(dng::seq::char_index('T') == 3);
//}
//
//BOOST_AUTO_TEST_CASE(Test_indexed_char)
//{
//    BOOST_CHECK(dng::seq::indexed_char(4) == 'N');
//    BOOST_CHECK(dng::seq::indexed_char(0) == 'A');
//    BOOST_CHECK(dng::seq::indexed_char(1) == 'C');
//    BOOST_CHECK(dng::seq::indexed_char(2) == 'G');
//    BOOST_CHECK(dng::seq::indexed_char(3) == 'T');
//
//}
