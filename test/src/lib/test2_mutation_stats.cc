//
// Created by steven on 2/9/16.
//


#define BOOST_TEST_MODULE dng::lib::mutation_stats

#include <boost/test/unit_test.hpp>
#include <iostream>

#include <dng/mutation_stats.h>

#include <boost_test_helper.h>
#include <fixture/fixture_trio_workspace.h>



const int NUM_TEST = 100;

struct Fixture {
    std::mt19937 random_gen_mt {1}; //seed = 1

    std::string fixture;
    std::uniform_real_distribution<double> rand_unif_log;
    std::uniform_real_distribution<double> rand_unif;

    double min_prob = 0.01;
    Fixture(std::string s = "TestMutationStats") : fixture(s) {

        rand_unif_log = std::uniform_real_distribution<double>(std::log(1e-20), 0);
        rand_unif = std::uniform_real_distribution<double>(0,1);

        BOOST_TEST_MESSAGE("set up fixture: " << fixture);
    }

    ~Fixture() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }

};

void setup() { BOOST_TEST_MESSAGE("set up fun"); }

void teardown() { BOOST_TEST_MESSAGE("tear down fun"); }



BOOST_FIXTURE_TEST_SUITE(test_mutation_stats_suite, Fixture )

BOOST_AUTO_TEST_CASE(test_set_mup, *utf::fixture(&setup, &teardown)) {

    MutationStats stats (min_prob);

    for (int t = 0; t <NUM_TEST; ++t) {

        double ln_mut = rand_unif_log(random_gen_mt);
        double ln_no_mut = rand_unif_log(random_gen_mt);

        double prob_mu = std::exp(ln_mut);
        double prob_no_mut = std::exp(ln_no_mut);

        double expected = 1 - (prob_no_mut / prob_mu);
        double expected_log =  - std::expm1(ln_no_mut - ln_mut);

        stats.set_mutation_prob(ln_no_mut, ln_mut);
        BOOST_CHECK_CLOSE(expected, expected_log, 0.01);
        BOOST_CHECK_CLOSE(expected_log, stats.get_mutation_prob(), BOOST_CLOSE_THRESHOLD);
    }

    for (int t = 0; t <NUM_TEST; ++t) {

        double prob_mu = rand_unif(random_gen_mt);
        double prob_no_mut = rand_unif(random_gen_mt);

        double ln_mu = std::log(prob_mu);
        double ln_no_mut = std::log(prob_no_mut);

        double expected = 1 - (prob_no_mut / prob_mu);

        stats.set_mutation_prob(ln_no_mut, ln_mu);
        BOOST_CHECK_CLOSE(expected, stats.get_mutation_prob(), BOOST_CLOSE_THRESHOLD);
    }

}


BOOST_AUTO_TEST_CASE(test_set_node, *utf::fixture(&setup, &teardown)) {

    MutationStats stats(min_prob);
    int number_of_event = 20;
    int first_non_founder = 5;

    std::vector<double> event (number_of_event, hts::bcf::float_missing);

    for (int t = first_non_founder; t < event.size(); ++t) {
        event[t] = t;
    }

    double total = 0;
    for (int i = first_non_founder; i < event.size(); ++i) {
        total += event[i];
    }

    std::vector<float> expected_mu1p (number_of_event, hts::bcf::float_missing);
    for (int i = first_non_founder; i < event.size(); ++i) {
        expected_mu1p[i] = static_cast<float>(event[i]/total);
    }

    stats.set_node_mup(event, first_non_founder);
    stats.set_node_mu1p(event, total, first_non_founder);

    for (int i = 0; i < first_non_founder; ++i) {
        BOOST_CHECK( bcf_float_is_missing(stats.node_mup[i]) );
        BOOST_CHECK( bcf_float_is_missing(stats.node_mu1p[i]));
    }

    for (int i = first_non_founder; i < event.size(); ++i) {
        BOOST_CHECK_EQUAL(i, stats.node_mup[i]);
        BOOST_CHECK_EQUAL(expected_mu1p[i], stats.node_mu1p[i]);
    }

}


BOOST_AUTO_TEST_CASE(test_record, *utf::fixture(&setup, &teardown)) {
    BOOST_TEST_MESSAGE("Not yet implemented!");
    //TODO: implement check on record_info/stats, hts::bcf::Variant
    // record_single_mutation_stats(hts::bcf::Variant &record){
}


BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE(test_mutation_stats_suite2, TrioWorkspace )

BOOST_AUTO_TEST_CASE(test_set_posterior_probabilities, *utf::fixture(&setup, &teardown)) {



    dng::IndividualVector expected_probs;
    expected_probs.reserve(workspace.num_nodes);

    for (int i = 0; i < workspace.num_nodes; ++i) {
        double sum = 0;
        for (int j = 0; j < 10; ++j) {
            expected_probs[i][j] = workspace.upper[i][j] * workspace.lower[i][j];
            sum += expected_probs[i][j];
        }
        expected_probs[i] /= sum;
    }

    MutationStats stats (min_prob);
    stats.set_posterior_probabilities(workspace);
    for (int i = 0; i < workspace.num_nodes; ++i) {
        boost_check_close_vector(expected_probs[i], stats.inspect_posterior_at(i));

    }

}


BOOST_AUTO_TEST_CASE(test_genotype_stats, *utf::fixture(&setup, &teardown)) {


    FindMutationsGetter find_mutation{min_prob, pedigree, test_param_1};
    MutationStats stats(min_prob);

    find_mutation.calculate_mutation(read_depths, ref_index, stats);

    const int acgt_to_refalt_allele[] = {-1, 2, 0, 1, -1};
    const int refalt_to_acgt_allele[5] = {2, 3, 1, -1, -1};
    const uint32_t n_alleles = 3;
    const std::size_t ref_index = 2;
    const std::size_t num_nodes = 5;
    const std::size_t library_start = 2;

    std::vector<float> expected_gp_scores{
            0.999833, 0.000167058, 1.78739e-12, 2.75325e-11, 3.57372e-16, 4.51356e-20,
            0.984574, 0.0154262, 1.66835e-10, 9.77602e-12, 1.67414e-14, 1.58815e-20,
            0.868225, 0.131775, 7.45988e-16, 1.99601e-09, 2.97932e-17, 3.704e-26,
            0.999837, 0.000163186, 2.11375e-15, 2.7779e-11, 4.22658e-19, 5.33792e-23,
            0.984578, 0.0154223, 6.97054e-14, 9.45085e-12, 6.98069e-18, 2.94856e-24
    };
    std::vector<float> expected_gl_scores{
            NAN, NAN, NAN, NAN, NAN, NAN,
            NAN, NAN, NAN, NAN, NAN, NAN,
            -6.74046,0,-6.07991,-7.58973,-7.5555,-15.93,
            0,-6.64246,-17.5301,-6.64246,-17.5301,-17.5301,
            0,-4.667,-16.0119,-7.10524,-16.3122,-18.7879
    };
    std::vector<int> expected_best_genotype {2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    std::vector<int> expected_genotype_qualities {38, 18, 9, 38, 18};


    stats.set_genotype_related_stats(acgt_to_refalt_allele, refalt_to_acgt_allele,
                                     n_alleles, ref_index, num_nodes, library_start);

    boost_check_close_vector(expected_gp_scores, stats.gp_scores);
    boost_check_equal_vector(expected_best_genotype, stats.best_genotypes);
    boost_check_equal_vector(expected_genotype_qualities, stats.genotype_qualities);

    BOOST_CHECK_EQUAL(expected_gl_scores.size(), stats.gl_scores.size());
    for (int i = 0; i < expected_gl_scores.size(); ++i) {
        if( isnan(expected_gl_scores[i]) ){
            BOOST_CHECK( bcf_float_is_missing(stats.gl_scores[i]) );
        }
        else{
            BOOST_CHECK_CLOSE(expected_gl_scores[i], stats.gl_scores[i], BOOST_CLOSE_THRESHOLD);
        }
    }


}

BOOST_AUTO_TEST_SUITE_END()

