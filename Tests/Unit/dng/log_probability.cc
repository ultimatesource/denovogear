/*
 * Copyright (c) 2016 Steven H. Wu
 * Copyright (c) 2017-2018 Reed A. Cartwright
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
 *           Reed A. Cartwright <reed@cartwrig.ht>
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

#define BOOST_TEST_MODULE dng::log_probability

#include <dng/probability.h>

#include <dng/stats.h>
#include <dng/detail/rangeio.h>
#include <dng/task/call.h>

#include <iostream>
#include <fstream>
#include <numeric>


#include "../testing.h"
#include "../xorshift64.h"

namespace dng {
    struct unittest_dng_log_probability {
        GETTER2(LogProbability, params, theta)
        GETTER2(LogProbability, params, ref_bias_hom)
        GETTER2(LogProbability, params, ref_bias_het)
        GETTER2(LogProbability, params, ref_bias_hap)
        GETTER2(LogProbability, params, over_dispersion_hom)
        GETTER2(LogProbability, params, over_dispersion_het)
        GETTER2(LogProbability, params, sequencing_bias)
        GETTER2(LogProbability, params, error_rate)
        GETTER2(LogProbability, params, lib_k_alleles)
        GETTER2(LogProbability, params, k_alleles)

        GETTERF(LogProbability, HaploidPrior)
        GETTERF(LogProbability, DiploidPrior)

        GETTER1(LogProbability, work)
    };
}

using u = dng::unittest_dng_log_probability;

using namespace dng;
using namespace dng::detail;
using Sex = dng::Pedigree::Sex;

int g_seed_counter = 0;

// Use a lambda function to construct a global relationship graph
const auto g_rel_graph = []() {
    libraries_t libs = {
        {"Mom", "Dad", "Eve"},
        {"Mom", "Dad", "Eve"}
    };
    Pedigree ped;
    ped.AddMember("Dad","0","0",Sex::Male,"");
    ped.AddMember("Mom","0","0",Sex::Female,"");
    ped.AddMember("Eve","Dad","Mom",Sex::Female,"");

    RelationshipGraph g;
    g.Construct(ped, libs, 1e-8, 1e-8, 1e-8,true);
    return g;
}();

const auto g_params = []() {
    LogProbability::params_t params;
    params.theta = 0.001;
    params.ref_bias_hom = 0.01;
    params.ref_bias_het = 0.011;
    params.ref_bias_hap = 0.012;
    params.over_dispersion_hom = 1e-4;
    params.over_dispersion_het = 1e-3;
    params.sequencing_bias = 1.1;
    params.error_rate = 2e-4;
    params.lib_k_alleles = 4;
    params.k_alleles = 5;

    return params;
}();

BOOST_AUTO_TEST_CASE(test_constructor_1) {
    LogProbability log_probability{g_rel_graph, g_params};

    BOOST_CHECK_EQUAL(u::theta(log_probability), g_params.theta);
    BOOST_CHECK_EQUAL(u::ref_bias_hom(log_probability), g_params.ref_bias_hom);
    BOOST_CHECK_EQUAL(u::ref_bias_het(log_probability), g_params.ref_bias_het);
    BOOST_CHECK_EQUAL(u::ref_bias_hap(log_probability), g_params.ref_bias_hap);
    
    BOOST_CHECK_EQUAL(u::over_dispersion_hom(log_probability), g_params.over_dispersion_hom);
    BOOST_CHECK_EQUAL(u::over_dispersion_het(log_probability), g_params.over_dispersion_het);
    BOOST_CHECK_EQUAL(u::sequencing_bias(log_probability), g_params.sequencing_bias);
    BOOST_CHECK_EQUAL(u::error_rate(log_probability), g_params.error_rate);
    BOOST_CHECK_EQUAL(u::lib_k_alleles(log_probability), g_params.lib_k_alleles);
    BOOST_CHECK_EQUAL(u::k_alleles(log_probability), g_params.k_alleles);
}

BOOST_AUTO_TEST_CASE(test_work) {
    LogProbability log_probability{g_rel_graph, g_params};
    BOOST_CHECK_EQUAL(&u::work(log_probability), &log_probability.work());
}

BOOST_AUTO_TEST_CASE(test_haploid_prior) {
    const double prec = 2.0*DBL_EPSILON;
    LogProbability log_probability{g_rel_graph, g_params};

    auto test = [prec, &log_probability](int num_obs_alleles) {
        BOOST_TEST_CONTEXT("num_obs_alleles=" << num_obs_alleles << " has_ref=true") {
            auto test_ = u::HaploidPrior(log_probability, num_obs_alleles, true);
            auto expected_ = population_prior_haploid(num_obs_alleles, g_params.theta,
                g_params.ref_bias_hap, g_params.k_alleles, true);
            auto test_range = make_test_range(test_);
            auto expected_range = make_test_range(expected_);
            CHECK_CLOSE_RANGES(test_range, expected_range, prec);
        }
        num_obs_alleles += 1;
        BOOST_TEST_CONTEXT("num_obs_alleles=" << num_obs_alleles << " has_ref=false") {
            auto test_ = u::HaploidPrior(log_probability, num_obs_alleles, false);
            auto expected_ = population_prior_haploid(num_obs_alleles, g_params.theta,
                g_params.ref_bias_hap, g_params.k_alleles, false);
            auto test_range = make_test_range(test_);
            auto expected_range = make_test_range(expected_);
            CHECK_CLOSE_RANGES(test_range, expected_range, prec);
        }
    };
    for(int i=1;i<=10;++i) {
        test(i);
    }
}

BOOST_AUTO_TEST_CASE(test_diploid_prior) {
    const double prec = 2.0*DBL_EPSILON;
    LogProbability log_probability{g_rel_graph, g_params};

    auto test = [prec, &log_probability](int num_obs_alleles) {
        BOOST_TEST_CONTEXT("num_obs_alleles=" << num_obs_alleles << " has_ref=true") {
            auto test_ = u::DiploidPrior(log_probability, num_obs_alleles, true);
            auto expected_ = population_prior_diploid(num_obs_alleles, g_params.theta,
                g_params.ref_bias_hom, g_params.ref_bias_het, g_params.k_alleles, true);
            auto test_range = make_test_range(test_);
            auto expected_range = make_test_range(expected_);
            CHECK_CLOSE_RANGES(test_range, expected_range, prec);
        }
        num_obs_alleles += 1;
        BOOST_TEST_CONTEXT("num_obs_alleles=" << num_obs_alleles << " has_ref=false") {
            auto test_ = u::DiploidPrior(log_probability, num_obs_alleles, false);
            auto expected_ = population_prior_diploid(num_obs_alleles, g_params.theta,
                g_params.ref_bias_hom, g_params.ref_bias_het, g_params.k_alleles, false);
            auto test_range = make_test_range(test_);
            auto expected_range = make_test_range(expected_);
            CHECK_CLOSE_RANGES(test_range, expected_range, prec);
        }
    };
    for(int i=1;i<=10;++i) {
        test(i);
    }
}

BOOST_AUTO_TEST_CASE(test_calcualte_ldd_trio_autosomal) {
    //using ad_t = dng::pileup::allele_depths_t;
    using ad_t = std::vector<std::vector<int>>;
    using dng::utility::make_array;

    xorshift64 xrand(++g_seed_counter);

    auto rexp = [&](double mean) {
        double d = -log(xrand.get_double52())*mean;
        if(xrand.get_double52() > 0.9) {
            d = 0.0;
        }
        return static_cast<int>(d);
    };

    const double prec = FLT_EPSILON;

    libraries_t libs = {
        {"Mom", "Dad", "Eve"},
        {"Mom", "Dad", "Eve"}
    };
    Pedigree ped;
    ped.AddMember("Dad","0","0",Sex::Male,"");
    ped.AddMember("Mom","0","0",Sex::Female,"");
    ped.AddMember("Eve","Dad","Mom",Sex::Female,"");

    const double mu_g = 1e-8;
    const double mu_s = 1e-8;
    const double mu_l = 1e-8;
    RelationshipGraph graph;
    graph.Construct(ped, libs, InheritanceModel::Autosomal, mu_g, mu_s, mu_l, true);

    LogProbability log_probability{graph, g_params};
    int counter = 0;
    auto test = [&](const ad_t& depths, int num_obs_alleles, bool has_ref) {
        BOOST_TEST_CONTEXT("mom=" << rangeio::wrap(depths[0]) << 
                         ", dad=" << rangeio::wrap(depths[1]) <<
                         ", eve=" << rangeio::wrap(depths[2]) <<
                         ", num_obs_alleles=" << num_obs_alleles <<
                         ", has_ref=" << (int)has_ref
        ) {
        int num_obs_alleles_old = num_obs_alleles;
        if(num_obs_alleles > 4) {
            num_obs_alleles = 4;
        }

        LogProbability::value_t test_value;
        double expected_log_scale, expected_log_data;

        Genotyper genotyper{
            g_params.over_dispersion_hom,
            g_params.over_dispersion_het,
            g_params.sequencing_bias,
            g_params.error_rate,
            g_params.lib_k_alleles
        };

        const int hap_sz = num_obs_alleles;
        const int gt_sz = hap_sz*(hap_sz+1)/2;

        auto hap_prior = population_prior_haploid(num_obs_alleles, g_params.theta,
                g_params.ref_bias_hap, g_params.k_alleles, has_ref);
        auto dip_prior = population_prior_diploid(num_obs_alleles, g_params.theta,
                g_params.ref_bias_hom, g_params.ref_bias_het, g_params.k_alleles, has_ref);

        GenotypeArray mom_lower, dad_lower, eve_lower;
        expected_log_scale = 0.0;
        expected_log_scale += genotyper(depths[0], num_obs_alleles, 2, &mom_lower);
        expected_log_scale += genotyper(depths[1], num_obs_alleles, 2, &dad_lower);
        expected_log_scale += genotyper(depths[2], num_obs_alleles, 2, &eve_lower);

        BOOST_REQUIRE_EQUAL(mom_lower.size(), gt_sz);
        BOOST_REQUIRE_EQUAL(dad_lower.size(), gt_sz);
        BOOST_REQUIRE_EQUAL(eve_lower.size(), gt_sz);

        const auto m_sl = Mk::matrix(hap_sz, mu_s+mu_l, g_params.k_alleles);
        const auto m_gsl = Mk::matrix(hap_sz, mu_g+mu_s+mu_l, g_params.k_alleles);

        const auto m_eve = meiosis_matrix(2, m_gsl, 2, m_gsl);
        const auto m_dad = mitosis_matrix(2, m_sl);
        const auto m_mom = mitosis_matrix(2, m_sl);

        dng::stats::ExactSum lld;
        for(int dad_g=0,par=0;dad_g<gt_sz;++dad_g) {
            for(int mom_g=0; mom_g<gt_sz;++mom_g,++par) {
                double lld_eve = 0.0, lld_mom = 0.0, lld_dad = 0.0;
                for(int eve=0;eve<gt_sz;++eve) {
                    lld_eve += eve_lower(eve)*m_eve(par,eve);
                }
                for(int mom=0;mom<gt_sz;++mom) {
                    lld_mom += mom_lower(mom)*m_mom(mom_g,mom);
                }
                for(int dad=0;dad<gt_sz;++dad) {
                    lld_dad += dad_lower(dad)*m_dad(dad_g,dad);
                }
                lld += lld_eve*lld_mom*lld_dad*dip_prior(mom_g)*dip_prior(dad_g);
            }
        }

        expected_log_scale /= M_LN10;
        expected_log_data = log(lld.result()) / M_LN10;

        test_value = log_probability.CalculateLLD(depths, num_obs_alleles_old, has_ref);
        BOOST_CHECK_CLOSE_FRACTION(test_value.log_scale, expected_log_scale, prec);
        BOOST_CHECK_CLOSE_FRACTION(test_value.log_data, expected_log_data, prec);

        test_value = log_probability(depths, num_obs_alleles_old, has_ref);
        BOOST_CHECK_CLOSE_FRACTION(test_value.log_scale, expected_log_scale, prec);
        BOOST_CHECK_CLOSE_FRACTION(test_value.log_data, expected_log_data, prec);
    }};

    test({{0},{0},{0}}, 1, true);
    test({{10},{5},{15}}, 1, true);
    
    test({{100,0},{100,0},{50,50}}, 2, true);
    test({{100,0},{50,50},{100,0}}, 2, true);
    test({{50,50},{100,0},{100,0}}, 2, true);
    
    test({{100,0,0},{50,50,0},{50,50,1}}, 3, true);
    
    test({{0,100},{0,100},{0,10}}, 2, false);

    for(int n=1; n<=6; ++n) {
        for(int i=0; i<100; ++i) {
            ad_t ad(3);
            for(auto &&a : ad) {
                for(int j=0;j<=n;++j) {
                    a.push_back(rexp(30));
                }
            }
            test(ad, n, true);
        }
    }

    for(int n=2; n<=6; ++n) {
        for(int i=0; i<100; ++i) {
            ad_t ad(3);
            for(auto &&a : ad) {
                a.push_back(0);
                for(int j=1;j<=n;++j) {
                    a.push_back(rexp(30));
                }
            }
            test(ad, n, false);
        }
    }

}

BOOST_AUTO_TEST_CASE(test_calcualte_ldd_trio_xlinked) {
    //using ad_t = dng::pileup::allele_depths_t;
    using ad_t = std::vector<std::vector<int>>;
    using dng::utility::make_array;

    xorshift64 xrand(++g_seed_counter);

    auto rexp = [&](double mean) {
        double d = -log(xrand.get_double52())*mean;
        if(xrand.get_double52() > 0.9) {
            d = 0.0;
        }
        return static_cast<int>(d);
    };

    const double prec = FLT_EPSILON;

    libraries_t libs = {
        {"Mom", "Dad", "Eve"},
        {"Mom", "Dad", "Eve"}
    };
    Pedigree ped;
    ped.AddMember("Dad","0","0",Sex::Male,"");
    ped.AddMember("Mom","0","0",Sex::Female,"");
    ped.AddMember("Eve","Dad","Mom",Sex::Female,"");

    const double mu_g = 1e-8;
    const double mu_s = 1e-8;
    const double mu_l = 1e-8;
    RelationshipGraph graph;
    graph.Construct(ped, libs, InheritanceModel::XLinked, mu_g, mu_s, mu_l, true);

    LogProbability log_probability{graph, g_params};
    int counter = 0;
    auto test = [&](const ad_t& depths, int num_obs_alleles, bool has_ref) {
        BOOST_TEST_CONTEXT("mom=" << rangeio::wrap(depths[0]) << 
                         ", dad=" << rangeio::wrap(depths[1]) <<
                         ", eve=" << rangeio::wrap(depths[2]) <<
                         ", num_obs_alleles=" << num_obs_alleles <<
                         ", has_ref=" << (int)has_ref
        ) {
        int num_obs_alleles_old = num_obs_alleles;
        if(num_obs_alleles > 4) {
            num_obs_alleles = 4;
        }

        LogProbability::value_t test_value;
        double expected_log_scale, expected_log_data;

        Genotyper genotyper{
            g_params.over_dispersion_hom,
            g_params.over_dispersion_het,
            g_params.sequencing_bias,
            g_params.error_rate,
            g_params.lib_k_alleles
        };

        const int hap_sz = num_obs_alleles;
        const int gt_sz = hap_sz*(hap_sz+1)/2;

        auto hap_prior = population_prior_haploid(num_obs_alleles, g_params.theta,
                g_params.ref_bias_hap, g_params.k_alleles, has_ref);
        auto dip_prior = population_prior_diploid(num_obs_alleles, g_params.theta,
                g_params.ref_bias_hom, g_params.ref_bias_het, g_params.k_alleles, has_ref);

        GenotypeArray mom_lower, dad_lower, eve_lower;
        expected_log_scale = 0.0;
        expected_log_scale += genotyper(depths[0], num_obs_alleles, 2, &mom_lower);
        expected_log_scale += genotyper(depths[1], num_obs_alleles, 1, &dad_lower);
        expected_log_scale += genotyper(depths[2], num_obs_alleles, 2, &eve_lower);

        BOOST_REQUIRE_EQUAL(mom_lower.size(), gt_sz);
        BOOST_REQUIRE_EQUAL(dad_lower.size(), hap_sz);
        BOOST_REQUIRE_EQUAL(eve_lower.size(), gt_sz);

        const auto m_sl = Mk::matrix(hap_sz, mu_s+mu_l, g_params.k_alleles);
        const auto m_gsl = Mk::matrix(hap_sz, mu_g+mu_s+mu_l, g_params.k_alleles);

        const auto m_eve = meiosis_matrix(1, m_gsl, 2, m_gsl);
        const auto m_dad = mitosis_matrix(1, m_sl);
        const auto m_mom = mitosis_matrix(2, m_sl);

        dng::stats::ExactSum lld;
        for(int dad_g=0,par=0;dad_g<hap_sz;++dad_g) {
            for(int mom_g=0; mom_g<gt_sz;++mom_g,++par) {
                double lld_eve = 0.0, lld_mom = 0.0, lld_dad = 0.0;
                for(int eve=0;eve<gt_sz;++eve) {
                    lld_eve += eve_lower(eve)*m_eve(par,eve);
                }
                for(int mom=0;mom<gt_sz;++mom) {
                    lld_mom += mom_lower(mom)*m_mom(mom_g,mom);
                }
                for(int dad=0;dad<hap_sz;++dad) {
                    lld_dad += dad_lower(dad)*m_dad(dad_g,dad);
                }
                lld += lld_eve*lld_mom*lld_dad*dip_prior(mom_g)*hap_prior(dad_g);
            }
        }

        expected_log_scale /= M_LN10;
        expected_log_data = log(lld.result()) / M_LN10;

        test_value = log_probability.CalculateLLD(depths, num_obs_alleles_old, has_ref);
        BOOST_CHECK_CLOSE_FRACTION(test_value.log_scale, expected_log_scale, prec);
        BOOST_CHECK_CLOSE_FRACTION(test_value.log_data, expected_log_data, prec);

        test_value = log_probability(depths, num_obs_alleles_old, has_ref);
        BOOST_CHECK_CLOSE_FRACTION(test_value.log_scale, expected_log_scale, prec);
        BOOST_CHECK_CLOSE_FRACTION(test_value.log_data, expected_log_data, prec);
    }};

    test({{0},{0},{0}}, 1, true);
    test({{10},{5},{15}}, 1, true);
    
    test({{100,0},{100,0},{50,50}}, 2, true);
    test({{100,0},{50,50},{100,0}}, 2, true);
    test({{50,50},{100,0},{100,0}}, 2, true);
    
    test({{100,0,0},{50,50,0},{50,50,1}}, 3, true);
    
    test({{0,100},{0,100},{0,10}}, 2, false);

    for(int n=1; n<=6; ++n) {
        for(int i=0; i<100; ++i) {
            ad_t ad(3);
            for(auto &&a : ad) {
                for(int j=0;j<=n;++j) {
                    a.push_back(rexp(30));
                }
            }
            test(ad, n, true);
        }
    }

    for(int n=2; n<=6; ++n) {
        for(int i=0; i<100; ++i) {
            ad_t ad(3);
            for(auto &&a : ad) {
                a.push_back(0);
                for(int j=1;j<=n;++j) {
                    a.push_back(rexp(30));
                }
            }
            test(ad, n, false);
        }
    }

}
