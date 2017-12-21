#define BOOST_TEST_MODULE dng::genotype

#include <dng/genotyper.h>

#include <dng/detail/rangeio.h>

#include "../testing.h"
#include "../xorshift64.h"

#include <boost/range/algorithm/fill.hpp>
#include <boost/range/algorithm/max_element.hpp>

using namespace dng;
using namespace dng::genotype;
using dng::detail::make_test_range;

int g_seed_counter = 0;

const double low_phi = DBL_EPSILON/2.0;
const double low_epsilon = DBL_EPSILON/2.0;

BOOST_AUTO_TEST_CASE(test_log_pochhammer) {
    double prec = 32*DBL_EPSILON;

    auto test = [prec](double a) -> void {
        BOOST_TEST_CONTEXT("a=" << a)
    {
        using dng::genotype::detail::log_pochhammer;
        log_pochhammer pochy(a);
        std::vector<double> test_values, expected_values;
        double d = 0.0;
        for(int i=0; i < 1000; ++i) {
            d += log(i+a);
            expected_values.push_back(d);
            test_values.push_back(pochy(i+1));
        }
        CHECK_CLOSE_RANGES( test_values, expected_values, prec);
    }
    };

    test(1e0);
    test(1e2);
    test(1e4);
    test(1e6);
    test(1e8);

    test(1e-2);
    test(1e-4);
    test(1e-6);
    test(1e-8);

    test((1.0 - low_phi) / low_phi);
}

BOOST_AUTO_TEST_CASE(test_make_alphas) {
    using dng::genotype::detail::make_alphas;

    double prec = DBL_EPSILON;

    auto test = [prec](double over_dispersion_hom, double over_dispersion_het, double ref_bias,
        double error_rate, double error_entropy) {
    BOOST_TEST_CONTEXT("over_dispersion_hom=" << over_dispersion_hom 
                  << ", over_dispersion_het=" << over_dispersion_het
                  << ", ref_bias=" << ref_bias
                  << ", error_rate=" << error_rate
                  << ", error_entropy=" << error_entropy
    ) {
        auto test_range = make_alphas(over_dispersion_hom, over_dispersion_het, ref_bias,
                                      error_rate, error_entropy);
        std::array<double, 8> expected_range;

        if(over_dispersion_hom < low_phi) {
            over_dispersion_hom = low_phi;
        }
        if(over_dispersion_het < low_phi) {
            over_dispersion_het = low_phi;
        }
        if(error_rate < low_epsilon) {
            error_rate = low_epsilon;
        }

        double a1 = (1.0-over_dispersion_hom)/over_dispersion_hom;
        double a2 = (1.0-over_dispersion_het)/over_dispersion_het;

        expected_range[0] = a1*(1.0-error_rate);
        expected_range[1] = a1*error_rate/exp(error_entropy);

        expected_range[2] = a2*(1.0-error_rate + error_rate/exp(error_entropy))*ref_bias/(1.0+ref_bias);
        expected_range[3] = a2*(1.0-error_rate + error_rate/exp(error_entropy))*1.0/(1.0+ref_bias);
        expected_range[4] = a2*error_rate/exp(error_entropy);
        expected_range[5] = a2*(1.0-error_rate + error_rate/exp(error_entropy))*0.5;

        expected_range[6] = a1;
        expected_range[7] = a2;

        CHECK_CLOSE_RANGES(test_range, expected_range, prec);

    }};

    test(0.0005, 0.0005, 1, 0.0005, log(3));
    test(0.0005, 0.001, 1.02, 1e-4, log(4));
    test(0, 0, 1, 0, log(4));
}

BOOST_AUTO_TEST_CASE(test_DirichletMultinomial) {
    using dng::genotype::detail::make_alphas;
    using dng::genotype::detail::log_pochhammer;
    using depths_t = std::vector< std::vector<int> >;

    double prec = 2*DBL_EPSILON;

    xorshift64 xrand(++g_seed_counter);

    auto test = [prec](const DirichletMultinomial& dm, const depths_t &depths) -> void {
    BOOST_TEST_CONTEXT("over_dispersion_hom=" << dm.over_dispersion_hom()
                  << ", over_dispersion_het=" << dm.over_dispersion_het()
                  << ", ref_bias=" << dm.ref_bias()
                  << ", error_rate=" << dm.error_rate()
                  << ", error_entropy=" << dm.error_entropy()
    ){
        auto alphas = make_alphas(dm.over_dispersion_hom(), dm.over_dispersion_het(),
            dm.ref_bias(), dm.error_rate(), dm.error_entropy());
        std::vector<log_pochhammer> f;
        for(auto &&a : alphas) {
            f.emplace_back(a);
        }
        for(int k : {0,1,2,3,4}) {         
            for(auto &&ad : depths) {
                // Wrap depths 
                boost::const_multi_array_ref<int, 2> r{ad.data(), boost::extents[1][ad.size()]};

                BOOST_TEST_CONTEXT("ploidy=1, depths=" << rangeio::wrap(ad) << ", k=" << k) {
                    std::vector<double> expected;
                    for(int g=0;g<=k;++g) {
                        double ll = 0.0;
                        int count = 0;
                        for(int d=0;d<ad.size();++d) {
                            if(d == g) {
                                ll += f[0](ad[d]);
                            } else {
                                ll += f[1](ad[d]);
                            }
                            count += ad[d];
                        }
                        ll -= f[6](count);
                        expected.push_back(ll);
                    }
                    double expected_scale = *boost::max_element(expected);
                    for(auto &&x : expected) {
                         x = exp(x-expected_scale);
                     }
                    auto results = dm(r[0],k,1);
                    double test_scale = results.second;
                    auto test = make_test_range(results.first);
                    BOOST_CHECK_CLOSE_FRACTION(test_scale, expected_scale, prec);
                    CHECK_CLOSE_RANGES(test, expected, prec);
                }
                BOOST_TEST_CONTEXT("ploidy=2, depths=" << rangeio::wrap(ad) << ", k=" << k) {
                    std::vector<double> expected;
                    for(int a=0;a<=k;++a) {
                        for(int b=0;b<=a;++b) {
                            double ll = 0.0;
                            int count = 0;
                            for(int d=0;d<ad.size();++d) {
                                if(a == b) {
                                    // HOM_MATCH or HOM_ERROR
                                    ll += (d == a) ? f[0](ad[d]) : f[1](ad[d]);
                                } else if(d != a && d != b) {
                                    // HET_ERROR
                                    ll += f[4](ad[d]);
                                } else if(a != 0 && b != 0) {
                                    // HET_ALTALT_MATCH
                                    ll += f[5](ad[d]);
                                } else {
                                    ll += (d == 0) ? f[2](ad[d]) : f[3](ad[d]);
                                }
                                count += ad[d];
                            }
                            ll -= (a == b) ? f[6](count) : f[7](count);
                            expected.push_back(ll);
                        }
                    }
                    double expected_scale = *boost::max_element(expected);
                    for(auto &&x : expected) {
                         x = exp(x-expected_scale);
                     }
                    auto results = dm(r[0],k,2);
                    double test_scale = results.second;
                    auto test = make_test_range(results.first);
                    BOOST_CHECK_CLOSE_FRACTION(test_scale, expected_scale, prec);
                    CHECK_CLOSE_RANGES(test, expected, prec);
                }
            }
        }
    }};

    depths_t ad;
    for(int i=0;i<1000;++i) {
        depths_t::value_type d;
        for(int j=0;j<=(i%5);++j) {
            if(xrand.get_double53() < 0.5) {
                d.push_back(xrand.get_uint64(100));
            } else if(xrand.get_double53() < 0.75) {
                d.push_back(xrand.get_uint64(1000));
            } else if(xrand.get_double53() < 0.9) {
                d.push_back(xrand.get_uint64(100000));
            } else {
                d.push_back(0);
            }
        }
        ad.push_back(d);
    }

    test({0.0005, 0.0005, 1, 0.0005, log(3)}, ad);
    test({0.0005, 0.001, 1.02, 1e-4, log(4)}, ad);
    test({0, 0, 1, 0, log(4)}, ad);
}
