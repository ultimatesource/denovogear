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
using d4 = std::array<double,4>;

int g_seed_counter = 0;

BOOST_AUTO_TEST_CASE(test_diploid_nucleotides) {
    for(int a : {0,1,2,3}) {
        for(int b : {0,1,2,3}) {
            BOOST_TEST_CONTEXT("a=" << a << ", b=" << b) {
                int gt = folded_diploid_genotypes_matrix[a][b];
                
                BOOST_CHECK((folded_diploid_nucleotides[gt][0] == a && folded_diploid_nucleotides[gt][1] == b) ||
                    folded_diploid_nucleotides[gt][0] == b && folded_diploid_nucleotides[gt][1] == a);
                int unfolded = a*4+b;
                int gt2 = folded_diploid_genotypes[unfolded];
                BOOST_CHECK((folded_diploid_nucleotides[gt2][0] == a && folded_diploid_nucleotides[gt2][1] == b) ||
                    folded_diploid_nucleotides[gt2][0] == b && folded_diploid_nucleotides[gt2][1] == a);

                if(a <= b) {
                    BOOST_CHECK_EQUAL(unfolded_diploid_genotypes_upper[gt2], unfolded);
                }
                if(b <= a) {
                    BOOST_CHECK_EQUAL(unfolded_diploid_genotypes_lower[gt2], unfolded);
                }
            }
        }
    }
}

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
    using boost::fill;

    double prec = DBL_EPSILON;

    auto test = [prec](double phi, double epsilon, double omega) -> void {
    BOOST_TEST_CONTEXT("phi=" << phi << ", epsilon=" << epsilon << ", omega=" << omega) {
        phi = (phi >= low_phi) ? phi : low_phi;
        epsilon = (epsilon >= low_epsilon) ? epsilon : low_epsilon;

        double alpha = (1.0 - phi) / phi;
        double e = epsilon/3.0;
        double m = alpha*(1.0 - 3.0 * e);
        double h = alpha*(1.0 - 2.0 * e);
        double href = h*omega/(1.0+omega);
        double halt = h*1.0/(1.0+omega);
        e = alpha*e;

        for(int ref : {0,1,2,3,4}) {
            for(int a=0,gt=0;a<4;++a) { 
                for(int b=0;b<=a;++b,++gt) {
                BOOST_TEST_CONTEXT("ref=" << ref << ", gt=" << gt) {
                    auto test_values = make_alphas(ref, gt, phi, epsilon, omega);
                    std::array<double,5> expected_values;
                    fill(expected_values, e);
                    if(a == b) {
                        expected_values[a] = m;
                    } else if(a == ref) {
                        expected_values[a] = href;
                        expected_values[b] = halt;
                    } else if(b == ref) {
                        expected_values[b] = href;
                        expected_values[a] = halt;                       
                    } else {
                        expected_values[a] = h/2.0;
                        expected_values[b] = h/2.0;
                    }
                    expected_values[4] = alpha;

                    CHECK_CLOSE_RANGES(test_values, expected_values, prec);
                }
            }
            }
        }
    }
    };

    test(0.0005, 0.0005, 1.02);
    test(0.0, 0.0, 1);
    test(0.1, 0.01, 1.1);
    test(0.001, 1e-4, 0.8);
}

BOOST_AUTO_TEST_CASE(test_DirichletMultinomial_RawDepths) {
    using dng::genotype::detail::make_alphas;
    using dng::genotype::detail::log_pochhammer;

    double prec = 2*DBL_EPSILON;

    xorshift64 xrand(++g_seed_counter);

    auto test = [prec](const DirichletMultinomial& dm, const pileup::RawDepths& depths) -> void {
    BOOST_TEST_CONTEXT("over_dispersion=" << dm.parameters().over_dispersion
         << ", error_rate=" << dm.parameters().error_rate
         << ", ref_bias=" << dm.parameters().ref_bias) {

        for(int ref : {0,1,2,3,4}) {
            std::vector<std::array<double,5>> alphas;
            for(int gt=0;gt<10;++gt) {
                alphas.push_back(make_alphas(ref, gt,
                    dm.parameters().over_dispersion, dm.parameters().error_rate,
                    dm.parameters().ref_bias));
            }
            for(size_t pos=0; pos < depths.size(); ++pos) {
                std::array<double,10> log_like;
                boost::fill(log_like, 0);
                for(int gt=0;gt<10;++gt) {
                    int total = 0.0;
                    for(size_t i=0;i<4;++i) {
                        int n = depths[pos].counts[i];
                        total += n;
                        log_like[gt] += log_pochhammer{alphas[gt][i]}(n);
                    }
                    log_like[gt] -= log_pochhammer{alphas[gt][4]}(total);
                }

                // ploidy = 2
                BOOST_TEST_CONTEXT("ploidy=2, ref=" << ref
                    << ", depths=" << rangeio::wrap(depths[pos].counts)) {
                    auto expected = log_like; 
                    double expected_scale = *boost::max_element(expected);
                    for(auto &&x : expected) {
                        x = exp(x-expected_scale);
                    }
                    auto results = dm(depths, pos, ref, 2);
                    double test_scale = results.second;
                    auto test = make_test_range(results.first);
                    BOOST_CHECK_CLOSE_FRACTION(test_scale, expected_scale, prec);
                    CHECK_CLOSE_RANGES(test, expected, prec);
                }
                // ploidy = 1
                BOOST_TEST_CONTEXT("ploidy=1, ref=" << ref
                    << ", depths=" << rangeio::wrap(depths[pos].counts)) {
                    std::array<double,4> expected = {
                        log_like[0], log_like[2],
                        log_like[5], log_like[9]
                    }; 
                    double expected_scale = *boost::max_element(expected);
                    for(auto &&x : expected) {
                        x = exp(x-expected_scale);
                    }
                    auto results = dm(depths, pos, ref, 1);
                    double test_scale = results.second;
                    auto test = make_test_range(results.first);
                    BOOST_CHECK_CLOSE_FRACTION(test_scale, expected_scale, prec);
                    CHECK_CLOSE_RANGES(test, expected, prec);
                }
            }
        }
    }
    };

    pileup::RawDepths data;
    for(int i=0;i<1000;++i) {
        pileup::depth_t d;
        for(int i=0;i<4;++i) {
            if(xrand.get_double53() < 0.5) {
                d.counts[i] = xrand.get_uint64(100);
            } else if(xrand.get_double53() < 0.75) {
                d.counts[i] = xrand.get_uint64(1000);
            } else if(xrand.get_double53() < 0.9) {
                d.counts[i] = xrand.get_uint64(100000);
            } else {
                d.counts[i] = 0;
            }
        }
        data.push_back(d);
    }

    test({0.0005, 0.0005, 1.02}, data);
    test({0.0, 0.0, 1}, data);
    test({0.1, 0.01, 1.1}, data);
    test({0.001, 1e-4, 0.8}, data);
}

BOOST_AUTO_TEST_CASE(test_DirichletMultinomial_AlleleDepths) {
    using dng::pileup::AlleleDepths;
    using dng::pileup::RawDepths;

    // The order of addition is different when working
    // on AlleleDepths than RawDepths. This causes loss
    // of precision when counts vary in magnitude.
    // Only validate to FLT_EPSILON precision.
    double prec = FLT_EPSILON;

    xorshift64 xrand(++g_seed_counter);

    auto test = [prec,&xrand](const DirichletMultinomial& dm) -> void {
    BOOST_TEST_CONTEXT("over_dispersion=" << dm.parameters().over_dispersion
         << ", error_rate=" << dm.parameters().error_rate
         << ", ref_bias=" << dm.parameters().ref_bias) {
        for(int color=0; color < AlleleDepths::type_info_table_length; ++color) {
            AlleleDepths data;
            RawDepths raw;
            data.resize(color, 10);
            raw.resize(10);
            for(auto &&a : data.data()) {
                if(xrand.get_double53() < 0.5) {
                    a = xrand.get_uint64(100);
                } else if(xrand.get_double53() < 0.75) {
                    a = xrand.get_uint64(1000);
                } else if(xrand.get_double53() < 0.9) {
                    a = xrand.get_uint64(10000);
                } else {
                    a = 0;
                }
            }
            for(size_t pos=0; pos < data.num_libraries(); ++pos) {
                boost::fill(raw[pos].counts, 0);
                for(int i=0; i < data.num_nucleotides(); ++i) {
                    raw[pos].counts[data.type_info().indexes[i]] = data(pos, i);
                }
            }

            for(size_t pos=0; pos < data.num_libraries(); ++pos) {
                // ploidy=2
                BOOST_TEST_CONTEXT("ploidy=2, color=" << color
                    << ", raw_depths=" << rangeio::wrap(raw[pos].counts)) {

                    auto test_results = dm(data, pos, 2);
                    double test_scale = test_results.second;
                    auto test = make_test_range(test_results.first);

                    auto expected_results = dm(raw, pos, data.type_info().reference, 2);
                    double expected_scale = expected_results.second;
                    std::vector<double> expected(data.type_gt_info().width);
                    for(int i=0; i < expected.size(); ++i) {
                        expected[i] = expected_results.first[data.type_gt_info().indexes[i]];
                    }
                    BOOST_CHECK_CLOSE_FRACTION(test_scale, expected_scale, prec);
                    CHECK_CLOSE_RANGES(test, expected, prec);
                }
                // ploidy=1
                BOOST_TEST_CONTEXT("ploidy=1, color=" << color
                    << ", raw_depths=" << rangeio::wrap(raw[pos].counts)) {

                    auto test_results = dm(data, pos, 1);
                    double test_scale = test_results.second;
                    auto test = make_test_range(test_results.first);

                    auto expected_results = dm(raw, pos, data.type_info().reference, 1);
                    double expected_scale = expected_results.second;
                    std::vector<double> expected(data.type_info().width);
                    for(int i=0; i < expected.size(); ++i) {
                        expected[i] = expected_results.first[data.type_info().indexes[i]];
                    }

                    BOOST_CHECK_CLOSE_FRACTION(test_scale, expected_scale, prec);
                    CHECK_CLOSE_RANGES(test, expected, prec);
                }
            }
        }
    }
    };

    test({0.0005, 0.0005, 1.02});
}
