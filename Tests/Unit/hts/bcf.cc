/*
 * Copyright (c) 2018 Reed A. Cartwright
 * Copyright (c) 2018 Juan J. Garcia Mesa
 * Authors: Reed A. Cartwright <reed@cartwrig.ht>
 *      Juan J. Garcia Mesa <jgarc111@asu.edu>
 *
 * This file is part of DeNovoGear.
 *
 * DeNovoGear is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distrubited in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY of FITNESS
 * FOR ANY PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_MODULE hts::bcf

#include <dng/hts/bcf.h>
#include <dng/hts/hts.h>

#include "../testing.h"

using namespace hts;
using namespace hts::bcf;
using hts::detail::make_data_url;
using dng::detail::make_test_range;
using dng::detail::make_unique_ptr;

BOOST_AUTO_TEST_CASE(test_vcf_special_values) {
    BOOST_CHECK_EQUAL(bcf::int8_missing, bcf_int8_missing);
    BOOST_CHECK_EQUAL(bcf::int8_vector_end, bcf_int8_vector_end);
    BOOST_CHECK_EQUAL(bcf::int16_missing, bcf_int16_missing);
    BOOST_CHECK_EQUAL(bcf::int16_vector_end, bcf_int16_vector_end);
    BOOST_CHECK_EQUAL(bcf::int32_missing, bcf_int32_missing);
    BOOST_CHECK_EQUAL(bcf::int32_vector_end, bcf_int32_vector_end);
    BOOST_CHECK_EQUAL(bcf::int32_missing, bcf_int32_missing);
    BOOST_CHECK_EQUAL(bcf::int32_vector_end, bcf_int32_vector_end);
    BOOST_CHECK(bcf_float_is_missing(bcf::float_missing));
    BOOST_CHECK(bcf_float_is_vector_end(bcf::float_vector_end));
}

const char vcfformatint[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=Integer>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1,2,3\t4,5,6\n"
;

BOOST_AUTO_TEST_CASE(test_vcf_read_int_PL) {
    int counter = 0;
    auto test = [&](const std::string &vcfstr, const std::vector<int> &expected) -> void {
    BOOST_TEST_CONTEXT("counter=" << ++counter) {
        // Set header and rec for reading PL values
        auto file = make_unique_ptr(hts_open(vcfstr.c_str(), "r"), &hts_close);
        auto header = make_unique_ptr(bcf_hdr_read(file.get()), &bcf_hdr_destroy);
        auto record = make_unique_ptr(bcf_init(), &bcf_destroy);
        bcf_read(file.get(), header.get(), record.get());

        int n_pl = 8;
        buffer_t<int32_t> buffer = hts::bcf::make_buffer<int32_t>(n_pl);

        // Get PL values
        int n = hts::bcf::get_format_int32(header.get(), record.get(), "PL", &buffer, &n_pl);
        BOOST_REQUIRE_GE(n,0);
        auto test = make_test_range(buffer.get(),buffer.get()+n);
        CHECK_EQUAL_RANGES(test, expected);
    }};

    test(make_data_url(vcfformatint), {1,2,3,4,5,6});
}

const char vcfformatfloat[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=Float>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1.1,2.1,3.1\t4.1,5.1,6.1\n"
;

BOOST_AUTO_TEST_CASE(test_vcf_read_float_PL) {
    int counter = 0;
    auto test = [&](const std::string &vcfstr, const std::vector<float> &expected) -> void {
    BOOST_TEST_CONTEXT("counter=" << ++counter) {
        // Set header and rec for reading PL values
        auto file = make_unique_ptr(hts_open(vcfstr.c_str(), "r"), &hts_close);
        auto header = make_unique_ptr(bcf_hdr_read(file.get()), &bcf_hdr_destroy);
        auto record = make_unique_ptr(bcf_init(), &bcf_destroy);
        bcf_read(file.get(), header.get(), record.get());

        int n_pl = 8;
        buffer_t<float> buffer = hts::bcf::make_buffer<float>(n_pl);

        // Get PL values
        int n = hts::bcf::get_format_float(header.get(), record.get(), "PL", &buffer, &n_pl);
        BOOST_REQUIRE_GE(n,0);
        auto test = make_test_range(buffer.get(),buffer.get()+n);
        CHECK_EQUAL_RANGES(test, expected);
    }};

    test(make_data_url(vcfformatfloat), {1.1,2.1,3.1,4.1,5.1,6.1});
}

const char vcfupdateint[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=Integer>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20\n"
;

const char vcfupdatefloat[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=Float>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0\n"
;

BOOST_AUTO_TEST_CASE(test_vcf_read_numeric_int_PL) {
    int counter = 0;
    auto test = [&](const std::string &vcfstr, const std::vector<int> &expected) -> void {
    BOOST_TEST_CONTEXT("counter=" << ++counter) {
        // Set header and rec for reading PL values
        auto file = make_unique_ptr(hts_open(vcfstr.c_str(), "r"), &hts_close);
        auto header = make_unique_ptr(bcf_hdr_read(file.get()), &bcf_hdr_destroy);
        auto record = make_unique_ptr(bcf_init(), &bcf_destroy);
        bcf_read(file.get(), header.get(), record.get());

        int n_pl = 8;
        buffer_t<int> buffer = hts::bcf::make_buffer<int>(n_pl);

        // Get PL values
        int n = hts::bcf::get_format_numeric(header.get(), record.get(), "PL", &buffer, &n_pl);
        BOOST_REQUIRE_GE(n,0);
        auto test = make_test_range(buffer.get(),buffer.get()+n);
        CHECK_EQUAL_RANGES(test, expected);
    }};

    test(make_data_url(vcfformatint), {1,2,3,4,5,6});
    test(make_data_url(vcfformatfloat), {1,2,3,4,5,6});
    test(make_data_url(vcfupdateint), {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20});
    test(make_data_url(vcfupdatefloat), {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20});
}

BOOST_AUTO_TEST_CASE(test_vcf_read_numeric_float_PL) {
    int counter = 0;
    auto test = [&](const std::string &vcfstr, const std::vector<float> &expected) -> void {
    BOOST_TEST_CONTEXT("counter=" << ++counter) {
        // Set header and rec for reading PL values
        auto file = make_unique_ptr(hts_open(vcfstr.c_str(), "r"), &hts_close);
        auto header = make_unique_ptr(bcf_hdr_read(file.get()), &bcf_hdr_destroy);
        auto record = make_unique_ptr(bcf_init(), &bcf_destroy);
        bcf_read(file.get(), header.get(), record.get());

        int n_pl = 8;
        buffer_t<float> buffer = hts::bcf::make_buffer<float>(n_pl);

        // Get PL values
        int n = hts::bcf::get_format_numeric(header.get(), record.get(), "PL", &buffer, &n_pl);
        BOOST_REQUIRE_GE(n,0);
        auto test = make_test_range(buffer.get(),buffer.get()+n);
        CHECK_EQUAL_RANGES(test, expected);
    }};

    test(make_data_url(vcfformatint), {1.0,2.0,3.0,4.0,5.0,6.0});
    test(make_data_url(vcfformatfloat), {1.1,2.1,3.1,4.1,5.1,6.1});
    test(make_data_url(vcfupdateint), {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0});
    test(make_data_url(vcfupdatefloat), {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0});
}

BOOST_AUTO_TEST_CASE(test_vcf_read_PL_tag_fail) {
    auto test = [](const std::string &vcfstr, const char *tag, int expected) -> void {
    BOOST_TEST_CONTEXT("tag=" << tag) {
    // Set header and rec for reading PL values
        auto file = make_unique_ptr(hts_open(vcfstr.c_str(), "r"), &hts_close);
        auto header = make_unique_ptr(bcf_hdr_read(file.get()), &bcf_hdr_destroy);
        auto record = make_unique_ptr(bcf_init(), &bcf_destroy);
        bcf_read(file.get(), header.get(), record.get());

        int n_pl = 8;//10*record->n_sample;
        buffer_t<float> buffer_int = hts::bcf::make_buffer<float>(n_pl);
        buffer_t<int> buffer_float = hts::bcf::make_buffer<int>(n_pl);
        // Get PL values
        int n = hts::bcf::get_format_numeric(header.get(), record.get(), tag, &buffer_int, &n_pl);
        BOOST_REQUIRE_EQUAL(n,expected);
        n = hts::bcf::get_format_numeric(header.get(), record.get(), tag, &buffer_float, &n_pl);
        BOOST_REQUIRE_EQUAL(n,expected);
    }};

    test(make_data_url(vcfformatint), "Pl",-1);
    test(make_data_url(vcfformatint), "pL",-1);
    test(make_data_url(vcfformatint), "pl",-1);
    test(make_data_url(vcfformatint), "Pl",-1);
}

const char vcftypechar[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=Character>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1,2,3,4\n"
;

const char vcftypestring[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=String>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1,2,3,4\n"
;

BOOST_AUTO_TEST_CASE(test_vcf_read_PL_type_fail) {
    auto test = [](const std::string &vcfstr, int expected) -> void {
    BOOST_TEST_CONTEXT(vcfstr) {
        // Set header and rec for reading PL values
        auto file = make_unique_ptr(hts_open(vcfstr.c_str(), "r"), &hts_close);
        auto header = make_unique_ptr(bcf_hdr_read(file.get()), &bcf_hdr_destroy);
        auto record = make_unique_ptr(bcf_init(), &bcf_destroy);
        bcf_read(file.get(), header.get(), record.get());

        int n_pl = 10*record->n_sample;
        buffer_t<float> buffer_int = hts::bcf::make_buffer<float>(n_pl);
        buffer_t<int> buffer_float = hts::bcf::make_buffer<int>(n_pl);
        // Get PL values
        int n = hts::bcf::get_format_numeric(header.get(), record.get(), "PL", &buffer_int, &n_pl);
        BOOST_REQUIRE_EQUAL(n,expected);
        n = hts::bcf::get_format_numeric(header.get(), record.get(), "PL", &buffer_float, &n_pl);
        BOOST_REQUIRE_EQUAL(n,expected);
    }};

    test(make_data_url(vcftypechar), -2);
    test(make_data_url(vcftypestring), -2);
}

const char vcfinfo[] =
    "##fileformat=VCFv4.2\n"
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n"
    "##contig=<ID=1,length=100>\n"
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"# high-quality bases\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\tDP=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20\t.\t.\n"
;

BOOST_AUTO_TEST_CASE(test_vcf_info_DP) {
    int counter = 0;
    auto test = [&](const std::string &vcfstr, const std::vector<int> &expected) -> void {
    BOOST_TEST_CONTEXT("counter=" << ++counter) {
        // Set header and rec for reading PL values
        auto file = make_unique_ptr(hts_open(vcfstr.c_str(), "r"), &hts_close);
        auto header = make_unique_ptr(bcf_hdr_read(file.get()), &bcf_hdr_destroy);
        auto record = make_unique_ptr(bcf_init(), &bcf_destroy);
        bcf_read(file.get(), header.get(), record.get());

        int n_pl = 8;
        buffer_t<int> buffer = hts::bcf::make_buffer<int>(n_pl);

        // Get PL info
        int n = hts::bcf::get_info_int32(header.get(), record.get(), "DP", &buffer, &n_pl);
        BOOST_REQUIRE_GE(n,0);
        auto test = make_test_range(buffer.get(),buffer.get()+n);
        CHECK_EQUAL_RANGES(test, expected);
    }};

    test(make_data_url(vcfinfo), {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20});
}


const char vcfgt[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
;

BOOST_AUTO_TEST_CASE(test_vcf_gt) {
    int counter = 0;
    auto test = [&](const std::string &vcfstr, int n_gt, const std::vector<int> &expected) -> void {
    BOOST_TEST_CONTEXT("counter=" << ++counter) {
        auto file = make_unique_ptr(hts_open(vcfstr.c_str(), "r"), &hts_close);
        auto header = make_unique_ptr(bcf_hdr_read(file.get()), &bcf_hdr_destroy);
        auto record = make_unique_ptr(bcf_init(), &bcf_destroy);
        bcf_read(file.get(), header.get(), record.get());

        buffer_t<int32_t> buffer = hts::bcf::make_buffer<int32_t>(n_gt);
        int n = hts::bcf::get_genotypes(header.get(), record.get(), &buffer, &n_gt);
        BOOST_REQUIRE_GE(n,0);
        auto test = make_test_range(buffer.get(),buffer.get()+n);
        CHECK_EQUAL_RANGES(test, expected);
    }};

    std::string h = vcfgt;
    test(make_data_url(h+"1\t1\t.\tG\tA\t.\t.\t.\tGT\t0/1\t0|1\t1/.\t1\n"), 8, {
        encode_allele_unphased(0), encode_allele_unphased(1),
        encode_allele_unphased(0), encode_allele_phased(1),
        encode_allele_unphased(1), encode_allele_missing(),
        encode_allele_unphased(1), bcf::int32_vector_end
    });
    auto z = encode_allele_unphased(0);
    test(make_data_url(h+"1\t1\t.\tG\tA\t.\t.\t.\tGT\t0/0/0/0/0/0\t0/0/0/0/0/0\t0/0/0/0/0/0\t0/0/0/0/0/0\n"), 4, {
        z,z,z,z,z,z, z,z,z,z,z,z, z,z,z,z,z,z, z,z,z,z,z,z
    });
}

const char vcftrim[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Posterior Probability\">\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
;

BOOST_AUTO_TEST_CASE(test_variant_trim) {
    int counter = 0;

    auto fm = -1027.1027f;
    auto fe = -1027.1028f;

    auto test = [&](const std::string &vcfstr, double af,
        const std::vector<std::string> &expected_alleles,
        const std::vector<int> &expected_gts,
        const std::vector<float> &expected_gps) -> void {
    BOOST_TEST_CONTEXT("counter=" << ++counter) {
        auto file = bcf::File(vcfstr.c_str(), "r");
        BOOST_REQUIRE(file.is_open());
        auto record = file.InitVariant();
        file.ReadRecord(&record);
        record.Unpack();
        bool ret = record.TrimAlleles(af);
        BOOST_CHECK_EQUAL(ret, true);
        std::vector<std::string> test_alleles;
        for(int i=0;i<record.num_alleles();++i) {
            test_alleles.emplace_back(record.allele(i));
        }
        CHECK_EQUAL_RANGES(test_alleles, expected_alleles);

        // Test GT field
        int n_gt = 4, n_gp = 4, n;
        auto int_buffer = hts::bcf::make_buffer<int32_t>(n_gt);
        n = record.get_genotypes(&int_buffer, &n_gt);
        if(!expected_gts.empty()) {
            BOOST_REQUIRE_GE(n, 0);
            auto test_gts = make_test_range(int_buffer.get(),int_buffer.get()+n);
            CHECK_EQUAL_RANGES(test_gts, expected_gts);            
        } else {
            BOOST_CHECK_LT(n,0);
        }

        // Test GP field
        auto flt_buffer = hts::bcf::make_buffer<float>(n_gt);
        n = record.get_format("GP",&flt_buffer, &n_gp);
        if(!expected_gps.empty()) {
            BOOST_REQUIRE_GE(n, 0);
            for(int i=0;i<n;++i) {
                if(bcf_float_is_missing(flt_buffer[i])) {
                    flt_buffer[i] = fm;
                } else if(bcf_float_is_vector_end(flt_buffer[i])) {
                    flt_buffer[i] = fe;
                }
            }
            auto test_gps = make_test_range(flt_buffer.get(),flt_buffer.get()+n);
            CHECK_EQUAL_RANGES(test_gps, expected_gps);
        } else {
            BOOST_CHECK_LT(n,0);
        }
    }};

    auto z = encode_allele_unphased(0);
    auto x = encode_allele_unphased(1);
    auto e = int32_vector_end;
    auto m = encode_allele_missing();

    std::string h = vcftrim;
    // check if alleles are preserved due to allele frequency
    std::string base1 = make_data_url(h+
        "1\t1\t.\tG\tA\t.\t.\t.\tGT:GP\t0/0:0.9,0.1,0\t0/0:0.9,0.1,0"
    );
    test(base1, 0.0,  {"G", "A"}, {z,z,z,z}, {0.9,0.1,0.0,0.9,0.1,0.0});
    test(base1, 0.05, {"G", "A"}, {z,z,z,z}, {0.9,0.1,0.0,0.9,0.1,0.0});
    test(base1, 0.06,  {"G"}, {z,z,z,z}, {0.9,0.9});

    // check if alleles are preserved due to gt presence
    std::string base2 = make_data_url(h+
        "1\t1\t.\tG\tA\t.\t.\t.\tGT:GP\t0/0:0.9,0.1,0\t0/1:0.9,0.1,0"
    );
    test(base2, 0.0, {"G", "A"}, {z,z,z,x}, {0.9,0.1,0.0,0.9,0.1,0.0});
    test(base2, 0.1, {"G", "A"}, {z,z,z,x}, {0.9,0.1,0.0,0.9,0.1,0.0});

    // check if missing values
    std::string base3 = make_data_url(h+
        "1\t1\t.\tG\tA\t.\t.\t.\tGT:GP\t0/0:0.9,0.1,0.0\t.:."
    );
    test(base3, 0.1, {"G"}, {z,z,m,e}, {0.9f,fm});

    // check if missing GP
    std::string base4 = make_data_url(h+
        "1\t1\t.\tG\tA\t.\t.\t.\tGT\t0/0\t0/0"
    );
    test(base4, 0.1, {"G"}, {z,z,z,z}, {});

    // check if missing GT
    std::string base5 = make_data_url(h+
        "1\t1\t.\tG\tA\t.\t.\t.\tGP\t0.9,0.1,0.0\t."
    );
    test(base5, 0.0, {"G","A"}, {}, {0.9,0.1,0.0,fm,fe,fe});
    test(base5, 0.1, {"G"}, {}, {0.9,fm});

    // check if missing GT and GP
    std::string base6 = make_data_url(h+
        "1\t1\t.\tG\tA\t.\t.\t.\t.\t\t"
    );
    test(base6, 0.0, {"G"}, {}, {});
}
