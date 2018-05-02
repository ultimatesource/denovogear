/*
 * Author: Juan J. Garcia Mesa <jgarc111@asu.edu>
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

const char vcfformatint[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=Integer>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1,2,3\n"
;

const char vcfformatfloat[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=Float>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1.1,2.1,3.1\n";

using namespace hts;
using namespace hts::bcf;
using hts::detail::make_data_url;
using dng::detail::make_test_range;

BOOST_AUTO_TEST_CASE(test_vcf_read_int_PL) {
    int counter = 0;
    auto test = [&](const std::string &vcfstr, const std::vector<int> &expected) -> void {
    BOOST_TEST_CONTEXT("counter=" << ++counter) {
        // Set header and rec for reading PL values
        auto file = hts_open(vcfstr.c_str(), "r");
        auto header = bcf_hdr_read(file);
        auto record = bcf_init();
        bcf_read(file, header, record);
        bcf_unpack(record, BCF_UN_STR);

        int n_pl = 10*record->n_sample;
        buffer_t<int> buffer = hts::bcf::make_buffer<int>(n_pl);

        // Get PL values
        int n = hts::bcf::get_format_numeric(header, record, "PL", &buffer, &n_pl);
        BOOST_REQUIRE_GE(n,0);
        auto test = make_test_range(buffer.get(),buffer.get()+n);
        CHECK_EQUAL_RANGES(test, expected);
    }};

    test(make_data_url(vcfformatint), {1,2,3});
    test(make_data_url(vcfformatfloat), {1,2,3});
}

BOOST_AUTO_TEST_CASE(test_vcf_read_floatPL) {
    auto test = [](std::string vcfstr, float expected_result[]) -> void {
    BOOST_TEST_CONTEXT("Reading float PL from VCF\n") {
        // Set header and rec for reading PL values
        bcf_srs_t *rec_reader = bcf_sr_init();
        int ret = bcf_sr_add_reader(rec_reader, vcfstr.c_str());
        BOOST_CHECK(ret != 0);
        const bcf_hdr_t *hdr = bcf_sr_get_header(rec_reader, 0);
        bcf_sr_next_line(rec_reader);
        bcf1_t *rec = bcf_sr_get_line(rec_reader, 0);
        bcf_unpack(rec,BCF_UN_STR);
        uint32_t n_alleles = rec->n_allele;
        int n_pl = n_alleles*10;
        buffer_t<float> buffer = hts::bcf::make_buffer<float>(n_pl);

        // Get PL values
        int n = hts::bcf::get_format_numeric(hdr, rec, "PL", &buffer, &n_pl);
        BOOST_CHECK(n > 0);
        for(int i=0;i<n;i++) {
        BOOST_CHECK_EQUAL(buffer[i],expected_result[i]);
        }

    }};

    std::string vcfstrfloat = make_data_url(vcfformatfloat);
    std::string vcfstrint = make_data_url(vcfformatint);
    float p[] = {114.0,0.0,90.0,177.0,172.0,215.0,0.0,229.0,233.0,229.0,233.0,233.0,0.0,102.0,142.0,102.0,142.0,142.0};
    test(vcfstrfloat,p);
    test(vcfstrint,p);
}
