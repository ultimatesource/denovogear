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
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"

const char vcfformatint[] =
    "##fileformat=VCFv4.2\n"
    "##INFO=<ID=I16,Number=.,Type=Integer,Description=\"\">\n"
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n"
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"# high-quality bases\">\n"
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n"
    "##contig=<ID=2,length=243199373>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878_vald-sorted.bam.bam\tNA12891_vald-sorted.bam.bam\tNA12892_vald-sorted.bam.bam\n"
    "2\t214668360\t.\tG\tA,N\t0\t.\tDP=158;I16=5,126,1,26,3276,83192,658,16086,7643,452287,1620,97200,2352,43808,491,8985\tPL:DP\t114,0,90,177,172,215:48\t0,229,233,229,233,233:76\t0,102,142,102,142,142:34\n"
;

const char vcfformatfloat[] =
    "##fileformat=VCFv4.2\n"
    "##INFO=<ID=I16,Number=.,Type=Integer,Description=\"\">\n"
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n"
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"# high-quality bases\">\n"
    "##FORMAT=<ID=PL,Number=G,Type=Float,Description=\"List of Phred-scaled genotype likelihoods\">\n"
    "##contig=<ID=2,length=243199373>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878_vald-sorted.bam.bam\tNA12891_vald-sorted.bam.bam\tNA12892_vald-sorted.bam.bam\n"
    "2\t214668360\t.\tG\tA,N\t0\t.\tDP=158;I16=5,126,1,26,3276,83192,658,16086,7643,452287,1620,97200,2352,43808,491,8985\tPL:DP\t114.0,0.0,90.0,177,172,215:48\t0.0,229.0,233.0,229,233,233:76\t0.0,102.0,142.0,102,142,142:34\n"
;

using namespace hts;
using namespace hts::bcf;
using hts::detail::make_data_url;

BOOST_AUTO_TEST_CASE(test_vcf_read_intPL) {
    auto test = [](std::string vcfstr, int expected_result[]) -> void {
	BOOST_TEST_CONTEXT("Reading int PL from VCF\n") {
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
	    buffer_t<int> buffer = hts::bcf::make_buffer<int>(n_pl);

	    // Get PL values
	    int n = hts::bcf::get_format_numeric(hdr, rec, "PL", &buffer, &n_pl);
	    BOOST_CHECK(n > 0);
	    for(int i=0;i<n;i++) {
		BOOST_CHECK_EQUAL(buffer[i],expected_result[i]);
	    }
    }};

    std::string vcfstr = make_data_url(vcfformatint);
    int p[] = {114,0,90,177,172,215,0,229,233,229,233,233,0,102,142,102,142,142};
    test(vcfstr,p);
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
