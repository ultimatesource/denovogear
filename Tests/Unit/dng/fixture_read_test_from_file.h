/*
 * Copyright (c) 2016 Steven H. Wu
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
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

#pragma once
#ifndef DNG_FIXTURE_READ_TEST_FROM_FILE_H
#define DNG_FIXTURE_READ_TEST_FROM_FILE_H

#include <fstream>

#include <dng/task/call.h>
#include <dng/hts/bcf.h>
#include <dng/vcfpileup.h>
#include <dng/seq.h>
#include <dng/io/utility.h>
#include <dng/find_mutations_abstract.h>

using namespace dng;

//// Helper function that mimics boost::istream_range
//template<class Elem, class Traits> inline boost::iterator_range<
//        std::istreambuf_iterator<Elem, Traits>> istreambuf_range(
//        std::basic_istream<Elem, Traits> &in) {
//    return boost::iterator_range<std::istreambuf_iterator<Elem, Traits>>(
//            std::istreambuf_iterator<Elem, Traits>(in),
//            std::istreambuf_iterator<Elem, Traits>());
//}

struct ReadFromFile {
    std::string fixture;

    dng::io::Pedigree io_pedigree;
    dng::ReadGroups rgs;

    double min_prob;
    int ref_index;
    std::vector<depth_t> read_depths;
    FindMutationsAbstract::params_t test_param_1 {0, { {0, 0, 0, 0}}, 0,
        std::string {"0,0,0,0"}, std::string {"0,0,0,0"}};


    typedef dng::task::Call task_type;
    struct arg_t : public task_type::argument_type {
        bool help;
        bool version;
        std::string arg_file;

        std::string run_name;
        std::string run_path;
    } arg;

    ReadFromFile(std::string s = "ReadFromFile") : fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture: ReadFromFile:" << fixture);

    }

    ~ReadFromFile() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }

    void ReadFile(std::string &ped_filename, std::string &vcf_filename){


        po::options_description ext_desc, int_desc;
        po::positional_options_description pos_desc;
        po::variables_map vm;

        int argc=4;
        char *argv[argc];
        argv[0] = (char*) "test";
        argv[1] = (char*) "-p";
        argv[2] = (char*) ped_filename.data();
        argv[3] = (char*) vcf_filename.data();

        add_app_args(ext_desc, static_cast<typename task_type::argument_type &>(arg));
        int_desc.add_options()
                ("input", po::value< std::vector<std::string> >(&arg.input), "input files")
                ;
        int_desc.add(ext_desc);
        pos_desc.add("input", -1);
        po::store(po::command_line_parser(argc, argv)
                          .options(int_desc).positional(pos_desc).run(), vm);
        po::notify(vm);

        // Parse pedigree from file

        std::ifstream ped_file(arg.ped);
        io_pedigree.Parse(io::istreambuf_range(ped_file));


        std::vector<hts::File> indata;
        std::vector<hts::bcf::File> bcfdata;
        for (auto &&str : arg.input) {
            indata.emplace_back(str.c_str(), "r");
            if (indata.back().is_open()) {
                continue;
            }
            throw std::runtime_error("unable to open input file '" + str + "'.");
        }
        bcfdata.emplace_back(std::move(indata[0]));
        rgs.ParseSamples(bcfdata[0]);

    }

    void InitFromDefaultArg(){
        min_prob = arg.min_prob;

        std::array<double, 4> freqs;
        auto f = dng::utility::parse_double_list(arg.nuc_freqs, ',', 4);
        std::copy(f.first.begin(), f.first.end(), &freqs[0]);
        test_param_1 = FindMutationsAbstract::params_t {arg.theta, freqs,
                arg.ref_weight, arg.gamma[0], arg.gamma[1]};

    }

};




struct ReadTrioFromFile : public ReadFromFile{
    std::string fixture;

    ReadTrioFromFile(std::string s = "ReadTrioFromFile")
            : ReadFromFile(" from "+s), fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture: ReadTrioFromFile: " << fixture);


        std::string ped_filename (TESTDATA_DIR);
        ped_filename.append("/sample_5_3/ceu.ped");
        std::string vcf_filename = TESTDATA_DIR;
        vcf_filename.append("/sample_5_3/test1.vcf");

        ReadFile(ped_filename, vcf_filename);
        InitFromDefaultArg();

        ref_index = 2;
        read_depths.resize(rgs.libraries().size());
        uint16_t cc[3][4] = {{0, 1, 25, 29},
                             {0, 0, 57, 0},
                             {0, 0, 76, 1}};
        for (int j = 0; j < 3; ++j) {
            std::copy(cc[j], cc[j] + 4, read_depths[j].counts);
        }
    }

    ~ReadTrioFromFile() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }
};


struct ReadM12FromFile : public ReadFromFile {
    std::string fixture;

    ReadM12FromFile(std::string s = "ReadM12FromFile") : ReadFromFile(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);

        std::string ped_filename (TESTDATA_DIR);
        ped_filename.append("/relationship_graph/relationship_graph.ped");

        std::string vcf_filename = TESTDATA_DIR;
        vcf_filename.append("/relationship_graph/relationship_graph.vcf");

        ReadFile(ped_filename, vcf_filename);
        InitFromDefaultArg();

        ref_index = 2;
        read_depths.resize(rgs.libraries().size());
        uint16_t cc[3][4] = {{0, 1, 25, 29},
                             {0, 0, 57, 0},
                             {0, 0, 76, 1}};
        for (int i = 0; i < 12;) {
            for (int j = 0; j < 3; ++j) {
                std::copy(cc[j], cc[j] + 4, read_depths[i++].counts);
            }
        }
    }

    ~ReadM12FromFile() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }
};


#endif //DNG_FIXTURE_READ_TEST_FROM_FILE_H
