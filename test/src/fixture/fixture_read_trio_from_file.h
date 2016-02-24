//
// Created by steven on 2/12/16.
//

#pragma once
#ifndef DENOVOGEAR_FIXTURE_READ_TRIO_FROM_FILE_H
#define DENOVOGEAR_FIXTURE_READ_TRIO_FROM_FILE_H

#include <fstream>

#include <dng/task/call.h>
#include <utils/boost_utils.h>
#include <helpers/find_mutation_helper.h>


using namespace dng;
namespace utf = boost::unit_test;

struct ReadTrioFromFile {

    std::string fixture;

    dng::io::Pedigree ped;
    dng::ReadGroups rgs;

    typedef dng::task::Call task_type;
    struct arg_t : public task_type::argument_type {
        bool help;
        bool version;
        std::string arg_file;

        std::string run_name;
        std::string run_path;
    } arg;

    ReadTrioFromFile(std::string s = "ReadTrioFromFile") : fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);

        po::options_description ext_desc_, int_desc_;
        po::positional_options_description pos_desc_;
        po::variables_map vm_;

        int argc=4;
        char *argv[argc];
        argv[0] = (char*) "test";
        argv[1] = (char*) "-p";

        std::string ped_filename (TESTDATA_DIR);
        ped_filename.append("/sample_5_3/ceu.ped");
        argv[2] = (char*) ped_filename.data();

        std::string vcf_filename = TESTDATA_DIR;
        vcf_filename.append("/sample_5_3/test1.vcf");
        argv[3] = (char*) vcf_filename.data();

        add_app_args(ext_desc_, static_cast<typename task_type::argument_type &>(arg));
        int_desc_.add_options()
                ("input", po::value< std::vector<std::string> >(&arg.input), "input files")
                ;
        int_desc_.add(ext_desc_);
        pos_desc_.add("input", -1);
        po::store(po::command_line_parser(argc, argv)
                          .options(int_desc_).positional(pos_desc_).run(), vm_);
        po::notify(vm_);

        // Parse pedigree from file

        std::ifstream ped_file(arg.ped);
        ped.Parse(utils::istreambuf_range(ped_file));


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



    ~ReadTrioFromFile() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }


};


#endif //DENOVOGEAR_FIXTURE_READ_TRIO_FROM_FILE_H
