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
#ifndef DENOVOGEAR_FIXTURE_TRIO_WORKSPACE_H
#define DENOVOGEAR_FIXTURE_TRIO_WORKSPACE_H

#include <algorithm>

#include <dng/utility.h>
#include <dng/find_mutations.h>
#include "fixture_read_trio_from_file.h"

struct TrioWorkspace : public  ReadTrioFromFile {

    double min_prob;

    dng::Pedigree pedigree;
    dng::peel::workspace_t workspace;

    int ref_index = 2;
    std::vector<depth_t> read_depths{3};

    FindMutations::params_t test_param_1 {0, {{0,0,0,0}}, 0,
                                          std::string{"0,0,0,0"},
                                          std::string{"0,0,0,0"} };

    TrioWorkspace(std::string s = "TrioWorkspace") : ReadTrioFromFile(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);


        pedigree.Construct(ped, rgs, arg.mu, arg.mu_somatic, arg.mu_library);

        std::array<double, 4> freqs;
        auto f = dng::utility::parse_double_list(arg.nuc_freqs, ',', 4);
        std::copy(f.first.begin(), f.first.end(), &freqs[0]);

        test_param_1 = FindMutations::params_t {arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1]};


        int min_qual = arg.min_basequal;
        min_prob = arg.min_prob;

        ref_index = 2;
        uint16_t cc[3][4] = {{0, 1, 25, 29},
                             {0, 0, 57, 0},
                             {0, 0, 76, 1}};
        for (int j = 0; j < 3; ++j) {
            std::copy(cc[j], cc[j] + 4, read_depths[j].counts);
        }
        setup_workspace(ref_index, read_depths);

    }

    double setup_workspace(int ref_index, std::vector<depth_t> &read_depths ){

        workspace.Resize(5);
        workspace.founder_nodes = std::make_pair(0, 2);
        workspace.germline_nodes = std::make_pair(0, 2);
        workspace.somatic_nodes = std::make_pair(2, 2);
        workspace.library_nodes = std::make_pair(2, 5);

        std::array<double, 4> prior {};
        prior.fill(0);
        prior[ref_index] = test_param_1.ref_weight;
        auto genotype_prior_prior = population_prior(test_param_1.theta, test_param_1.nuc_freq, prior);
        workspace.SetFounders(genotype_prior_prior);

        std::vector<std::string> expect_gamma{"0.98, 0.0005, 0.0005, 1.04",
                                              "0.02, 0.075,  0.005,  1.18"};
        dng::genotype::DirichletMultinomialMixture genotype_likelihood_{
                dng::genotype::DirichletMultinomialMixture::params_t {expect_gamma[0]},
                dng::genotype::DirichletMultinomialMixture::params_t {expect_gamma[1]}  };
        double scale = 0.0, stemp;
        for (std::size_t u = 0; u < read_depths.size(); ++u) {
            std::tie(workspace.lower[2 + u], stemp) =
                    genotype_likelihood_(read_depths[u], ref_index);
            scale += stemp;
        }
        return scale;
    }

    ~TrioWorkspace() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }


};


#endif //DENOVOGEAR_FIXTURE_TRIO_WORKSPACE_H
