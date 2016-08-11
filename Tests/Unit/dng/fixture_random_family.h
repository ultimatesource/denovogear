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
#ifndef DENOVOGEAR_FIXTURE_RANDOM_FAMILY_H
#define DENOVOGEAR_FIXTURE_RANDOM_FAMILY_H

#include <dng/peeling.h>


struct RandomFamily {
    const int CHILD_OFFSET = 2;
    std::string fixture;

    std::mt19937 random_gen_mt {1}; //seed = 1

    dng::peel::family_members_t family;
    std::vector<dng::TransitionMatrix> trans_matrix;
    std::vector<dng::GenotypeArray> upper_array;
    std::vector<dng::GenotypeArray> lower_array;

    int num_child;
    int total_family_size;
    std::uniform_int_distribution<> rand_unif;

    RandomFamily(std::string s = "RandomFamily") : fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture:  " << fixture);
        rand_unif = std::uniform_int_distribution<>(1,10);
        init_family();
    }

    void init_family(){

        int num_child = rand_unif(random_gen_mt);
        total_family_size = num_child + CHILD_OFFSET;
        family.clear();
        trans_matrix.resize(total_family_size);
        upper_array.resize(total_family_size);
        lower_array.resize(total_family_size);

        for (int k = 0; k < total_family_size; ++k) {
            family.push_back(k);
            lower_array[k] = dng::GenotypeArray::Random(10);
            if (k < CHILD_OFFSET) {
                trans_matrix[k] = dng::TransitionMatrix::Random(10, 10);
                upper_array[k] = dng::GenotypeArray::Random(10);
            }
            else {
                trans_matrix[k] = dng::TransitionMatrix::Random(100, 10);
            }

        }
    }

    void init_family_parent_child_only(){

        total_family_size = 2;

        family.clear();
        trans_matrix.resize(total_family_size);
        upper_array.resize(total_family_size);
        lower_array.resize(total_family_size);

        for (int k = 0; k < total_family_size; ++k) {
            family.push_back(k);
            trans_matrix[k] = dng::TransitionMatrix::Random(10, 10);
            lower_array[k] = dng::GenotypeArray::Random(10);
//            upper_array[k] = GenotypeArray::Random();

        }
    }

    ~RandomFamily() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }
};

void copy_family_to_workspace(dng::peel::workspace_t &workspace,
		dng::TransitionVector &full_matrix, int total_family_size,
		const std::vector<dng::GenotypeArray> &lower,
		const std::vector<dng::GenotypeArray> &upper,
		const std::vector<dng::TransitionMatrix> &trans_mat) {
	workspace.Resize(total_family_size);
    full_matrix.resize(total_family_size);
    for (int l = 0; l < total_family_size; ++l) {
        workspace.lower[l] = lower[l];
        workspace.upper[l] = upper[l];
        full_matrix[l] = trans_mat[l];
    }
}


struct RandomWorkspace {

	std::string fixture;

	dng::peel::workspace_t workspace;
	dng::TransitionVector full_matrix;

	RandomWorkspace(std::string s = "RandomWorkspace") : fixture(s) {
		BOOST_TEST_MESSAGE("set up fixture:  " << s);
	}

	void init_workspace() {
		RandomFamily family;

		family.init_family();
		copy_family_to_workspace(workspace, full_matrix,
				family.total_family_size, family.lower_array,
				family.upper_array, family.trans_matrix);

	}

    ~RandomWorkspace() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }

};



#endif //DENOVOGEAR_FIXTURE_RANDOM_FAMILY_H
