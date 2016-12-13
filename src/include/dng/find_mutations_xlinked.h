/*
 * Copyright (c) 2016 Steven H. Wu
 * Copyright (c) 2016 Reed A. Cartwright
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

#pragma once
#ifndef DNG_FIND_MUTATIONS_XLINKED_H_
#define DNG_FIND_MUTATIONS_XLINKED_H_

#include <dng/find_mutations_abstract.h>
#include <dng/detail/unit_test.h>

namespace dng {

class FindMutationsXLinked : public FindMutationsAbstract {

public:

    FindMutationsXLinked(double min_prob, const RelationshipGraph &graph,
            params_t params);

    ~FindMutationsXLinked();

    bool operator()(const std::vector<depth_t> &depths, int ref_index,
                    stats_t *stats);


protected:
    void SetupTransitionMatrix();


    DNG_UNIT_TEST(test_constructor);
    DNG_UNIT_TEST(test_prior);
    DNG_UNIT_TEST(test_genotype);
    DNG_UNIT_TEST(test_full_transition);
    DNG_UNIT_TEST(test_operator);
    DNG_UNIT_TEST(test_operator_m12);
};
} // namespace dng

#endif /* DNG_FIND_MUTATIONS_XLINKED_H_ */
