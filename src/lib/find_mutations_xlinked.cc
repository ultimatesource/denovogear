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


#include <dng/find_mutations_xlinked.h>


using namespace dng;

FindMutationsXLinked::~FindMutationsXLinked() {
	// TODO Auto-generated destructor stub
}

FindMutationsXLinked::FindMutationsXLinked(double min_prob,
        const RelationshipGraph &graph, params_t params)
        : FindMutationsAbstract(min_prob, graph, params) {

    std::cerr << "This is class is not implemented yet!!" << std::endl;
    SetupPopulationPriorHaploid();
    SetupTransitionMatrix();


#if CALCULATE_ENTROPY == 1
    std::cerr << "Entropy will not be calculated for X-linked model!!" << std::endl;
#endif

}

// Returns true if a mutation was found and the record was modified
bool FindMutationsXLinked::operator()(const std::vector<depth_t> &depths,
                               int ref_index, stats_t *stats) {
    //PR_NOTE(SW): Implement this at next PR
    return false;
}



void FindMutationsXLinked::SetupTransitionMatrix(){
    //PR_NOTE(SW): Implement this at next PR
}
