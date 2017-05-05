/*
 * Copyright (c) 2015 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
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
#ifndef DNG_PEELING_H
#define DNG_PEELING_H

#include <vector>
#include <utility>

#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/fill.hpp>
#include <boost/range/algorithm/fill_n.hpp>

#include <dng/matrix.h>
#include <dng/genotyper.h>
#include <dng/detail/unit_test.h>

namespace dng {
namespace peel {

namespace op {
enum Type {
    UP, DOWN, TOFATHER, TOMOTHER, TOCHILD,
    UPFAST, DOWNFAST, TOFATHERFAST, TOMOTHERFAST,
    TOCHILDFAST,
    NUM // Total number of possible forward operations
};
} // namespace dng::peel::op

struct workspace_t {
    GenotypeArrayVector upper; // Holds P(~Descendent_Data & G=g)
    GenotypeArrayVector lower; // Holds P( Descendent_Data | G=g)
    ParentArrayVector super; // Holds P(~Descendent_Data & G=g) for parent nodes

    bool dirty_lower = false;
    double forward_result;
    // Temporary data used by some peeling ops
    TemporaryMatrix temp_buffer;

    // Information about the pedigree
    std::size_t num_nodes = 0;
    typedef std::pair<std::size_t, std::size_t> node_range_t;
    node_range_t founder_nodes,
                 germline_nodes,
                 somatic_nodes,
                 library_nodes;

    std::vector<int> ploidies;

    // Resize the workspace to fit a pedigree with sz nodes
    void Resize(std::size_t sz) {
        num_nodes = sz;
        upper.resize(num_nodes);
        super.resize(num_nodes);
        lower.assign(num_nodes, DNG_INDIVIDUAL_BUFFER_ONES);
        temp_buffer.resize(100, 1);
        dirty_lower = false;
    }

    // Cleanup after a backwards peeling algorithm
    void Cleanup() {
        boost::fill(lower, DNG_INDIVIDUAL_BUFFER_ONES);
        dirty_lower = false;
    }

    // Cleanup after a backwards peeling algorithm
    // Do not update libraries since they might be set in another operation
    void CleanupFast() {
        // TODO: create a check that sees if this has been done before the
        // forward algorithm.
        boost::fill_n(lower, somatic_nodes.second, DNG_INDIVIDUAL_BUFFER_ONES);
        dirty_lower = false;
    }

    // Set the prior probability of the founders given the reference
    void SetFounders(const GenotypeArray &prior) {
        assert(founder_nodes.first <= founder_nodes.second);
        std::fill(upper.begin() + founder_nodes.first,
                  upper.begin() + founder_nodes.second, prior);
    }

    void SetFounders(const GenotypeArray &diploid_prior, const GenotypeArray &haploid_prior) {
        for(auto i = founder_nodes.first; i < founder_nodes.second; ++i) {
            if(ploidies[i] == 2) {
                upper[i] = diploid_prior;
            } else if(ploidies[i] == 1) {
                upper[i] = haploid_prior;                
            } else {
                assert(false); // should not be here
            }
        } 
    }

    inline
    double SetGenotypeLikelihoods(const Genotyper &gt,
        const pileup::RawDepths &depths, const int ref_index) {
        double scale = 0.0, stemp;
        for(std::size_t u = 0; u < depths.size(); ++u) {
            auto pos = library_nodes.first + u;
            std::tie(lower[pos], stemp) =
                gt(depths, u, ref_index, ploidies[pos]);
            scale += stemp;
        }
        return scale;
    }

    inline
    double SetGenotypeLikelihoods(const Genotyper &gt,
        const pileup::AlleleDepths &depths) {
        double scale=0.0, stemp;
        for(std::size_t u = 0; u < depths.num_libraries(); ++u) {
            auto pos = library_nodes.first + u;
            assert(pos < lower.size());
            std::tie(lower[pos], stemp) =
                gt(depths, u, ploidies[pos]);
            scale += stemp;
        }
        return scale;
}


};

typedef std::vector<std::size_t> family_members_t;

// Basic peeling operations
void up(workspace_t &work, const family_members_t &family,
        const TransitionMatrixVector &mat);
void down(workspace_t &work, const family_members_t &family,
          const TransitionMatrixVector &mat);
void to_father(workspace_t &work, const family_members_t &family,
               const TransitionMatrixVector &mat);
void to_mother(workspace_t &work, const family_members_t &family,
               const TransitionMatrixVector &mat);
void to_child(workspace_t &work, const family_members_t &family,
              const TransitionMatrixVector &mat);

// Fast versions that can be used on "dirty" workspaces
void up_fast(workspace_t &work, const family_members_t &family,
             const TransitionMatrixVector &mat);
void down_fast(workspace_t &work, const family_members_t &family,
               const TransitionMatrixVector &mat);
void to_father_fast(workspace_t &work, const family_members_t &family,
                    const TransitionMatrixVector &mat);
void to_mother_fast(workspace_t &work, const family_members_t &family,
                    const TransitionMatrixVector &mat);
void to_child_fast(workspace_t &work, const family_members_t &family,
                   const TransitionMatrixVector &mat);

// Reverse versions that can be used for backwards algorithms
void up_reverse(workspace_t &work, const family_members_t &family,
                const TransitionMatrixVector &mat);
void down_reverse(workspace_t &work, const family_members_t &family,
                  const TransitionMatrixVector &mat);
void to_father_reverse(workspace_t &work, const family_members_t &family,
                       const TransitionMatrixVector &mat);
void to_mother_reverse(workspace_t &work, const family_members_t &family,
                       const TransitionMatrixVector &mat);
void to_child_reverse(workspace_t &work, const family_members_t &family,
                      const TransitionMatrixVector &mat);

typedef decltype(&down) function_t;

struct info_t {
    bool writes_lower;
    int writes_to;
};

constexpr info_t info[op::NUM] = {
    /* Up           */ {true,  0},
    /* Down         */ {false, 1},
    /* ToFather     */ {true,  0},
    /* ToMother     */ {true,  1},
    /* ToChild      */ {false, 2},
    /* UpFast       */ {true,  0},
    /* DownFast     */ {false, 1},
    /* ToFatherFast */ {true,  0},
    /* ToMotherFast */ {true,  1},
    /* ToChildFast  */ {false, 2}
};

// TODO: Write test case to check that peeling ops are in the right order.
constexpr function_t functions[op::NUM] = {
    &up, &down, &to_father, &to_mother, &to_child,
    &up_fast, &down_fast, &to_father_fast, &to_mother_fast,
    &to_child_fast
};

constexpr function_t reverse_functions[op::NUM] = {
    &up_reverse, &down_reverse, &to_father_reverse,
    &to_mother_reverse, &to_child_reverse,
    &up_reverse, &down_reverse, &to_father_reverse,
    &to_mother_reverse, &to_child_reverse
};

} // namespace dng::peel
} // namespace dng

#endif // DNG_PEELING_H
