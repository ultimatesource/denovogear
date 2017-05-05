/*
 * Copyright (c) 2014-2015 Reed A. Cartwright
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

#include <dng/peeling.h>

#define DEBUG_PEELING 0
#if DEBUG_PEELING == 1
#include <iostream>
#endif
// Family Order: Parent, Child
void dng::peel::down(workspace_t &work, const family_members_t &family,
                     const TransitionMatrixVector &mat) {
    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];
    work.upper[child] = (mat[child].transpose() * (work.upper[parent] *
                         work.lower[parent]).matrix()).array();
}

// Family Order: Parent, Child
void dng::peel::down_fast(workspace_t &work, const family_members_t &family,
                          const TransitionMatrixVector &mat) {
    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];
    work.upper[child] = (mat[child].transpose() *
                         work.upper[parent].matrix()).array();
}

// Family Order: Parent, Child
void dng::peel::up(workspace_t &work, const family_members_t &family,
                   const TransitionMatrixVector &mat) {
    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];
    work.lower[parent] *= (mat[child] * work.lower[child].matrix()).array();
}

// Family Order: Parent, Child
void dng::peel::up_fast(workspace_t &work, const family_members_t &family,
                        const TransitionMatrixVector &mat) {
    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];
    work.lower[parent] = (mat[child] * work.lower[child].matrix()).array();
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_father(workspace_t &work, const family_members_t &family,
                          const TransitionMatrixVector &mat) {
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.temp_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.temp_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    // Include Mom
    auto mom_width = work.upper[mom].size();
    auto dad_width = work.temp_buffer.size()/mom_width;
    work.temp_buffer.resize(mom_width, dad_width);
    work.lower[dad] *= (work.temp_buffer.matrix().transpose() * (work.upper[mom] *
                        work.lower[mom]).matrix()).array();

    //work.temp_buffer.resize(work.temp_buffer.size(), 1);
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_father_fast(workspace_t &work,
                               const family_members_t &family, const TransitionMatrixVector &mat) {
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.temp_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.temp_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    // Include Mom
    auto mom_width = work.upper[mom].size();
    auto dad_width = work.temp_buffer.size()/mom_width;
    work.temp_buffer.resize(mom_width, dad_width);
    work.lower[dad] = (work.temp_buffer.matrix().transpose() * (work.upper[mom] *
                       work.lower[mom]).matrix()).array();
    //work.temp_buffer.resize(work.temp_buffer.size(), 1);
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_mother(workspace_t &work, const family_members_t &family,
                          const TransitionMatrixVector &mat) {
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.temp_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(size_t i = 3; i < family.size(); ++i) {
        work.temp_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    // Include Mom
    auto dad_width = work.upper[dad].size();
    auto mom_width = work.temp_buffer.size()/dad_width;
    work.temp_buffer.resize(mom_width, dad_width);
    work.lower[mom] *= (work.temp_buffer.matrix() * (work.upper[dad] *
                        work.lower[dad]).matrix()).array();
    //work.temp_buffer.resize(work.temp_buffer.size(), 1);
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_mother_fast(workspace_t &work,
                               const family_members_t &family, const TransitionMatrixVector &mat) {
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.temp_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.temp_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    // Include Dad
    auto dad_width = work.upper[dad].size();
    auto mom_width = work.temp_buffer.size()/dad_width;
    work.temp_buffer.resize(mom_width, dad_width);
    work.lower[mom] = (work.temp_buffer.matrix() * (work.upper[dad] *
                       work.lower[dad]).matrix()).array();
    //work.temp_buffer.resize(work.temp_buffer.size(), 1);
}

// Family Order: Father, Mother, Child, Child2, ....
void dng::peel::to_child(workspace_t &work, const family_members_t &family,
                         const TransitionMatrixVector &mat) {
    assert(family.size() >= 4);
    auto dad = family[0];
    auto mom = family[1];
    auto child = family[2];
    // Parents
    work.temp_buffer = kroneckerProduct((work.lower[dad] *
                                           work.upper[dad]).matrix(),
                                          (work.lower[mom] * work.upper[mom]).matrix()).array();
    // Sum over fullsibs
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.temp_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }

    work.upper[child] = (mat[child].transpose() *
                         work.temp_buffer.matrix()).array();
}

// Family Order: Father, Mother, CHild
void dng::peel::to_child_fast(workspace_t &work, const family_members_t &family,
                              const TransitionMatrixVector &mat) {
    assert(family.size() == 3);
    auto dad = family[0];
    auto mom = family[1];
    auto child = family[2];
    work.upper[child] = (mat[child].transpose() * kroneckerProduct(
                             (work.lower[dad] * work.upper[dad]).matrix(),
                             (work.lower[mom] * work.upper[mom]).matrix())).array();
}

// Family Order: Parent, Child
void dng::peel::down_reverse(workspace_t &work, const family_members_t &family,
                             const TransitionMatrixVector &mat) {
    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];

    work.super[child] = work.upper[parent];
    work.lower[parent] *= (mat[child] * work.lower[child].matrix()).array();
}

// Family Order: Parent, Child
// prevent divide by zero errors by adding a minor offset to one of the calculations
void dng::peel::up_reverse(workspace_t &work, const family_members_t &family,
                           const TransitionMatrixVector &mat) {
    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];

    work.super[child] = work.upper[parent] * (work.lower[parent] /
                        ((mat[child] * work.lower[child].matrix()).array() +
                         DNG_INDIVIDUAL_BUFFER_MIN));
    work.upper[child] = (mat[child].transpose() *
                         work.super[child].matrix()).array();
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_father_reverse(workspace_t &work,
                                  const family_members_t &family, const TransitionMatrixVector &mat) {
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    // work.temp_buffer will contain P(child data | mom & dad)
    work.temp_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.temp_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }

    // Calculate P(mom-only data & mom = g)
    GenotypeArrayVector::value_type mom_v = work.upper[mom] * work.lower[mom];

    // Calculate P(dad-only data & dad = g)
    auto mom_width = work.upper[mom].size();
    auto dad_width = work.upper[dad].size();
    work.temp_buffer.resize(mom_width, dad_width);
    GenotypeArrayVector::value_type dad_v = work.upper[dad] * (work.lower[dad] /
                                         ((work.temp_buffer.matrix().transpose() * mom_v.matrix()).array() +
                                          DNG_INDIVIDUAL_BUFFER_MIN));

    // Calculate P(dependent data | mom = g)
    work.lower[mom] *= (work.temp_buffer.matrix() *
                        dad_v.matrix()).array();
    work.temp_buffer.resize(work.temp_buffer.size(), 1);

    // Calculate P(data & dad & mom)
    work.temp_buffer *= kroneckerProduct(dad_v.matrix(), mom_v.matrix()).array();

    for(std::size_t i = 2; i < family.size(); ++i) {
        auto child = family[i];
        work.super[child] = work.temp_buffer / ((mat[child] *
                            work.lower[child].matrix()).array() + DNG_INDIVIDUAL_BUFFER_MIN);
        work.upper[child] = mat[child].transpose() * work.super[child].matrix();
    }
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_mother_reverse(workspace_t &work,
                                  const family_members_t &family, const TransitionMatrixVector &mat) {
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.temp_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.temp_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    GenotypeArrayVector::value_type dad_v = work.upper[dad] * work.lower[dad];

    auto mom_width = work.upper[mom].size();
    auto dad_width = work.upper[dad].size();
    work.temp_buffer.resize(mom_width, dad_width); //TODO: FIXME: add transpose()
    GenotypeArrayVector::value_type mom_v = work.upper[mom] * (work.lower[mom] /
                                         ((work.temp_buffer.matrix().transpose() * dad_v.matrix()).array() +
                                          DNG_INDIVIDUAL_BUFFER_MIN));
    work.lower[dad] *= (work.temp_buffer.matrix() * mom_v.matrix()).array(); //TODO: FIXME: add transpose()
    work.temp_buffer.resize(work.temp_buffer.size(), 1);

    work.temp_buffer *= kroneckerProduct(dad_v.matrix(), mom_v.matrix()).array();

    for(std::size_t i = 2; i < family.size(); ++i) {
        auto child = family[i];
        work.super[child] = work.temp_buffer / ((mat[child] *
                            work.lower[child].matrix()).array() + DNG_INDIVIDUAL_BUFFER_MIN);
        work.upper[child] = mat[child].transpose() * work.super[child].matrix();
    }
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_child_reverse(workspace_t &work,
                                 const family_members_t &family, const TransitionMatrixVector &mat) {
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    auto child = family[2];
    // Sum over children
    work.temp_buffer = (mat[child] * work.lower[child].matrix()).array();
    for(std::size_t i = 3; i < family.size(); i++) {
        work.temp_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    // Update Parents
    GenotypeArrayVector::value_type dad_v = work.upper[dad] * work.lower[dad];
    GenotypeArrayVector::value_type mom_v = work.upper[mom] * work.lower[mom];

    auto mom_width = work.upper[mom].size();
    auto dad_width = work.upper[dad].size();
    work.temp_buffer.resize(mom_width, dad_width); //TODO: FIXME: add transpose()
    work.lower[dad] *= (work.temp_buffer.matrix() * mom_v.matrix()).array(); //TODO: FIXME: add transpose()
    work.lower[mom] *= (work.temp_buffer.matrix().transpose() *
                        dad_v.matrix()).array(); //TODO: FIXME: remove transpose()
    work.temp_buffer.resize(work.temp_buffer.size(), 1);

    // Update Siblings
    work.temp_buffer *= kroneckerProduct(dad_v.matrix(), mom_v.matrix()).array();
    work.super[child] = work.temp_buffer / ((mat[child] *
                        work.lower[child].matrix()).array() + DNG_INDIVIDUAL_BUFFER_MIN);
    for(std::size_t i = 3; i < family.size(); ++i) {
        child = family[i];
        work.super[child] = work.temp_buffer / ((mat[child] *
                            work.lower[child].matrix()).array() + DNG_INDIVIDUAL_BUFFER_MIN);
        work.upper[child] = mat[child].transpose() * work.super[child].matrix();
    }
}
