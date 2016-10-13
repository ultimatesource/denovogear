/*1
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
#ifndef INCLUDE_DNG_INHERITANCE_MODEL_H_
#define INCLUDE_DNG_INHERITANCE_MODEL_H_

#include <string>

namespace dng{

enum class InheritancePattern : int {
    DEFAULT = 0,
    AUTOSOMAL = 0, // default
    MITOCHONDRIA = 1, // transmitted by mother to child
    MATERNAL = 1, // transmitted by mother to child
    PATERNAL = 2, // transmitter by father to child
    X_LINKED = 3, // females have 2 copies, males have 1; males transmit to daughters, not to sons
    Y_LINKED = 4, // males have 1 copy, only transmits it to sons
    W_LINKED = 5, // females have 1 copy, only transmited to daughters
    Z_LINKED = 6  // males have 2 copies, females have 1; females transmit to sons, not to daughters

};



class InheritanceModel{
public:
    InheritanceModel();
    ~InheritanceModel();

    void parse_model(std::string &pattern_model);

    InheritancePattern GetInheritancePattern();


    const std::pair<std::string, InheritancePattern> INHERITANCE_KEYS[13] = {
        {"DEFAULT", InheritancePattern::AUTOSOMAL},
        {"AUTOSOMAL", InheritancePattern::AUTOSOMAL},
        {"MITOCHONDRIA", InheritancePattern::MATERNAL},
        {"MATERNAL", InheritancePattern::MATERNAL},
        {"PATERNAL", InheritancePattern::PATERNAL},
        {"X_LINKED", InheritancePattern::X_LINKED},
        {"Y_LINKED", InheritancePattern::Y_LINKED},
        {"W_LINKED", InheritancePattern::W_LINKED},
        {"Z_LINKED", InheritancePattern::Z_LINKED},
        {"XLINKED", InheritancePattern::X_LINKED},
        {"YLINKED", InheritancePattern::Y_LINKED},
        {"WLINKED", InheritancePattern::W_LINKED},
        {"ZLINKED", InheritancePattern::Z_LINKED}
    };



private:
    InheritancePattern pattern = InheritancePattern::AUTOSOMAL;
};



} // namespace dng


#endif /* INCLUDE_DNG_INHERITANCE_MODEL_H_ */
