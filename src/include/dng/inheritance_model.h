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
    AUTOSOMAL = 0,
    DEFAULT = 0,
    MATERNAL = 1,
    PATERNAL = 2,
    X_LINKED = 3,
    Y_LINKED = 4,
    W_LINKED = 5,
    Z_LINKED = 6,
    MITOCHONDRIA = 1
//        autosomal (the default)
//        xlinked (females have 2 copies, males have 1; males transmit to daughters, not to sons)
//        ylinked (males have 1 copy, only transmits it to sons)
//        wlinked (females have 1 copy, only transmited to daughters)
//        zlinked (males have 2 copies, females have 1; females transmit to sons, not to daughters)
//        maternal (transmitted by mother to child)
//        paternal (transmitter by father to child)
//        mitochondria  (transmitted by mother to child)
};

class InheritanceModel{
public:
    InheritanceModel();
    ~InheritanceModel();

    void parse_model(std::string &pattern_model);

    InheritancePattern GetInheritancePattern();

private:
    InheritancePattern pattern = InheritancePattern::AUTOSOMAL;
};



} // namespace dng


#endif /* INCLUDE_DNG_INHERITANCE_MODEL_H_ */
