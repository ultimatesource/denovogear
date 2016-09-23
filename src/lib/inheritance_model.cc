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

#include <dng/inheritance_model.h>
#include <iostream>
#include <stdexcept>

dng::InheritanceModel::~InheritanceModel() {
}

dng::InheritanceModel::InheritanceModel() {
}

void dng::InheritanceModel::parse_model(std::string &model_string) {

    char init = tolower(model_string.at(0));
    switch (init) {
        case 'a':
        case 'd':
            pattern = InheritancePattern::AUTOSOMAL;
            break;
        case 'm':
            pattern = InheritancePattern::MATERNAL;
            break;
        case 'p':
            pattern = InheritancePattern::PATERNAL;
            break;
        case 'x':
            pattern = InheritancePattern::X_LINKED;
            break;
        case 'y':
            pattern = InheritancePattern::Y_LINKED;
            break;
        case 'w':
            pattern = InheritancePattern::W_LINKED;
            break;
        case 'z':
            pattern = InheritancePattern::Z_LINKED;
            break;
        default:
            throw std::runtime_error(
                    "ERROR!! Inheritance model (" + model_string
                            + ") is not supported.\nSupported values are: "
                            + "[autosomal, default, xlinked, ylinked, wlinked, zlinked, maternal, paternal, mitochondria]");
            break;
    }
}

dng::InheritancePattern dng::InheritanceModel::GetInheritancePattern(){
    return pattern;
}
