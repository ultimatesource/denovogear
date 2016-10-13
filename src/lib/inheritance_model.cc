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

#include <dng/utility.h>

dng::InheritanceModel::~InheritanceModel() {
}

dng::InheritanceModel::InheritanceModel() {
}

void dng::InheritanceModel::parse_model(std::string &model_string) {

    auto match = dng::utility::key_switch_tuple(model_string, INHERITANCE_KEYS,
            std::pair<std::string, InheritancePattern>{"UNMATCHED", InheritancePattern::DEFAULT});
    if ("UNMATCHED" == match.first ){
        std::cerr << "Warning!! Inheritance model (" + model_string
                + ") is not supported.\nDefault inheritance model is used.\nSupported values are: "
                + "[default, autosomal, mitochondria, maternal, paternal, x_linked, y_linked, w_linked, z_linked]"
                << std::endl;
    }
    pattern = match.second;
}

dng::InheritancePattern dng::InheritanceModel::GetInheritancePattern(){
    return pattern;
}
