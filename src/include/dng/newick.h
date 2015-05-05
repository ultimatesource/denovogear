/*
 * Copyright (c) 2014 Reed A. Cartwright
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
#ifndef DNG_NEWICK_H
#define DNG_NEWICK_H

#include <dng/graph.h>

namespace dng {
namespace newick {

typedef dng::Graph tree_t;
typedef boost::graph_traits<tree_t>::vertex_descriptor node_t;
typedef boost::graph_traits<tree_t>::edge_descriptor branch_t;

int parse(const std::string &text, node_t root, tree_t &tree);

}
};

#endif // DNG_NEWICK_H

