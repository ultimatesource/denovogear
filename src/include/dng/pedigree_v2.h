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

#pragma once
#ifndef DNG_PEDIGREE_V2_H
#define DNG_PEDIGREE_V2_H

//#include <functional>
//#include <cmath>
//#include <array>

//#include <dng/matrix.h>
//#include <dng/io/ped.h>
//#include <dng/newick.h>
//#include <dng/read_group.h>
//#include <dng/peeling.h>

#include <dng/pedigree.h>
namespace dng {

class PedigreeV2 : public Pedigree {

    typedef boost::property_map<Graph, boost::edge_type_t>::type PropEdgeType;
    typedef boost::property_map<Graph, boost::edge_length_t>::type PropEdgeLength;
    typedef boost::property_map<Graph, boost::vertex_label_t>::type PropVertexLabel;
    typedef boost::property_map<Graph, boost::vertex_group_t>::type PropVertexGroup;
    typedef boost::property_map<Graph, boost::vertex_index_t>::type PropVertexIndex;

    typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;


    typedef std::vector<std::vector<boost::graph_traits<dng::Graph>::edge_descriptor>> family_labels_t;
//    auto edge_types = get(edge_type, pedigree_graph);
//    auto lengths = get(edge_length, pedigree_graph);
//    auto labels = get(vertex_label, pedigree_graph);
//    auto groups = get(vertex_group, pedigree_graph);
//    auto families = get(edge_family, pedigree_graph);
//    // Add the labels for the germline nodes
//


    enum class InheritancePattern : int {
        AUTOSOMAL = 0,
        DEFAULT = 0,
        MATERNAL = 1,
        PATERNAL = 2,
        X_LINKED = 3,
        Y_LINKED = 4,
        W_LINKED = 5,
        Z_LINKED = 6,

//        autosomal (the default)
//        xlinked (females have 2 copies, males have 1; males transmit to daughters, not to sons)
//        ylinked (males have 1 copy, only transmits it to sons)
//        wlinked (females have 1 copy, only transmited to daughters)
//        zlinked (males have 2 copies, females have 1; females transmit to sons, not to daughters)
//        maternal (transmitted by mother to child)
//        paternal (transmitter by father to child)
    };

    //TODO: struct FamilyInfo/Family structure.
    //Op1: A struct to record info in each family. family_t and ops
    //Op2: Another struct to group families together, include pivots and root?
    enum class FamilyType : int {
        PAIR = 2,
        TRIO = 3
    };
    struct FamilyInfo {
        FamilyType family_type;
        family_labels_t family_labels;//(num_families);
        std::vector<vertex_t> pivots;//(num_families, dummy_index);
    };

public:
    bool Construct(const io::Pedigree &pedigree, dng::ReadGroups &rgs,
                   double mu, double mu_somatic, double mu_library);

    //TODO: Eventually replace with this, or pass inheritance with a different method
    bool Construct(const io::Pedigree &pedigree, dng::ReadGroups &rgs, const InheritancePattern &pattern,
                   double mu, double mu_somatic, double mu_library);

    bool Equal(Pedigree &other_ped);


    const std::vector<peel::family_members_t>& inspect_family_members() const {
        return family_members_;
    }
    const std::vector<decltype(peel::op::NUM)>& inspect_peeling_ops() const {
        return peeling_ops_;
    };
    const std::vector<decltype(peel::op::NUM)>& inspect_peeling_functions_ops() const {
        return peeling_functions_ops_;
    };
//    // The original, simplified peeling operations
//    std::vector<decltype(peel::op::NUM)> peeling_ops_;
//    // The modified, "faster" operations
//    std::vector<decltype(peel::op::NUM)> peeling_functions_ops_;
//    // Array of functions that will be called to perform the peeling
//    std::vector<peel::function_t> peeling_functions_;
//    std::vector<peel::function_t> peeling_reverse_functions_;
//
//    // The arguments to a peeling operation
//    std::vector<peel::family_members_t> family_members_;

protected:


    void ParseIoPedigree(dng::Graph &pedigree_graph, const dng::io::Pedigree &pedigree);

    void AddLibrariesFromReadGroups(dng::Graph &pedigree_graph, const dng::ReadGroups &rgs,
                                                PropVertexLabel &labels);

    void ConnectSomaticToLibraries(dng::Graph &pedigree_graph, const ReadGroups &rgs,
                                   const PropVertexLabel &labels);

    void UpdateEdgeLengths(dng::Graph &pedigree_graph, double mu, double mu_somatic, double mu_library);

    void SimplifyPedigree(dng::Graph &pedigree_graph, const PropEdgeType &edge_types, const PropEdgeLength &lengths);

    void UpdateLabelsNodeIds(dng::Graph &pedigree_graph, dng::ReadGroups rgs, std::vector<size_t> &node_ids);

    void EraseRemovedLibraries(dng::ReadGroups &rgs, std::vector<size_t> &node_ids);

    void CreateFamiliesInfo(dng::Graph &pedigree_graph, family_labels_t &family_labels, std::vector<vertex_t> &pivots);

    void CreatePeelingOps(dng::Graph &pedigree_graph, const std::vector<size_t> &node_ids,
                          family_labels_t &family_labels, std::vector<vertex_t> &pivots);

private:
    void PrintDebugEdges(const std::string &prefix, const dng::Graph &pedigree_graph);


    const vertex_t DUMMY_INDEX = 0;
};

}// namespace dng


#endif // DNG_PEDIGREE_H
