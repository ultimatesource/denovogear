/*
 * Copyright (c) 2014-2015 Reed A. Cartwright
 * Copyright (c) 2016 Steven H. Wu
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *           Steven H. Wu <stevenwu@asu.edu>
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
#ifndef DNG_RELATIONSHIP_GRAPH_H
#define DNG_RELATIONSHIP_GRAPH_H

#include <cmath>
#include <string>
#include <vector>

#include <dng/matrix.h>
#include <dng/library.h>
#include <dng/peeling.h>
#include <dng/pedigree.h>
#include <dng/detail/graph.h>
#include <dng/detail/unit_test.h>

#include <dng/io/ped.h>

namespace dng {

enum struct InheritanceModel {
    Autosomal,    // default option
    Maternal,     // transmitted by mother to child
    Paternal,     // transmitter by father to child
    XLinked,      // females have 2 copies, males have 1; males transmit to daughters, not to sons
    YLinked,      // males have 1 copy, only transmits it to sons
    WLinked,      // females have 1 copy, only transmitted to daughters
    ZLinked,      // males have 2 copies, females have 1; females transmit to sons, not to daughters
    Unknown
};

InheritanceModel inheritance_model(const std::string &pattern);

std::string to_string(InheritanceModel model);

class RelationshipGraph {
public:
    template<typename T>
    using property_t = typename boost::property_map<detail::graph::Graph, T>::type;
    using PropEdgeType = property_t<boost::edge_type_t>; 
    using PropEdgeLength = property_t<boost::edge_length_t>;
    using PropVertexLabel = property_t<boost::vertex_label_t>;
    using PropVertexGroup = property_t<boost::vertex_group_t>;
    using PropVertexIndex = property_t<boost::vertex_index_t>;
    using PropVertexSex = property_t<boost::vertex_sex_t>;
    using IndexMap = property_t<boost::vertex_index_t>;

    using family_labels_t = std::vector<std::vector<boost::graph_traits<detail::graph::Graph>::edge_descriptor>>;

    enum struct TransitionType {
        Founder, Pair, Trio
    };

    struct transition_t {
        TransitionType type;
        std::size_t parent1;
        std::size_t parent2;
        double length1;
        double length2;

        bool is_germline;
        bool is_somatic;
        bool is_library;
    };

    bool Construct(const Pedigree& pedigree, const libraries_t& libs,
            InheritanceModel inheritance_model,
            double mu, double mu_somatic, double mu_library,
            bool normalize_somatic_trees);

    bool Construct(const Pedigree& pedigree, const libraries_t& libs,
            double mu, double mu_somatic, double mu_library,
            bool normalize_somatic_trees);

    double PeelForwards(peel::workspace_t &work,
                        const TransitionMatrixVector &mat) const {
        if(work.dirty_lower) {
            work.CleanupFast();
        }

        // Peel pedigree one family at a time
        for(std::size_t i = 0; i < peeling_functions_.size(); ++i) {
            (*peeling_functions_[i])(work, family_members_[i], mat);
        }

        // Sum over roots
        double ret = 0.0;
        for(auto r : roots_) {
            ret += log((work.lower[r] * work.upper[r]).sum());
        }
        
        return ret;
    }

    double PeelBackwards(peel::workspace_t &work,
                         const TransitionMatrixVector &mat) const {
        double ret = 0.0;
        // Divide by the likelihood
        for(auto r : roots_) {
            double sum = (work.lower[r] * work.upper[r]).sum();
            ret += log(sum);
            work.lower[r] /= sum;
        }

        for(std::size_t i = peeling_reverse_functions_.size(); i > 0; --i) {
            (*peeling_reverse_functions_[i - 1])(work, family_members_[i - 1], mat);
        }
        work.dirty_lower = true;
        return ret;
    }

    peel::workspace_t CreateWorkspace() const {
        peel::workspace_t work;
        work.Resize(num_nodes_);
        work.founder_nodes = std::make_pair(first_founder_, first_nonfounder_);
        work.germline_nodes = std::make_pair(first_founder_, first_somatic_);
        work.somatic_nodes = std::make_pair(first_somatic_, first_library_);
        work.library_nodes = std::make_pair(first_library_, num_nodes_);

        work.ploidies = ploidies_;

        return work;
    }

    std::vector<std::string> BCFHeaderLines() const;

    const std::vector<transition_t> &transitions() const { return transitions_; }
    const std::vector<std::string> &labels() const { return labels_; }
    const std::vector<int> &ploidies() const { return ploidies_; }

    const transition_t & transition(size_t pos) const { return transitions_[pos]; }
    const std::string & label(size_t pos) const { return labels_[pos]; }
    int ploidy(size_t pos) const { return ploidies_[pos]; }

    size_t num_nodes() const { return num_nodes_; }
    std::pair<size_t, size_t> library_nodes() const { return {first_library_, num_nodes_}; }

    const std::vector<std::string> &library_names() const { return library_names_; }

protected:
    using Graph = dng::detail::graph::Graph;
    using vertex_t = dng::detail::graph::vertex_t;

    InheritanceModel inheritance_model_{InheritanceModel::Autosomal};

    // node structure:
    // founder germline, non-founder germline, somatic, library
    std::size_t num_nodes_{0};        // total number of nodes
    std::size_t first_founder_{0};    // start of founder germline
    std::size_t first_nonfounder_{0}; // start of non-founder germline
    std::size_t first_somatic_{0};    // start of somatic nodes
    std::size_t first_library_{0};    // start of libraries

    std::vector<std::size_t> roots_;

    // Pedigree Structure
    std::vector<std::string> labels_;
    std::vector<int> ploidies_;
    std::vector<transition_t> transitions_;

    // The original, simplified peeling operations
    std::vector<peel::Op> peeling_ops_;
    // The modified, "faster" operations
    std::vector<peel::Op> peeling_functions_ops_;
    // Array of functions that will be called to perform the peeling
    std::vector<peel::function_t> peeling_functions_;
    std::vector<peel::function_t> peeling_reverse_functions_;

    // The arguments to a peeling operation
    std::vector<peel::family_members_t> family_members_;

    std::vector<std::string> library_names_;

    void ConstructPeelingMachine();

    std::vector<size_t> ConstructNodes(const Graph &pedigree_graph);

    using pivots_t = std::vector<boost::optional<vertex_t>>;

    void CreateFamiliesInfo(Graph &pedigree_graph,
            family_labels_t *family_labels, pivots_t *pivots);

    void CreatePeelingOps(const Graph &pedigree_graph,
            const std::vector<size_t> &node_ids,
            family_labels_t &family_labels,
            const pivots_t &pivots);

private:
    void ClearFamilyInfo();

    DNG_UNIT_TEST_CLASS(unittest_dng_relationship_graph);
};

template<typename A, typename M>
inline
RelationshipGraph create_relationship_graph(const A& arg, M* mpileup) {
    assert(mpileup != nullptr);
    // Parse Pedigree from File
    Pedigree ped = io::parse_ped(arg.ped);
    // Construct peeling algorithm from parameters and pedigree information
    RelationshipGraph relationship_graph;
    if (!relationship_graph.Construct(ped, mpileup->libraries(), inheritance_model(arg.model),
                                      arg.mu, arg.mu_somatic, arg.mu_library,
                                      arg.normalize_somatic_trees)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                 "possible non-zero-loop relationship_graph.");
    }
    // Select libraries in the input that are used in the pedigree
    mpileup->SelectLibraries(relationship_graph.library_names());

    return relationship_graph;
}

}; // namespace dng

#endif // DNG_RELATIONSHIP_GRAPH_H
