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

#include <dng/pedigree.h>
#include <dng/graph.h>
#include <dng/mutation.h>

#include <algorithm>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/range/algorithm/replace.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>

#define DNG_GL_PREFIX "GL-"
#define DNG_SM_PREFIX "SM-" // define also in newick.cc
#define DNG_LB_PREFIX "LB-"

/*
RULES FOR LINKING READ GROUPS TO PEOPLE.

0) Read all read groups in bam files. Map RG to Library to Sample.

1) If tissue information is present, build a tissue tree connecting samples
   to zygotic genotypes.  Build mutation matrices for every unique branch
   length.

2) If no tissue information is present in pedigree, check to see if there is a
   sample with the same label as the individual.  If so, connect that sample to
   the individual with a somatic branch of length 1.

3) If a sample has multiple libraries, connect the libraries to sample with a
   library branch.

4) If a library has multiple read-groups, concat the read-groups.

*/

bool dng::Pedigree::Construct(const io::Pedigree &pedigree,
                              dng::ReadGroups &rgs,
                              double mu, double mu_somatic, double mu_library) {
    using namespace boost;
    using namespace std;

    // Setup counts
    std::size_t num_members = pedigree.member_count();
    // Count the founders
    first_founder_ = 0;
    first_somatic_ = num_members;
    for(first_nonfounder_ = first_founder_; first_nonfounder_ < num_members;
            ++first_nonfounder_) {
        if(pedigree.table()[first_nonfounder_].dad != 0 &&
                pedigree.table()[first_nonfounder_].mom != 0) {
            break;
        }
    }
    std::size_t num_libraries = rgs.libraries().size();

    // Construct a graph of the pedigree and somatic information
    Graph pedigree_graph(num_members);
    graph_traits<Graph>::edge_iterator ei, ei_end;
    graph_traits<Graph>::vertex_iterator vi, vi_end;

    auto edge_types = get(edge_type, pedigree_graph);
    auto lengths = get(edge_length, pedigree_graph);
    auto labels = get(vertex_label, pedigree_graph);
    auto groups = get(vertex_group, pedigree_graph);
    auto families = get(edge_family, pedigree_graph);

    // Add the labels for the germline nodes
    labels[0] = DNG_GL_PREFIX "unknown";
    for(size_t i = 1; i < num_members; ++i) {
        labels[i] = DNG_GL_PREFIX + pedigree.name(i);
    }

    // Go through rows and construct the pedigree part.
    vertex_t dummy_index = 0;
    for(auto &row : pedigree.table()) {
        vertex_t child = row.child;
        vertex_t dad = row.dad;
        vertex_t mom = row.mom;

        if(child == 0) {
            continue;
        }
        if(dad == mom && dad != 0) {
            // Selfing is not supported
            return false;
        }

        // check to see if mom and dad have been seen before
        auto id = edge(dad, mom, pedigree_graph);
        if(!id.second) {
            add_edge(dad, mom, EdgeType::Spousal, pedigree_graph);
        }
        // add the meiotic edges
        add_edge(mom, child, {EdgeType::Meiotic, 1.0f}, pedigree_graph);
        add_edge(dad, child, {EdgeType::Meiotic, 1.0f}, pedigree_graph);

        // Process newick file
        int res = newick::parse(row.sample_tree, child, pedigree_graph);
        if(res == 0) {
            // this line has a blank somatic line, so use the name from the pedigree
            vertex_t v = add_vertex(DNG_SM_PREFIX + pedigree.name(child), pedigree_graph);
            add_edge(child, v, {EdgeType::Mitotic, 1.0f}, pedigree_graph);
        } else if(res == -1) {
            throw std::runtime_error(
                "unable to parse somatic data for individual '" +
                pedigree.name(child) + "'.");
        }
    }
    // Remove the dummy individual from the graph
    clear_vertex(dummy_index, pedigree_graph);

    first_library_ = num_vertices(pedigree_graph);

    // Add library nodes to graph
    for(auto && a : rgs.libraries()) {
        vertex_t v = add_vertex(pedigree_graph);
        labels[v] = DNG_LB_PREFIX + a;
    }

    num_nodes_ = num_vertices(pedigree_graph);

    // Connect somatic samples to libraries
    for(vertex_t v = first_somatic_; v < first_library_; ++v) {
        if(labels[v].empty()) {
            continue;
        }

        auto r = rgs.data().get<rg::sm>().equal_range(labels[v].c_str() + strlen(
                     DNG_SM_PREFIX));
        for(; r.first != r.second; ++r.first) {
            vertex_t w = first_library_ + rg::index(rgs.libraries(),
                                                    r.first->library);
            if(!edge(v, w, pedigree_graph).second) {
                add_edge(v, w, {EdgeType::Library, 1.0f}, pedigree_graph);
            }
        }
    }

    // Update edge lengths
    for(tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
        if(edge_types[*ei] == EdgeType::Meiotic) {
            lengths[*ei] *= mu;
        } else if(edge_types[*ei] == EdgeType::Mitotic) {
            lengths[*ei] *= mu_somatic;
        } else if(edge_types[*ei] == EdgeType::Library) {
            lengths[*ei] *= mu_library;
        }
    }

    // Simplify Pedigree
    for(vertex_t w = first_library_; w > first_founder_; --w) {
        vertex_t v = w - 1;
        auto rng = out_edges(v, pedigree_graph);
        size_t children = 0, ancestors = 0, spouses = 0;
        auto it = rng.first;
        for(; it != rng.second; ++it) {
            edge_t e = *it;
            if(edge_types[e] == EdgeType::Spousal) {
                spouses += 1;
            } else if(target(e, pedigree_graph) > v) {
                children += 1;
            } else {
                ancestors += 1;
            }
        }
        if(children == 0) {
            // this node has no descendants
            clear_vertex(v, pedigree_graph);
        } else if(children >= 2 || spouses != 0) {
            /*noop*/;
        } else if(ancestors == 1) {
            // this node has 1 ancestor and 1 descendant
            edge_t e1 = *rng.first;
            edge_t e2 = *(rng.first + 1);
            vertex_t parent = target(e1, pedigree_graph);
            vertex_t child = target(e2, pedigree_graph);
            if(child < parent) {
                swap(child, parent);
                swap(e1, e2);
            }
            add_edge(parent, child, {edge_types[e1], lengths[e1] + lengths[e2]},
                     pedigree_graph);
            clear_vertex(v, pedigree_graph);
        } else if(ancestors == 2) {
            edge_t e1 = *rng.first;
            edge_t e2 = *(rng.first + 1);
            edge_t e3 = *(rng.first + 2);
            vertex_t parent1 = target(e1, pedigree_graph);
            vertex_t parent2 = target(e2, pedigree_graph);
            vertex_t child = target(e3, pedigree_graph);
            if(child < parent1) {
                swap(child, parent1);
                swap(e3, e1);
            }
            if(child < parent2) {
                swap(child, parent2);
                swap(e3, e2);
            }
            add_edge(parent1, child, {edge_types[e1], lengths[e1] + lengths[e3]},
                     pedigree_graph);
            add_edge(parent2, child, {edge_types[e2], lengths[e2] + lengths[e3]},
                     pedigree_graph);
            clear_vertex(v, pedigree_graph);
        }
    }

    vector<size_t> node_ids(num_nodes_, -1);
    labels_.clear();
    labels_.reserve(128);
    size_t vid = 0;
    for(std::size_t u = 0; u < num_nodes_ ; ++u) {
        if(out_degree(u, pedigree_graph) == 0) {
            continue;
        }
        if(!labels[u].empty()) {
            labels_.push_back(labels[u]);
        } else {
            labels_.push_back(DNG_SM_PREFIX "unnamed_node_" + util::to_pretty(vid));
        }
        node_ids[u] = vid++;
    }
    num_nodes_ = vid;

    // Erase Libraries that have been removed
    {
        vector<string> bad_libraries;
        bad_libraries.reserve(num_libraries);
        auto it = rgs.libraries().begin();
        for(size_t u = first_library_; u < node_ids.size(); ++u, ++it) {
            if(node_ids[u] != -1) {
                continue;
            }
            bad_libraries.push_back(*it);
        }
        rgs.EraseLibraries(bad_libraries);
    }

    // Update indexes
    auto update_position = [&node_ids, vid](size_t pos) -> size_t {
        for(; pos < node_ids.size() && node_ids[pos] == -1; ++pos)
            /*noop*/;
        return (pos < node_ids.size()) ? node_ids[pos] : vid;
    };

    first_founder_ = update_position(first_founder_);
    first_nonfounder_ = update_position(first_nonfounder_);
    first_somatic_ = update_position(first_somatic_);
    first_library_ = update_position(first_library_);

    // Reset Family Information
    roots_.clear();
    roots_.reserve(16);
    family_members_.clear();
    family_members_.reserve(128);
    peeling_ops_.clear();
    peeling_ops_.reserve(128);

    // Calculate the connected components.  This defines independent sections
    // of the graph.
    std::size_t num_groups = connected_components(pedigree_graph, groups);

    // Calculate the biconnected components and articulation points.
    // This defines "nuclear" families and pivot individuals.
    // Nodes which have no edges will not be part of any family.
    vector<vertex_t> articulation_vertices;
    std::size_t num_families = biconnected_components(pedigree_graph, families,
                               back_inserter(articulation_vertices)).first;

    // Determine which edges belong to which nuclear families.
    typedef vector<vector<graph_traits<Graph>::edge_descriptor>>
            family_labels_t;
    family_labels_t family_labels(num_families);
    for(tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
        family_labels[families[*ei]].push_back(*ei);
    }

    // Determine the last family in each group.  All singleton groups will have
    // a value of -1 since they have no family assignment.
    typedef deque<std::size_t> root_families_t;
    root_families_t root_families(num_groups, -1);
    for(std::size_t f = 0; f < family_labels.size(); ++f) {
        // last one wins
        auto first_edge = family_labels[f][0];
        auto src_vertex = source(first_edge, pedigree_graph);
        root_families[groups[src_vertex]] = f;
    }

    // Identify the pivot for each family.
    // The pivot will be the last art. point that has an edge in
    // the group.  The pivot of the last group doesn't matter.
    vector<vertex_t> pivots(num_families, dummy_index);
    for(auto a : articulation_vertices) {
        graph_traits<Graph>::out_edge_iterator ei, ei_end;
        for(tie(ei, ei_end) = out_edges(a, pedigree_graph); ei != ei_end; ++ei) {
            // Just overwrite existing value so that the last one wins.
            pivots[families[*ei]] = a;
        }
    }
    // Root Pivots are special
    for(auto f : root_families) {
        if(f == -1) { // skip all singleton groups
            continue;
        }
        pivots[f] = dummy_index;
    }

    // Resize the information in the pedigree
    transitions_.resize(num_nodes_);

    for(std::size_t i = first_founder_; i < first_nonfounder_; ++i) {
        transitions_[i] = {TransitionType::Founder, static_cast<size_t>(-1), static_cast<size_t>(-1), 0, 0};
    }

    // Detect Family Structure and pivot positions
    for(std::size_t k = 0; k < family_labels.size(); ++k) {
        auto &family_edges = family_labels[k];
        // Sort edges based on type and target
        boost::sort(family_edges, [&](edge_t x, edge_t y) -> bool { return
                        (edge_types(x) < edge_types(y)) &&
                        (target(x, pedigree_graph) < target(y, pedigree_graph)); });

        // Find the range of the parent types
        auto pos = boost::find_if(family_edges, [&](edge_t x) -> bool {
            return (edge_types(x) != EdgeType::Spousal); });
        size_t num_parent_edges = distance(family_edges.begin(), pos);

        // Check to see what type of graph we have
        if(num_parent_edges == 0) {
            // If we do not have a parent-child single branch,
            // we can't construct the pedigree.
            // TODO: Write Error message
            if(family_edges.size() != 1) {
                return false;
            }
            // Create a mitotic peeling operation.
            size_t parent = node_ids[source(*pos, pedigree_graph)];
            size_t child = node_ids[target(*pos, pedigree_graph)];

            TransitionType tt = (edge_types(*pos) == EdgeType::Library) ?
                                TransitionType::Library : TransitionType::Somatic;
            transitions_[child] = {tt, parent, static_cast<size_t>(-1), lengths[*pos], 0};

            family_members_.push_back({parent, child});
            if(node_ids[pivots[k]] == child) {
                peeling_ops_.push_back(peel::op::DOWN);
            } else {
                peeling_ops_.push_back(peel::op::UP);
                if(pivots[k] == dummy_index) {
                    roots_.push_back(parent);
                }
            }
        } else if(num_parent_edges == 1) {
            // If this family contains no children, skip it
            if(pos == family_edges.end()) {
                continue;
            }
            // We have a nuclear family with 1 or more children
            size_t dad = node_ids[source(family_edges.front(), pedigree_graph)];
            size_t mom = node_ids[target(family_edges.front(), pedigree_graph)];
            family_members_.push_back({dad, mom});
            auto &family_members = family_members_.back();
            while(pos != family_edges.end()) {
                vertex_t child = node_ids[target(*pos, pedigree_graph)];
                transitions_[child] = { TransitionType::Germline, dad, mom,
                                        lengths[*pos], lengths[*(pos + 1)]
                                      };
                family_members.push_back(child); // Child
                // child edges come in pairs
                ++pos;
                assert(node_ids[target(*pos, pedigree_graph)] == child);
                ++pos;
            }
            if(node_ids[pivots[k]] == node_ids[dummy_index]) {
                // A family without a pivot is a root family
                peeling_ops_.push_back(peel::op::TOFATHER);
                roots_.push_back(family_members[0]);
            } else {
                auto pivot_pos = boost::range::find(family_members, node_ids[pivots[k]]);
                size_t p = distance(family_members.begin(), pivot_pos);
                if(p == 0) {
                    peeling_ops_.push_back(peel::op::TOFATHER);
                } else if(p == 1) {
                    peeling_ops_.push_back(peel::op::TOMOTHER);
                } else if(p == 2) {
                    peeling_ops_.push_back(peel::op::TOCHILD);
                } else {
                    peeling_ops_.push_back(peel::op::TOCHILD);
                    swap(family_members[p], family_members[2]);
                }
            }
        } else {
            // Not a zero-loop pedigree
            // TODO: write error message
            return false;
        }
    }

    ConstructPeelingMachine();

#ifdef DNG_DEVEL
    PrintMachine(cerr);
#else
    //PrintMachine(cerr);
#endif

    return true;
}

void dng::Pedigree::ConstructPeelingMachine() {
    using namespace dng::peel;
    peeling_functions_.clear();
    peeling_functions_.reserve(peeling_ops_.size());
    peeling_functions_ops_.reserve(peeling_ops_.size());
    peeling_reverse_functions_.reserve(peeling_ops_.size());
    std::vector<std::size_t> lower_written(num_nodes_, -1);
    for(std::size_t i = 0 ; i < peeling_ops_.size(); ++i) {
        auto a = peeling_ops_[i];
        const auto &fam = family_members_[i];
        auto w = fam[info[a].writes_to];
        int b = a;
        switch(a) {
        case op::DOWN:
            // If the lower of the parent has never been written to, we can use the fast version
            b = (lower_written[fam[0]] == -1) ? op::UPFAST + b : b;
            break;
        case op::TOCHILD:
            // If the we only have one child, we can use the fast version
            b = (fam.size() == 3) ? op::UPFAST + b : b;
            break;
        case op::TOMOTHER:
        case op::TOFATHER:
        case op::UP:
            // If the lower of the destination has never been written to, we can use the fast version
            b = (lower_written[w] == -1) ? op::UPFAST + b : b;
            break;
        default:
            assert(false); // should never get here
            break;
        }
        peeling_functions_ops_.push_back(static_cast<decltype(op::NUM)>(b));
        peeling_functions_.push_back(functions[b]);
        peeling_reverse_functions_.push_back(reverse_functions[b]);

        // If the operation writes to a lower value, make note of it
        if(info[a].writes_lower) {
            lower_written[w] = i;
        }
    }
}


void dng::Pedigree::PrintMachine(std::ostream &os) {
    os << "Init Op\n";
    for(int i = first_founder_; i < first_nonfounder_; ++i) {
        os << "\tw\tupper[" << i << "] // " << labels_[i] << "\n";
    }
    for(int i = first_library_; i < num_nodes_; ++i) {
        os << "\tw\tlower[" << i << "] // "
           << labels_[i]
           << "\n";
    }
    for(int i = 0; i < peeling_functions_ops_.size(); ++i) {
        const auto &fam = family_members_[i];
        os << "Peeling Op " << i + 1;
        switch(peeling_functions_ops_[i]) {
        case peel::op::DOWNFAST:
            os << " (DownFast)\n";
            os << "\tw\tupper[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            os << "\tr\tupper[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            break;
        case peel::op::DOWN:
            os << " (Down)\n";
            os << "\tw\tupper[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            os << "\tr\tupper[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            os << "\tr\tlower[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            break;
        case peel::op::UPFAST:
            os << " (UpFast)\n";
            os << "\tw\tlower[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            os << "\tr\tlower[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            break;
        case peel::op::UP:
            os << " (Up)\n";
            os << "\trw\tlower[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            os << "\tr\tlower[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            break;
        case peel::op::TOFATHERFAST:
            os << " (ToFatherFast)\n";
            os << "\tw\tlower[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            os << "\tr\tupper[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            os << "\tr\tlower[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            for(int j = 2; j < fam.size(); ++j) {
                os << "\tr\tlower[" << fam[j] << "] // " << labels_[fam[j]] << "\n";
            }
            break;
        case peel::op::TOFATHER:
            os << " (ToFather)\n";
            os << "\trw\tlower[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            os << "\tr\tupper[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            os << "\tr\tlower[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            for(int j = 2; j < fam.size(); ++j) {
                os << "\tr\tlower[" << fam[j] << "] // " << labels_[fam[j]] << "\n";
            }
            break;
        case peel::op::TOMOTHERFAST:
            os << " (ToMotherFast)\n";
            os << "\tw\tlower[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            os << "\tr\tupper[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            os << "\tr\tlower[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            for(int j = 2; j < fam.size(); ++j) {
                os << "\tr\tlower[" << fam[j] << "] // " << labels_[fam[j]] << "\n";
            }
            break;
        case peel::op::TOMOTHER:
            os << " (ToMother)\n";
            os << "\trw\tlower[" << fam[0] << "] // " << labels_[fam[1]] << "\n";
            os << "\tr\tupper[" << fam[1] << "] // " << labels_[fam[0]] << "\n";
            os << "\tr\tlower[" << fam[1] << "] // " << labels_[fam[0]] << "\n";
            for(int j = 2; j < fam.size(); ++j) {
                os << "\tr\tlower[" << fam[j] << "] // " << labels_[fam[j]] << "\n";
            }
            break;
        case peel::op::TOCHILDFAST:
            os << " (ToChildFast)\n";
            os << "\tw\tupper[" << fam[2] << "] // " << labels_[fam[2]] << "\n";
            os << "\tr\tupper[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            os << "\tr\tlower[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            os << "\tr\tupper[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            os << "\tr\tlower[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            break;
        case peel::op::TOCHILD:
            os << " (ToChild)\n";
            os << "\tw\tupper[" << fam[2] << "] // " << labels_[fam[2]] << "\n";
            os << "\tr\tupper[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            os << "\tr\tlower[" << fam[0] << "] // " << labels_[fam[0]] << "\n";
            os << "\tr\tupper[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            os << "\tr\tlower[" << fam[1] << "] // " << labels_[fam[1]] << "\n";
            for(int j = 3; j < fam.size(); ++j) {
                os << "\tr\tlower[" << fam[j] << "] // " << labels_[fam[j]] << "\n";
            }
            break;
        default:
            os << " (Unknown " << peeling_functions_ops_[i] << ")\n";
        }
    }

    os << "Root Op\n";
    for(int i = 0; i < roots_.size(); ++i) {
        os << "\tr\tupper[" << roots_[i] << "] // " << labels_[roots_[i]] << "\n";
        os << "\tr\tlower[" << roots_[i] << "] // " << labels_[roots_[i]] << "\n";
    }
}

void dng::Pedigree::PrintTable(std::ostream &os) {
    std::vector<int> write_low(num_nodes_, -1);
    std::vector<int> write_up(num_nodes_, -1);

    for(int i = first_founder_; i < first_nonfounder_; ++i) {
        write_up[i] = 0;
    }
    for(int i = first_library_; i < num_nodes_; ++i) {
        write_low[i] = 0;
    }
    for(int i = 0; i < peeling_ops_.size(); ++i) {
        const auto &fam = family_members_[i];
        switch(peeling_ops_[i]) {
        case peel::op::DOWN:
        case peel::op::DOWNFAST:
            write_up[fam[0]] = i + 1;
            break;
        case peel::op::TOCHILD:
        case peel::op::TOCHILDFAST:
            write_up[fam[2]] = i + 1;
            break;
        case peel::op::UP:
        case peel::op::UPFAST:
        case peel::op::TOFATHER:
        case peel::op::TOFATHERFAST:
            write_low[fam[0]] = i + 1;
            break;
        case peel::op::TOMOTHER:
        case peel::op::TOMOTHERFAST:
            write_low[fam[1]] = i + 1;
        default:
            break;
        }
    }
    os << "Node\tLower\tUpper\n";
    for(int i = 0; i < num_nodes_; ++i) {
        os << i << "\t" << write_low[i] << "\t" << write_up[i] << "\n";
    }
}

// void dng::Pedigree::PrintStates(std::ostream &os, double scale) {
//     for(int i = 0; i < lower_.size(); ++i) {
//         os << "Node " << i << " // " << labels_[i] << "\n";
//         os << "  Upper\t" << upper_[i].transpose() << "\n";
//         os << "  Lower\t" << lower_[i].transpose() << "\n";
//         auto p = upper_[i] * lower_[i];
//         auto s = p.sum();
//         os << "  Prod\t" << p.transpose() << "\n";
//         os << "  Cond\t" << p.transpose() / s << "\n";
//         os << "  LogSum\t" << log(s) + scale << "\n";
//         os << "\n";
//     }
// }

std::vector<std::string> dng::Pedigree::BCFHeaderLines() const {
    using namespace std;
    vector<string> ret = { "##META=<ID=FatherMR,Type=Float,Number=1,Description=\"Paternal mutation rate parameter\">",
                           "##META=<ID=MotherMR,Type=Float,Number=1,Description=\"Maternal mutation rate parameter\">",
                           "##META=<ID=OriginalMR,Type=Float,Number=1,Description=\"Somatic or library mutation rate parameter\">"
                         };

    string head{"##PEDIGREE=<ID="};
    for(size_t child = first_nonfounder_; child < num_nodes_; ++child) {
        auto parents = transitions_[child];
        if(parents.parent2 != -1) {
            ret.push_back(head + labels_[child]
                          + ",Father=" + labels_[parents.parent1]
                          + ",Mother=" + labels_[parents.parent2]
                          + ",FatherMR=" + util::to_pretty(parents.length1)
                          + ",MotherMR=" + util::to_pretty(parents.length2)
                          + ">"
                         );
        } else {
            ret.push_back(head + labels_[child]
                          + ",Original=" + labels_[parents.parent1]
                          + ",OriginalMR=" + util::to_pretty(parents.length1)
                          + ">"
                         );
        }
    }
    return ret;
}
