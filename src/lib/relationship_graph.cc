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
#include <dng/relationship_graph.h>
#include <dng/utility.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/range/algorithm/find.hpp>

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


bool dng::RelationshipGraph::Construct(const io::Pedigree& pedigree,
        dng::ReadGroups& rgs, double mu, double mu_somatic, double mu_library) {
    return Construct(pedigree, rgs, InheritancePattern::AUTOSOMAL, mu,
                     mu_somatic, mu_library);
}

bool dng::RelationshipGraph::Construct(const io::Pedigree& pedigree,
        dng::ReadGroups& rgs, InheritancePattern inheritance_pattern,
        double mu, double mu_somatic, double mu_library) {

    using namespace std;


    //PR_NOTE(SW): overkill, but it forced to be Autosomal for now, other model are not fully implemented
    inheritance_pattern = InheritancePattern::AUTOSOMAL;//

    SetupFirstNodeIndex(pedigree);

    // Construct a graph of the pedigree and somatic information
    Graph pedigree_graph(first_somatic_);
    PrintDebugEdges("========== VERISON 2 =========\ninit pedigree", pedigree_graph);

    auto edge_types = get(boost::edge_type, pedigree_graph);
    auto lengths = get(boost::edge_length, pedigree_graph);
    auto labels = get(boost::vertex_label, pedigree_graph);
    auto sex  = get(boost::vertex_sex, pedigree_graph);

    ParseIoPedigree(pedigree_graph, pedigree);
    AddLibrariesFromReadGroups(pedigree_graph, rgs);

    UpdateEdgeLengths(pedigree_graph, mu, mu_somatic, mu_library);
    SimplifyPedigree(pedigree_graph);

    std::vector<size_t> node_ids(num_nodes_, -1);//num_nodes used twice for different contex
    UpdateLabelsNodeIds(pedigree_graph, rgs, node_ids);

    family_labels_t family_labels;//(num_families);
    std::vector<vertex_t> pivots;//(num_families, dummy_index);
    CreateFamiliesInfo(pedigree_graph, family_labels, pivots);

    if (inheritance_pattern == InheritancePattern::AUTOSOMAL) {
        CreatePeelingOps(pedigree_graph, node_ids, family_labels, pivots);
        ConstructPeelingMachine();
    }
    else {
        std::cerr << "You should never get here, other inheritance pattern not implemented yet!! EXTI!" << std::endl;
        std::exit(9);
    }


#ifdef DEBUG_RGRAPH
    PrintMachine(cout);
    PrintTable(cout);
#endif


#ifdef DNG_DEVEL
    PrintMachine(cerr);
#else
    //PrintMachine(cerr);
#endif

    return true;
}


void dng::RelationshipGraph::ConstructPeelingMachine() {
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


void dng::RelationshipGraph::PrintMachine(std::ostream &os) {
    os << "Init Op\n";
    for(int i = first_founder_; i < first_nonfounder_; ++i) {
        os << "\tw\tupper[" << i << "] // " << labels_[i] << "\n";
    }
    for(int i = first_library_; i < num_nodes_; ++i) {
        os << "\tw\tlower[" << i << "] // "
           << labels_[i] << "\n";
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

void dng::RelationshipGraph::PrintTable(std::ostream &os) {
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

std::vector<std::string> dng::RelationshipGraph::BCFHeaderLines() const {
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
                          + ",FatherMR=" + utility::to_pretty(parents.length1)
                          + ",MotherMR=" + utility::to_pretty(parents.length2)
                          + ">"
                         );
        } else {
            ret.push_back(head + labels_[child]
                          + ",Original=" + labels_[parents.parent1]
                          + ",OriginalMR=" + utility::to_pretty(parents.length1)
                          + ">"
                         );
        }
    }
    return ret;
}



void dng::RelationshipGraph::SetupFirstNodeIndex(const io::Pedigree &pedigree){
    first_founder_ = 0; //TODO(SW): Can this be something other than 0??
    first_somatic_ = pedigree.member_count();

    for (first_nonfounder_ = first_founder_; first_nonfounder_ < first_somatic_;
            ++first_nonfounder_) {
        if (pedigree.table()[first_nonfounder_].dad != 0
                && pedigree.table()[first_nonfounder_].mom != 0) {
            break;
        }
    }

}
void dng::RelationshipGraph::ParseIoPedigree(dng::Graph &pedigree_graph,
        const dng::io::Pedigree &pedigree) {

    auto labels = get(boost::vertex_label, pedigree_graph);
    auto sex  = get(boost::vertex_sex, pedigree_graph);
    labels[0] = DNG_GL_PREFIX "unknown";
    for (size_t i = 1; i < first_somatic_; ++i) {
        labels[i] = DNG_GL_PREFIX + pedigree.name(i);
    }

    PrintDebugEdges("Start", pedigree_graph);
    for (auto &row : pedigree.table()) {
        vertex_t child = row.child;
        vertex_t dad = row.dad;
        vertex_t mom = row.mom;
        sex[child] = row.sex;
#ifdef DEBUG_RGRAPH
        std::cout << "\n==START\n===Child_Dad_MoM_Sex: " << child << "\t"
            << dad << "\t" << mom << "\t"
            << (int) row.sex << std::endl;
#endif

        if (child == 0) {
            continue;
        }
        if (dad == mom && dad != 0) {
            // Selfing is not supported
            throw std::runtime_error("Unable to construct peeler for pedigree; selfing is not supported");
//                                             "possible non-zero-loop pedigree.");
        }

        if ((dad != 0 && sex[dad] != dng::io::Pedigree::Sex::Male)
                || (mom != 0 && sex[mom] != dng::io::Pedigree::Sex::Female)) {
            throw std::runtime_error(
                    "Error: Sex error, dad!=Male || mom!=Female");
        }

        // check to see if mom and dad have been seen before
        auto id = edge(dad, mom, pedigree_graph);
        if (!id.second) { //Connect dad-mom to make a dependend trio
            add_edge(dad, mom, EdgeType::Spousal, pedigree_graph);
        }

#ifdef DEBUG_RGRAPH
    std::cout << "===after spoesal: V E:" << num_vertices(pedigree_graph) << "\t" << num_edges(pedigree_graph) << std::endl;
    int sex2 = static_cast<int>( sex[child] );
        std::cout << "\n===Child_Dad_MoM_sex: " << child << "\t" << dad
                << "\t" << mom << "\tsex: " << (int) sex[child] << "\t" << sex2
                << std::endl;
#endif
        // add the meiotic edges
        add_edge(mom, child, {EdgeType::Meiotic, 1.0f}, pedigree_graph);
        add_edge(dad, child, {EdgeType::Meiotic, 1.0f}, pedigree_graph);

#ifdef DEBUG_RGRAPH
    std::cout << "===after add_edge: V E:" << num_vertices(pedigree_graph) << "\t" << num_edges(pedigree_graph) << std::endl;
    std::cout << "===Child_Dad_MoM: " << child << "\t" << dad << "\t" << mom << std::endl;
#endif

        // Process newick file
        // TODO(SW): What should newick::parse do? parse "tree" only? Or add vertex/edge/sex as well?
        // TODO(SW): Refactor newick::parse and the if_else(res==?) condition. The same operation performed at two different places
        //HACK(SW): add sex without touch newick::parse at this stage
        int current_index = num_vertices(pedigree_graph);
        int res = newick::parse(row.sample_tree, child, pedigree_graph);

#ifdef DEBUG_RGRAPH
    std::cout << "===after parse: V E:" << num_vertices(pedigree_graph) << "\t" << num_edges(pedigree_graph) <<
    std::endl;
#endif

        for (int i = current_index; i < num_vertices(pedigree_graph); ++i) {
//            std::cout << i << "\t" << num_vertices(pedigree_graph) << std::endl;
            sex[i] = sex[child];
        }

        if (res == 0) {
            // this line has a blank somatic line, so use the name from the pedigree
            vertex_t v = add_vertex(DNG_SM_PREFIX + pedigree.name(child), pedigree_graph);
            add_edge(child, v, {EdgeType::Mitotic, 1.0f}, pedigree_graph);
            sex[v] = sex[child];

#ifdef DEBUG_RGRAPH
    std::cout << "===add res==0: " << child << "\t" << v << "\t" << pedigree.name(child) << std::endl;
#endif

        } else if (res == -1) {
            throw std::runtime_error("unable to parse somatic data for individual '" +
                                     pedigree.name(child) + "'.");
        }

#ifdef DEBUG_RGRAPH
    std::cout << "Loop END: V E:" << num_vertices(pedigree_graph) << "\t" << num_edges(pedigree_graph) << std::endl;
#endif

    }
    PrintDebugEdges("Parsed io::pedigree, before cleas_vertex dummy", pedigree_graph);
    clear_vertex(DUMMY_INDEX, pedigree_graph);
    PrintDebugEdges("Parsed io::pedigree, after cleas_vertex dummy", pedigree_graph);
}

void dng::RelationshipGraph::AddLibrariesFromReadGroups(
        Graph &pedigree_graph, const dng::ReadGroups &rgs) {

    auto labels = get(boost::vertex_label, pedigree_graph);
    first_library_ = num_vertices(pedigree_graph);

    // Add library nodes to graph
    for (auto &&a : rgs.libraries()) {
        vertex_t v = add_vertex(pedigree_graph);
        labels[v] = DNG_LB_PREFIX + a;

#ifdef DEBUG_RGRAPH
        std::cout << "Add lib: " << labels[v] << std::endl;
#endif
    }
    PrintDebugEdges("After add libs", pedigree_graph);

    ConnectSomaticToLibraries(pedigree_graph, rgs, labels);
    num_nodes_ = num_vertices(pedigree_graph);

    PrintDebugEdges("After connect somatic", pedigree_graph);
}

void dng::RelationshipGraph::ConnectSomaticToLibraries(
        dng::Graph &pedigree_graph, const ReadGroups &rgs,
        const PropVertexLabel &labels) {

    auto sex  = get(boost::vertex_sex, pedigree_graph);
    const size_t STRLEN_DNG_SM_PREFIX = strlen(DNG_SM_PREFIX);

    for (vertex_t v = (vertex_t) first_somatic_; v < first_library_; ++v) {
        if (labels[v].empty()) {
            continue;
        }

        //using orders from vcf files rgs
        auto r = rgs.data().get<rg::sm>().equal_range(
                labels[v].c_str() + STRLEN_DNG_SM_PREFIX);
        for (; r.first != r.second; ++r.first) {
            vertex_t w = first_library_ + rg::index(rgs.libraries(), r.first->library);
            sex[w] = sex[v];
            if (!edge(v, w, pedigree_graph).second) {
                add_edge(v, w, {EdgeType::Library, 1.0f}, pedigree_graph);
            }
        }
    }
}



void dng::RelationshipGraph::UpdateEdgeLengths(dng::Graph &pedigree_graph,
        double mu_meiotic, double mu_somatic, double mu_library) {
    boost::graph_traits<dng::Graph>::edge_iterator ei, ei_end;
    auto edge_types = get(boost::edge_type, pedigree_graph);
    auto lengths = get(boost::edge_length, pedigree_graph);

    for (tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
        if (edge_types[*ei] == dng::graph::EdgeType::Meiotic) {
            lengths[*ei] *= mu_meiotic;
        } else if (edge_types[*ei] == dng::graph::EdgeType::Mitotic) {
            lengths[*ei] *= mu_somatic;
        } else if (edge_types[*ei] == dng::graph::EdgeType::Library) {
            lengths[*ei] *= mu_library;
        }
    }

#ifdef DEBUG_RGRAPH
    PrintDebugEdges("After UpdateEdgeLengths", pedigree_graph);
#endif
}

void dng::RelationshipGraph::SimplifyPedigree(dng::Graph &pedigree_graph) {

    auto edge_types = get(boost::edge_type, pedigree_graph);
    auto lengths = get(boost::edge_length, pedigree_graph);

    for (auto w = first_library_; w > first_founder_; --w) {
        vertex_t v = (vertex_t) (w - 1);
        size_t children = 0, ancestors = 0, spouses = 0;

        auto rng = out_edges(v, pedigree_graph);
        for (auto it = rng.first; it != rng.second; ++it) {
            edge_t e = *it;
            if (edge_types[e] == EdgeType::Spousal) {
                spouses += 1;
            } else if (target(e, pedigree_graph) > v) {
                children += 1;
            } else {
                ancestors += 1;
            }
        }
        if (children == 0) {
            // this node has no descendants
            clear_vertex(v, pedigree_graph);
#ifdef DEBUG_RGRAPH
    std::cout << "Remove child==0:" << v << std::endl;
#endif
        } else if (children >= 2 || spouses != 0) {
            /*noop*/;
        }
        else if (ancestors > 0) {
            assert(ancestors < 3); // TODO: How "wrong" does the graph has to be to have ancestor>2? should be impossible (at least logically)
            edge_t edge_trio[3];
            vertex_t vertex_trio[3];//
            int child_index = ancestors;
            for (size_t j = 0; j <= ancestors; ++j) {
                edge_trio[j] = *(rng.first + j);
                vertex_trio[j] = target(edge_trio[j], pedigree_graph);
            }
            for (size_t p = 0; p < child_index; ++p) {
                if (vertex_trio[child_index] < vertex_trio[p]) {
                    boost::swap(vertex_trio[child_index], vertex_trio[p]);
                    boost::swap(edge_trio[child_index], edge_trio[p]);
                }
                add_edge(vertex_trio[p], vertex_trio[child_index],
                         {edge_types[edge_trio[p]], lengths[edge_trio[p]] + lengths[edge_trio[child_index]]},
                         pedigree_graph);
            }
            clear_vertex(v, pedigree_graph);

#ifdef DEBUG_RGRAPH
    std::cout << "CleanV2 anc==" << ancestors << " remove: " << v << std::endl;
#endif
        }

    }
    PrintDebugEdges("After SimplifyPedigree", pedigree_graph);

}

void dng::RelationshipGraph::UpdateLabelsNodeIds(dng::Graph &pedigree_graph,
        dng::ReadGroups &rgs, std::vector<size_t> &node_ids) {

    auto labels = get(boost::vertex_label, pedigree_graph);

    labels_.clear();
    labels_.reserve(128);
    std::size_t vid = 0;
    for (std::size_t u = 0; u < num_nodes_; ++u) {
        if (out_degree(u, pedigree_graph) == 0) {
            continue;
        }
        if (!labels[u].empty()) {
            labels_.push_back(labels[u]);
        } else {
            labels_.push_back(DNG_SM_PREFIX "unnamed_node_" + utility::to_pretty(vid));
        }
        node_ids[u] = vid++;
    }
    num_nodes_ = vid;

    // Update rgs so we know what libraries to filter out when building the vcf.
    EraseRemovedLibraries(rgs, node_ids);

    auto update_position = [&node_ids, vid](size_t pos) -> size_t {
        for (; pos < node_ids.size() && node_ids[pos] == -1; ++pos)
            /*noop*/;
        return (pos < node_ids.size()) ? node_ids[pos] : vid;
    };

    first_founder_ = update_position(first_founder_);
    first_nonfounder_ = update_position(first_nonfounder_);
    first_somatic_ = update_position(first_somatic_);
    first_library_ = update_position(first_library_);



#ifdef DEBUG_RGRAPH
    auto sex  = get(boost::vertex_sex, pedigree_graph);
    std::cout << num_nodes_ << std::endl;
    std::cout << "After remove some nodes: new V labels_sizse: " << vid << "\t" << labels_.size() << std::endl;
    for (auto item : labels_) {
        std::cout << item << std::endl;
    }
    std::cout << "node_ids:" << std::endl;
    for (int j = 0; j < node_ids.size(); ++j) {
        std::cout << j << " -> " << node_ids[j] << "\tsex:"
                << (int) sex[j] << std::endl;
    }
#endif
    PrintDebugEdges("END UpdateLabelsNodeIds", pedigree_graph);


}

void dng::RelationshipGraph::EraseRemovedLibraries(dng::ReadGroups &rgs,
        std::vector<size_t> &node_ids) {

    std::vector<std::string> bad_libraries;
    std::size_t num_libraries = rgs.libraries().size();
    bad_libraries.reserve(num_libraries);
    auto it = rgs.libraries().begin();
    for (size_t u = first_library_; u < node_ids.size(); ++u, ++it) {
        if (node_ids[u] != -1) {
            continue;
        }

        bad_libraries.push_back(*it);
    }
    rgs.EraseLibraries(bad_libraries);

}


void dng::RelationshipGraph::CreateFamiliesInfo(dng::Graph &pedigree_graph,
        family_labels_t &family_labels, std::vector<vertex_t> &pivots) {

//    //PR_NOTE:: only need 2 variables here
//    family_labels_t family_labels(num_families);
//    vector<vertex_t> pivots(num_families, dummy_index);

    auto groups = get(boost::vertex_group, pedigree_graph);
    auto families = get(boost::edge_family, pedigree_graph);

    // Calculate the connected components.  This defines independent sections
    // of the graph.
    std::size_t num_groups = connected_components(pedigree_graph, groups);

    // Calculate the biconnected components and articulation points.
    // This defines "nuclear" families and pivot individuals.
    // Nodes which have no edges will not be part of any family.
    std::vector<vertex_t> articulation_vertices;
    std::size_t num_families =
            biconnected_components(pedigree_graph, families,
                                   back_inserter(articulation_vertices)).first;

    family_labels = family_labels_t(num_families);
    pivots = std::vector<vertex_t>(num_families, DUMMY_INDEX);

    // Determine which edges belong to which nuclear families.
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
        family_labels[families[*ei]].push_back(*ei);
    }

    // Determine the last family in each group.  All singleton groups will have
    // a value of -1 since they have no family assignment.
    //XXX: strictly size_t != -1, is it good to overflow trick here?
    typedef std::deque<std::size_t> root_families_t;
    root_families_t root_families(num_groups, -1);
    for (std::size_t f = 0; f < family_labels.size(); ++f) {
        // last one wins
        auto first_edge = family_labels[f][0];
        auto src_vertex = source(first_edge, pedigree_graph);
        root_families[groups[src_vertex]] = f;

    }

    // Identify the pivot for each family.
    // The pivot will be the last art. point that has an edge in
    // the group.  The pivot of the last group doesn't matter.
    for (auto a : articulation_vertices) {
        boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = out_edges(a, pedigree_graph); ei != ei_end; ++ei) {
            // Just overwrite existing value so that the last one wins.
            pivots[families[*ei]] = a;
        }
    }

    // Root Pivots are special
    for (auto f : root_families) {
        if (f != -1) { //Assign to non-singleton groups
            pivots[f] = DUMMY_INDEX;
        }
    }


#ifdef DEBUG_RGRAPH
    for (int l = 0; l < family_labels.size(); ++l) {
        std::cout << "Families : " << l << "\tPivots: " << pivots[l] << "\t";
        for (auto f : family_labels[l]) {
            std::cout << f << " ";
        }
        std::cout << std::endl;
    }
#endif
}

void dng::RelationshipGraph::CreatePeelingOps(
        dng::Graph &pedigree_graph, const std::vector<size_t> &node_ids,
        family_labels_t &family_labels, std::vector<vertex_t> &pivots) {

    auto edge_types = get(boost::edge_type, pedigree_graph);
    auto lengths = get(boost::edge_length, pedigree_graph);
    auto sex  = get(boost::vertex_sex, pedigree_graph);

    ResetFamilyInfo();
    for (std::size_t i = first_founder_; i < first_nonfounder_; ++i) {

//        auto k = std::find(node_ids.begin(), node_ids.end(), i) - node_ids.begin();
        auto it = std::find(node_ids.begin(), node_ids.end(), i);
        assert (it != node_ids.end());
        auto k = std::distance(node_ids.begin(), it);

        transitions_[i] = {TransitionType::Founder, static_cast<size_t>(-1),
                          static_cast<size_t>(-1), 0, 0, sex[k]};
        //PR_NOTE(RC): I'm pretty sure this can be optimized.
        //PR_NOTE(SW): I think the k==i+1 is true all the time
        // But, let's redo this part later, after all inheritance models are implemented.
        // Maybe non of these are needed since we might never use these info.
        // The original code might just work fine.
//        for (std::size_t i = first_founder_; i < first_nonfounder_; ++i) {      +
//            transitions_[i] = {TransitionType::Founder, static_cast<size_t>(-1), static_cast<size_t>(-1), 0, 0};
//        }
        //TODO: faster hack! k=i+1, needs to double check if k==i+1 is always true
    }


//    index = get(vertex_index, pedigree_graph);
    // Detect Family Structure and pivot positions
    for (std::size_t k = 0; k < family_labels.size(); ++k) {
        auto &family_edges = family_labels[k];

#ifdef DEBUG_RGRAPH
        std::cout << "\nStart family: " << k << "\nEdges: ";
        for (auto a : family_edges) {
            std::cout << a << " ";
        }std::cout << "" << std::endl;
#endif
        // Sort edges based on type and target
        boost::sort(family_edges, [&](edge_t x, edge_t y) -> bool {
            return (edge_types(x) < edge_types(y)) && (target(x, pedigree_graph) < target(y, pedigree_graph));
        });

#ifdef DEBUG_RGRAPH
    std::cout << "Sorted edges: ";
    for (auto a : family_edges) {
        std::cout << a << " ";
    }
    std::cout << "" << std::endl;
#endif
        // Find the range of the parent types
        auto pos = boost::find_if(family_edges, [&](edge_t x) -> bool {
            return (edge_types(x) != EdgeType::Spousal);
        });
        size_t num_parent_edges = distance(family_edges.begin(), pos);

#ifdef DEBUG_RGRAPH
    std::cout << "family_edge.size(): " << family_edges.size() << "\tnum_parent_E: " << num_parent_edges << "\t"
            << "\tpos!=(EdgeType::Spousal): " << *pos << "\tpivot_V: "<< pivots[k] << std::endl;
#endif
        // Check to see what type of graph we have
        if (num_parent_edges == 0) {
            // If we do not have a parent-child single branch,
            // we can't construct the pedigree.
            // TODO: Write Error message
            if (family_edges.size() != 1) {
//                return false;
                throw std::runtime_error("Unable to construct peeler for pedigree;  "
                        "do not have a parent-child single branch");
            }
            // Create a mitotic peeling operation.
            auto child_index = target(*pos, pedigree_graph);
            size_t parent = node_ids[source(*pos, pedigree_graph)];
            size_t child = node_ids[child_index];

            TransitionType tt = (edge_types(*pos) == EdgeType::Library) ?
                                TransitionType::Library : TransitionType::Somatic;
            transitions_[child] = {tt, parent, static_cast<size_t>(-1),
                    lengths[*pos], 0, sex[child_index]};//TODO(SW): Double check the index is correct
            family_members_.push_back({parent, child});

#ifdef DEBUG_RGRAPH
    std::cout << "=numParentEdge==0: parent " << parent << "\tChild " << child << std::endl;
    std::cout << "===pivots[k]: " << pivots[k] << "\t\tnode_id:" << node_ids[pivots[k]] << std::endl;
#endif
            if (node_ids[pivots[k]] == child) {
                peeling_ops_.push_back(peel::op::DOWN);
//                std::cout << "=====ADD OP: down" << std::endl;
            } else {
                peeling_ops_.push_back(peel::op::UP);
//                std::cout << "=====ADD OP: up" << std::endl;
                if (pivots[k] == DUMMY_INDEX) {
                    roots_.push_back(parent);
                }
            }

        } else if (num_parent_edges == 1) {
            // If this family contains no children, skip it
            if (pos == family_edges.end()) {
                continue;
            }
            // We have a nuclear family with 1 or more children
            size_t dad = node_ids[source(family_edges.front(), pedigree_graph)];
            size_t mom = node_ids[target(family_edges.front(), pedigree_graph)];

            family_members_.push_back({dad, mom});
            auto &family_members = family_members_.back();

#ifdef DEBUG_RGRAPH
    std::cout << "==numParentEdge==1: dad: " << dad << "\tmom: " << mom << std::endl;
#endif
            while (pos != family_edges.end()) {
                auto child_index = target(*pos, pedigree_graph);
                vertex_t child = node_ids[target(*pos, pedigree_graph)];
                transitions_[child] = {TransitionType::Germline, dad, mom,
                    lengths[*pos], lengths[*(pos + 1)], sex[child_index]}; //TODO(SW): Double check the index is correct
                family_members.push_back(child); // Child

#ifdef DEBUG_RGRAPH
    std::cout << "==Add child to family: " << child << std::endl;
#endif
                // child edges come in pairs
                ++pos;
                assert(node_ids[target(*pos, pedigree_graph)] == child);
                ++pos;
            }
            if (node_ids[pivots[k]] == node_ids[DUMMY_INDEX]) {
                // A family without a pivot is a root family
                peeling_ops_.push_back(peel::op::TOFATHER);
                roots_.push_back(family_members[0]);

#ifdef DEBUG_RGRAPH
    std::cout << "=====ADD OP: toFather Root" << std::endl;
#endif
            } else {
                auto pivot_pos = boost::range::find(family_members, node_ids[pivots[k]]);
                size_t p = distance(family_members.begin(), pivot_pos);

#ifdef DEBUG_RGRAPH
    std::cout << "===pivot_pos: " << *pivot_pos << "," << node_ids[pivots[k]] <<
            "\tdistance_pivot_to_family_member_begin(): " << p << std::endl;
#endif
                if (p == 0) {
                    peeling_ops_.push_back(peel::op::TOFATHER);
//                    std::cout << "=====ADD OP: toFather" << std::endl;
                } else if (p == 1) {
                    peeling_ops_.push_back(peel::op::TOMOTHER);
//                    std::cout << "=====ADD OP: toMother" << std::endl;
                } else if (p == 2) {
                    peeling_ops_.push_back(peel::op::TOCHILD);
//                    std::cout << "=====ADD OP: toChild" << std::endl;
                } else {
                    peeling_ops_.push_back(peel::op::TOCHILD);
//                    std::cout << "=====ADD OP: toChild swap" << std::endl;
                    boost::swap(family_members[p], family_members[2]);
                }
            }

        } else {
            throw std::runtime_error("Unable to construct peeler for pedigree; Not a zero-loop pedigree");
            // TODO: write error message
        }
    }


}

void dng::RelationshipGraph::ResetFamilyInfo(){
    roots_.clear();
    roots_.reserve(16);
    family_members_.clear();
    family_members_.reserve(128);
    peeling_ops_.clear();
    peeling_ops_.reserve(128);

    // Resize the information in the pedigree
    transitions_.resize(num_nodes_);

}

void dng::RelationshipGraph::PrintDebugEdges(const std::string &prefix,
                                             const dng::Graph &pedigree_graph) {

#ifdef DEBUG_RGRAPH

//    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    int verbose_level = 2;
    PropVertexIndex index = get(boost::vertex_index, pedigree_graph);
    auto sex = get(boost::vertex_sex, pedigree_graph);

    boost::graph_traits<Graph>::edge_iterator ei2, ei_end2;

    std::cout << "==DEBUG: " << prefix << ": V: " << num_vertices(pedigree_graph) << "\tE: " << num_edges(pedigree_graph) << std::endl;
    for (tie(ei2, ei_end2) = edges(pedigree_graph); ei2 != ei_end2; ++ei2) {
        std::cout << "(" << index[source(*ei2, pedigree_graph)] << ","
                << index[target(*ei2, pedigree_graph)] << ") ";
    }
    std::cout << std::endl;
    if (verbose_level > 1) {
        boost::graph_traits<Graph>::vertex_iterator vi, vi_end;
        for (boost::tie(vi, vi_end) = vertices(pedigree_graph); vi != vi_end; ++vi) {
            std::cout << "Vertex_Sex:" << *vi << "_" << (int) sex[*vi] << ". ";
        }
        std::cout << "\t\t==" << std::endl;
    }

    std::cout << "Founder, Non_F, Somatic, Lib: " << first_founder_ << "\t"
            << first_nonfounder_ << "\t" << first_somatic_ << "\t"
            << first_library_ << std::endl;
    std::cout << "==END==\n" << std::endl;

#endif

}
